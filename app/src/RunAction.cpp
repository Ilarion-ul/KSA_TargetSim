#include "RunAction.h"

#include <G4Run.hh>
#include <G4SystemOfUnits.hh>
#include <G4Threading.hh>

#if __has_include(<filesystem>)
#include <filesystem>
#define KSA_HAVE_STD_FILESYSTEM 1
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
#define KSA_HAVE_STD_FILESYSTEM 0
#else
#error "No filesystem support available"
#endif

#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <utility>

#ifdef KSA_USE_ROOT
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TTree.h>
#endif

namespace {
constexpr int kPlateHeatmapBinsX = 50;
constexpr int kPlateHeatmapBinsY = 50;
constexpr int kEdepBinsX = 40;
constexpr int kEdepBinsY = 40;
constexpr int kEdepBinsZ = 80;
constexpr int kSurfBinsX = 80;
constexpr int kSurfBinsY = 80;
constexpr int kSurfBinsZ = 120;

size_t DeterminePlateCount(const AppConfig& config) {
  if (config.target.type == "U-Mo" || (config.target.type == "W-Ta" && config.geometry.simpleCylinder)) {
    return config.target.plate_thicknesses_mm.size();
  }
  return 1;
}

RunAction::HeatmapBounds DeterminePlateBounds(const AppConfig& config) {
  double halfXY = 0.0;
  if (config.target.type == "U-Mo" || (config.target.type == "W-Ta" && config.geometry.simpleCylinder)) {
    halfXY = 0.5 * config.target.plate_xy_mm;
  } else {
    halfXY = config.target.radius_mm;
  }
  return {-halfXY, halfXY, -halfXY, halfXY, 0.0, 0.0};
}

RunAction::HeatmapBounds DetermineEdepBounds(const AppConfig& config) {
  const auto plateBounds = DeterminePlateBounds(config);
  double zMin = 0.0;
  double zMax = 0.0;
  if (config.target.type == "U-Mo") {
    const double totalAssemblyLen = config.geometry.total_assembly_len_mm;
    const double beamlineVacLen = config.geometry.beamline_vacuum_len_mm;
    const double targetStartZ = -0.5 * totalAssemblyLen + beamlineVacLen;
    const double targetRegionLen = totalAssemblyLen - beamlineVacLen;
    zMin = targetStartZ;
    zMax = targetStartZ + targetRegionLen;
  } else if (config.target.type == "W-Ta" && config.geometry.simpleCylinder) {
    zMin = -0.5 * config.target.assembly_thickness_mm;
    zMax = 0.5 * config.target.assembly_thickness_mm;
  } else {
    const double totalTargetT = config.target.substrate_thickness_mm + config.target.coating_thickness_mm;
    zMin = -0.5 * totalTargetT;
    zMax = 0.5 * totalTargetT;
  }
  return {plateBounds.xMinMm, plateBounds.xMaxMm, plateBounds.yMinMm, plateBounds.yMaxMm, zMin, zMax};
}


std::vector<double> BuildUMoInterPlateGapsMmForBounds(const AppConfig& config) {
  const auto& t = config.target;
  const size_t n = t.plate_thicknesses_mm.size();
  if (n <= 1) return {};
  if (!t.inter_plate_gaps_mm.empty() && t.inter_plate_gaps_mm.size() == n - 1) {
    return t.inter_plate_gaps_mm;
  }
  std::vector<double> gaps(n - 1, t.gap_rear_mm);
  const size_t split = static_cast<size_t>(std::min<int>(t.gap_split_index, static_cast<int>(n - 1)));
  for (size_t i = 0; i < split; ++i) gaps[i] = t.gap_front_mm;
  if (n == 12 && t.gap_mid_mm > 0.0 && t.gap_front_mm == 3.0 && t.gap_rear_mm == 1.75) {
    std::fill(gaps.begin(), gaps.end(), t.gap_mid_mm);
  }
  return gaps;
}
std::pair<double, double> DeterminePlateStackZBounds(const AppConfig& config) {
  if (config.target.type == "U-Mo") {
    const double totalAssemblyLen = config.geometry.total_assembly_len_mm;
    const double beamlineVacLen = config.geometry.beamline_vacuum_len_mm;
    const double targetRegionLen = totalAssemblyLen - beamlineVacLen;
    const double clearance = config.geometry.target_region_extra_clearance_mm;
    const double targetStartZ = -0.5 * totalAssemblyLen + beamlineVacLen;
    const auto& plateTs = config.target.plate_thicknesses_mm;
    const double cladFront = config.target.clad_thickness_front_mm;
    const double cladRest = config.target.clad_thickness_rest_mm;
    const double gapInOut = config.target.gap_inout_mm;
    const auto interGaps = BuildUMoInterPlateGapsMmForBounds(config);
    double stackLen = gapInOut;
    for (size_t i = 0; i < plateTs.size(); ++i) {
      const double coreT = plateTs[i];
      const double clad = (i < 4) ? cladFront : cladRest;
      stackLen += coreT + 2.0 * clad;
      if (i + 1 < plateTs.size()) {
        stackLen += (i < interGaps.size()) ? interGaps[i] : gapInOut;
      }
    }
    stackLen += gapInOut;
    const double stackStartZ = targetStartZ + clearance + config.target.entrance_window_mm;
    const double stackEndZ = stackStartZ + stackLen;
    const double minZ = std::max(stackStartZ, targetStartZ);
    const double maxZ = std::min(stackEndZ, targetStartZ + targetRegionLen);
    return {minZ, maxZ};
  }
  if (config.target.type == "W-Ta" && config.geometry.simpleCylinder) {
    const auto& wPlateTs = config.target.plate_thicknesses_mm;
    const double waterGap = config.target.water_gap_mm;
    const double ta = config.target.clad_ta_mm;
    const double ti = config.target.buffer_ti_mm;
    double stackT = 8.0 * waterGap;
    for (double tmm : wPlateTs) {
      stackT += tmm + 2.0 * ti + 2.0 * ta;
    }
    return {-0.5 * stackT, 0.5 * stackT};
  }
  const double totalTargetT = config.target.substrate_thickness_mm + config.target.coating_thickness_mm;
  return {-0.5 * totalTargetT, 0.5 * totalTargetT};
}

constexpr int kNeutronFluxGroups = 20;
constexpr double kElementaryChargeC = 1.602176634e-19;

double MeVToJ(double mev) { return mev * 1.602176634e-13; }
double MmToCm(double mmVal) { return mmVal * 0.1; }

enum class MaterialId : int { Unknown = 0, W = 1, Ta = 2, U = 3, Al = 4, Water = 5, Ti = 6 };
enum class LayerId : int { Unknown = 0, Substrate = 1, Coating = 2, Gap = 3, Buffer = 4 };

struct VoxelLabel {
  int material_id{0};
  std::string material_name{"Unknown"};
  int volume_id{0};
  std::string volume_name{"UnknownVolume"};
  int layer_id{0};
  std::string layer_name{"unknown"};
};

struct MeshVoxelRow {
  int mesh_id{0};
  long long voxel_id{0};
  int ix{0};
  int iy{0};
  int iz{0};
  double x_mm{0.0};
  double y_mm{0.0};
  double z_mm{0.0};
  double volume_cm3{0.0};
  VoxelLabel label{};
  double edep_MeV_per_primary{0.0};
  double edep_J_per_primary{0.0};
  double damage_energy_eV_per_primary{0.0};
  int n_events_scored{0};
  double flux_n_total_per_cm2_per_primary{0.0};
  std::vector<double> flux_n_groups_per_cm2_per_primary;
  double h_prod_per_primary{0.0};
  double he_prod_per_primary{0.0};
};

std::vector<double> BuildNeutronEnergyGroupEdgesMeV() {
  std::vector<double> edges;
  edges.reserve(kNeutronFluxGroups + 1);
  const double eMin = 2.5e-9; // 0.0025 eV
  const double eMax = 5.0;
  const double logMin = std::log10(eMin);
  const double logMax = std::log10(eMax);
  for (int i = 0; i <= kNeutronFluxGroups; ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(kNeutronFluxGroups);
    edges.push_back(std::pow(10.0, logMin + t * (logMax - logMin)));
  }
  return edges;
}


std::vector<double> BuildNeutronEnergyGroupEdgesLinearMeV() {
  std::vector<double> edges;
  edges.reserve(kNeutronFluxGroups + 1);
  const double eMin = 2.5e-9; // 0.0025 eV
  const double eMax = 5.0;
  for (int i = 0; i <= kNeutronFluxGroups; ++i) {
    const double t = static_cast<double>(i) / static_cast<double>(kNeutronFluxGroups);
    edges.push_back(eMin + t * (eMax - eMin));
  }
  return edges;
}

VoxelLabel BuildLabel(MaterialId m, LayerId l, int volumeId, const std::string& volumeName) {
  VoxelLabel out;
  out.material_id = static_cast<int>(m);
  out.layer_id = static_cast<int>(l);
  out.volume_id = volumeId;
  out.volume_name = volumeName;
  switch (m) {
    case MaterialId::W: out.material_name = "W"; break;
    case MaterialId::Ta: out.material_name = "Ta"; break;
    case MaterialId::U: out.material_name = "U"; break;
    case MaterialId::Al: out.material_name = "Al"; break;
    case MaterialId::Water: out.material_name = "Water"; break;
    case MaterialId::Ti: out.material_name = "Ti"; break;
    default: out.material_name = "Unknown"; break;
  }
  switch (l) {
    case LayerId::Substrate: out.layer_name = "substrate"; break;
    case LayerId::Coating: out.layer_name = "coating"; break;
    case LayerId::Gap: out.layer_name = "gap"; break;
    case LayerId::Buffer: out.layer_name = "buffer"; break;
    default: out.layer_name = "unknown"; break;
  }
  return out;
}

VoxelLabel ClassifyTargetByZ(const AppConfig& config, double zMm, double zMin, double zMax) {
  if (zMm < zMin || zMm > zMax) {
    return BuildLabel(MaterialId::Unknown, LayerId::Unknown, 0, "OutsideMesh");
  }

  if (config.target.type == "W-Ta" && config.geometry.simpleCylinder) {
    const auto& plateTs = config.target.plate_thicknesses_mm;
    const double gap = config.target.water_gap_mm;
    const double ta = config.target.clad_ta_mm;
    const double ti = config.target.buffer_ti_mm;
    double cursor = zMin;
    auto inSeg = [&](double a, double b) { return zMm >= a && zMm < b; };
    for (size_t i = 0; i < plateTs.size(); ++i) {
      double next = cursor + gap;
      if (inSeg(cursor, next)) return BuildLabel(MaterialId::Water, LayerId::Gap, 300, "TargetWaterGap");
      cursor = next;
      next = cursor + ta;
      if (inSeg(cursor, next)) return BuildLabel(MaterialId::Ta, LayerId::Coating, 201, "TargetCoatingTa");
      cursor = next;
      next = cursor + ti;
      if (inSeg(cursor, next)) return BuildLabel(MaterialId::Ti, LayerId::Buffer, 202, "TargetBufferTi");
      cursor = next;
      next = cursor + plateTs[i];
      if (inSeg(cursor, next)) return BuildLabel(MaterialId::W, LayerId::Substrate, 101, "TargetSubstrateW");
      cursor = next;
      next = cursor + ti;
      if (inSeg(cursor, next)) return BuildLabel(MaterialId::Ti, LayerId::Buffer, 202, "TargetBufferTi");
      cursor = next;
      next = cursor + ta;
      if (inSeg(cursor, next)) return BuildLabel(MaterialId::Ta, LayerId::Coating, 201, "TargetCoatingTa");
      cursor = next;
    }
    if (zMm >= cursor && zMm <= zMax) return BuildLabel(MaterialId::Water, LayerId::Gap, 300, "TargetWaterGap");
  }

  if (config.target.type == "U-Mo") {
    const auto& plateTs = config.target.plate_thicknesses_mm;
    const double cladFront = config.target.clad_thickness_front_mm;
    const double cladRest = config.target.clad_thickness_rest_mm;
    const double gapInOut = config.target.gap_inout_mm;
    const double gapMid = config.target.gap_mid_mm;
    double cursor = zMin;
    auto inSeg = [&](double a, double b) { return zMm >= a && zMm < b; };
    for (size_t i = 0; i < plateTs.size(); ++i) {
      const double gapBefore = (i == 0) ? gapInOut : gapMid;
      double next = cursor + gapBefore;
      if (inSeg(cursor, next)) return BuildLabel(MaterialId::Water, LayerId::Gap, 300, "WaterGap");
      cursor = next;
      const double clad = (i < 4) ? cladFront : cladRest;
      next = cursor + clad;
      if (inSeg(cursor, next)) return BuildLabel(MaterialId::Al, LayerId::Coating, 210, "PlateCladAl");
      cursor = next;
      next = cursor + plateTs[i];
      if (inSeg(cursor, next)) return BuildLabel(MaterialId::U, LayerId::Substrate, 110, "PlateU");
      cursor = next;
      next = cursor + clad;
      if (inSeg(cursor, next)) return BuildLabel(MaterialId::Al, LayerId::Coating, 210, "PlateCladAl");
      cursor = next;
    }
    if (zMm >= cursor && zMm <= zMax) return BuildLabel(MaterialId::Water, LayerId::Gap, 300, "WaterGap");
  }

  if (config.target.type == "U-Al") {
    const double totalTargetT = config.target.substrate_thickness_mm + config.target.coating_thickness_mm;
    const double z0 = -0.5 * totalTargetT;
    const double z1 = z0 + config.target.substrate_thickness_mm;
    if (zMm < z1) return BuildLabel(MaterialId::U, LayerId::Substrate, 110, "TargetSubstrate");
    return BuildLabel(MaterialId::Al, LayerId::Coating, 210, "TargetCoating");
  }

  return BuildLabel(MaterialId::Unknown, LayerId::Unknown, 0, "UnknownVolume");
}

std::vector<MeshVoxelRow> BuildMeshRows(const AppConfig& config,
                                        const RunAction::HeatmapBounds& bounds,
                                        int binsX,
                                        int binsY,
                                        int binsZ,
                                        int nEvents,
                                        const std::vector<double>& edep3dMeV) {
  std::vector<MeshVoxelRow> rows;
  if (binsX <= 0 || binsY <= 0 || binsZ <= 0) return rows;
  rows.reserve(static_cast<size_t>(binsX) * static_cast<size_t>(binsY) * static_cast<size_t>(binsZ));

  const auto zBounds = DeterminePlateStackZBounds(config);
  const double zMin = zBounds.first;
  const double zMax = zBounds.second;
  const double xMin = bounds.xMinMm;
  const double xMax = bounds.xMaxMm;
  const double yMin = bounds.yMinMm;
  const double yMax = bounds.yMaxMm;
  const double zOrigMin = bounds.zMinMm;
  const double zOrigMax = bounds.zMaxMm;
  const double zOrigSpan = zOrigMax - zOrigMin;
  const double zSpan = zMax - zMin;

  const double dx_mm = (xMax - xMin) / static_cast<double>(binsX);
  const double dy_mm = (yMax - yMin) / static_cast<double>(binsY);
  const double dz_mm = zSpan / static_cast<double>(binsZ);
  const double voxel_cm3 = MmToCm(dx_mm) * MmToCm(dy_mm) * MmToCm(dz_mm);

  std::vector<double> remapped(static_cast<size_t>(binsX) * static_cast<size_t>(binsY) * static_cast<size_t>(binsZ), 0.0);
  if (zOrigSpan > 0.0 && zSpan > 0.0) {
    for (int iz = 0; iz < binsZ; ++iz) {
      const double zCenterOrig = zOrigMin + (static_cast<double>(iz) + 0.5) * zOrigSpan / static_cast<double>(binsZ);
      if (zCenterOrig < zMin || zCenterOrig > zMax) continue;
      int izPlate = static_cast<int>(std::floor((zCenterOrig - zMin) / zSpan * binsZ));
      if (izPlate < 0) izPlate = 0;
      if (izPlate >= binsZ) izPlate = binsZ - 1;
      for (int iy = 0; iy < binsY; ++iy) {
        for (int ix = 0; ix < binsX; ++ix) {
          const size_t src = (static_cast<size_t>(iz) * static_cast<size_t>(binsY) + static_cast<size_t>(iy)) *
                                 static_cast<size_t>(binsX) +
                             static_cast<size_t>(ix);
          const size_t dst = (static_cast<size_t>(izPlate) * static_cast<size_t>(binsY) + static_cast<size_t>(iy)) *
                                 static_cast<size_t>(binsX) +
                             static_cast<size_t>(ix);
          if (src < edep3dMeV.size() && dst < remapped.size()) remapped[dst] += edep3dMeV[src];
        }
      }
    }
  }

  long long voxelId = 0;
  for (int iz = 0; iz < binsZ; ++iz) {
    const double zCenter = zMin + (static_cast<double>(iz) + 0.5) * dz_mm;
    const auto label = ClassifyTargetByZ(config, zCenter, zMin, zMax);
    for (int iy = 0; iy < binsY; ++iy) {
      const double yCenter = yMin + (static_cast<double>(iy) + 0.5) * dy_mm;
      for (int ix = 0; ix < binsX; ++ix) {
        const double xCenter = xMin + (static_cast<double>(ix) + 0.5) * dx_mm;
        const size_t idx = (static_cast<size_t>(iz) * static_cast<size_t>(binsY) + static_cast<size_t>(iy)) *
                               static_cast<size_t>(binsX) +
                           static_cast<size_t>(ix);
        MeshVoxelRow row;
        row.mesh_id = 1;
        row.voxel_id = voxelId++;
        row.ix = ix;
        row.iy = iy;
        row.iz = iz;
        row.x_mm = xCenter;
        row.y_mm = yCenter;
        row.z_mm = zCenter;
        row.volume_cm3 = voxel_cm3;
        row.label = label;
        row.edep_MeV_per_primary = (idx < remapped.size()) ? remapped[idx] : 0.0;
        row.edep_J_per_primary = MeVToJ(row.edep_MeV_per_primary);
        row.damage_energy_eV_per_primary = row.edep_MeV_per_primary * 1e6; // proxy fallback
        row.n_events_scored = nEvents;
        row.flux_n_total_per_cm2_per_primary = 0.0; // placeholder for dedicated flux scorer
        row.flux_n_groups_per_cm2_per_primary.assign(kNeutronFluxGroups, 0.0);
        row.h_prod_per_primary = 0.0; // placeholder for per-voxel scorer
        row.he_prod_per_primary = 0.0; // placeholder for per-voxel scorer
        rows.push_back(std::move(row));
      }
    }
  }
  return rows;
}

} // namespace

RunAction::RunAction(const AppConfig& config)
    : config_(config),
      plateCount_(DeterminePlateCount(config)),
      plateHeatmapBinsX_(kPlateHeatmapBinsX),
      plateHeatmapBinsY_(kPlateHeatmapBinsY),
      edepBinsX_(kEdepBinsX),
      edepBinsY_(kEdepBinsY),
      edepBinsZ_(kEdepBinsZ),
      plateHeatmapBounds_(DeterminePlateBounds(config)),
      edepBounds_(DetermineEdepBounds(config)) {}
RunAction::~RunAction() = default;

namespace {
#if KSA_HAVE_STD_FILESYSTEM
namespace fs = std::filesystem;
#else
namespace fs = std::experimental::filesystem;
#endif

// Shared accumulators across worker/master RunAction instances.
// In MT Geant4 creates separate action objects per thread, so member counters
// are not automatically merged.
std::mutex gRunSummaryMutex;
double gTotalEdepSubstrate = 0.0;
double gTotalEdepCoating = 0.0;
long long gTotalGamma = 0;
long long gTotalNeutron = 0;
long long gTotalNeutronExit = 0;
long long gTotalNeutronModelExit = 0;
std::vector<double> gTotalPlateEdep;
std::vector<double> gTotalPlateNeutronTrackLen;
std::vector<double> gTotalPlateNeutronHeatmap;
std::vector<double> gTotalEdep3d;
std::vector<RunAction::NeutronSurfaceHit> gNeutronSurfaceHits;
std::vector<RunAction::PhotonSurfaceHit> gPhotonSurfaceHits;
std::vector<double> gTotalPlateNiel;
std::vector<double> gTotalPlateGasH;
std::vector<double> gTotalPlateGasHe;
}

void RunAction::BeginOfRunAction(const G4Run*) {
  // Reset global counters once at run start by master only.
  if (G4Threading::IsMasterThread()) {
    std::lock_guard<std::mutex> lock(gRunSummaryMutex);
    gTotalEdepSubstrate = 0.0;
    gTotalEdepCoating = 0.0;
    gTotalGamma = 0;
    gTotalNeutron = 0;
    gTotalNeutronExit = 0;
    gTotalNeutronModelExit = 0;
    gTotalPlateEdep.assign(plateCount_, 0.0);
    gTotalPlateNeutronTrackLen.assign(plateCount_, 0.0);
    gTotalPlateNeutronHeatmap.assign(static_cast<size_t>(plateHeatmapBinsX_) * static_cast<size_t>(plateHeatmapBinsY_) *
                                         plateCount_,
                                     0.0);
    gTotalEdep3d.assign(static_cast<size_t>(edepBinsX_) * static_cast<size_t>(edepBinsY_) * static_cast<size_t>(edepBinsZ_),
                        0.0);
    gNeutronSurfaceHits.clear();
    gPhotonSurfaceHits.clear();
    gTotalPlateNiel.assign(plateCount_, 0.0);
    gTotalPlateGasH.assign(plateCount_, 0.0);
    gTotalPlateGasHe.assign(plateCount_, 0.0);
  }

  const auto base = fs::path(config_.run.outputDir.empty() ? "results" : config_.run.outputDir);
  fs::create_directories(base / "logs");
  fs::create_directories(base / "root");
  fs::create_directories(base / "vis");

#ifdef KSA_USE_ROOT
  // ROOT global state is not thread-safe for opening/writing one file from
  // all workers; keep file creation/writing on master only.
  if (G4Threading::IsMasterThread()) {
    const auto rootPath = base / "root" / config_.run.outputRootFile;
    rootFile_ = TFile::Open(rootPath.string().c_str(), "RECREATE");
    runTree_ = new TTree("run_summary", "KSA run summary");
  }
#endif
}

void RunAction::EndOfRunAction(const G4Run* run) {
  // Master writes final summary once. Workers only contribute counters.
  if (!G4Threading::IsMasterThread()) {
    return;
  }

  const auto base = fs::path(config_.run.outputDir.empty() ? "results" : config_.run.outputDir);

  double edepSub = 0.0;
  double edepCoat = 0.0;
  long long nGammaTotal = 0;
  long long nNeutronTotal = 0;
  long long nNeutronExitTotal = 0;
  long long nNeutronModelExitTotal = 0;
  std::vector<double> plateEdep;
  std::vector<double> plateNeutronTrackLen;
  std::vector<double> plateNeutronHeatmap;
  std::vector<double> edep3d;
  std::vector<NeutronSurfaceHit> neutronSurfaceHits;
  std::vector<PhotonSurfaceHit> photonSurfaceHits;
  std::vector<double> plateNiel;
  std::vector<double> plateGasH;
  std::vector<double> plateGasHe;
  {
    std::lock_guard<std::mutex> lock(gRunSummaryMutex);
    edepSub = gTotalEdepSubstrate;
    edepCoat = gTotalEdepCoating;
    nGammaTotal = gTotalGamma;
    nNeutronTotal = gTotalNeutron;
    nNeutronExitTotal = gTotalNeutronExit;
    nNeutronModelExitTotal = gTotalNeutronModelExit;
    plateEdep = gTotalPlateEdep;
    plateNeutronTrackLen = gTotalPlateNeutronTrackLen;
    plateNeutronHeatmap = gTotalPlateNeutronHeatmap;
    edep3d = gTotalEdep3d;
    neutronSurfaceHits = gNeutronSurfaceHits;
    photonSurfaceHits = gPhotonSurfaceHits;
    plateNiel = gTotalPlateNiel;
    plateGasH = gTotalPlateGasH;
    plateGasHe = gTotalPlateGasHe;
  }

  std::vector<double> plateEdepMeV;
  plateEdepMeV.reserve(plateEdep.size());
  for (double value : plateEdep) {
    plateEdepMeV.push_back(value / MeV);
  }
  std::vector<double> plateNeutronTrackLenMM;
  plateNeutronTrackLenMM.reserve(plateNeutronTrackLen.size());
  for (double value : plateNeutronTrackLen) {
    plateNeutronTrackLenMM.push_back(value / mm);
  }
  std::vector<double> plateNeutronHeatmapMM;
  plateNeutronHeatmapMM.reserve(plateNeutronHeatmap.size());
  for (double value : plateNeutronHeatmap) {
    plateNeutronHeatmapMM.push_back(value / mm);
  }
  std::vector<double> edep3dMeV;
  edep3dMeV.reserve(edep3d.size());
  for (double value : edep3d) {
    edep3dMeV.push_back(value / MeV);
  }
  std::vector<double> plateNielMeV;
  plateNielMeV.reserve(plateNiel.size());
  for (double value : plateNiel) {
    plateNielMeV.push_back(value / MeV);
  }

  const int nEventsForNorm = run ? run->GetNumberOfEvent() : config_.run.nEvents;
  const double beamCurrentA = config_.beam.I_avg_A;
  const double irradiationTimeS = config_.run.irradiation_time_s;
  const double electronsPerSecond = beamCurrentA / kElementaryChargeC;
  const double electronsTotal = electronsPerSecond * irradiationTimeS;
  const auto neutronGroupEdgesMeV = BuildNeutronEnergyGroupEdgesMeV();
  const auto neutronGroupEdgesLinearMeV = BuildNeutronEnergyGroupEdgesLinearMeV();
  std::vector<MeshVoxelRow> meshRows;
  if (config_.run.enableSwellingOutput) {
    meshRows = BuildMeshRows(config_, edepBounds_, edepBinsX_, edepBinsY_, edepBinsZ_, nEventsForNorm, edep3dMeV);
  }

#ifdef KSA_USE_ROOT
  if (rootFile_ && runTree_) {
    double edep_substrate_MeV = edepSub / MeV;
    double edep_coating_MeV = edepCoat / MeV;
    long long nGamma = nGammaTotal;
    long long nGammaAbove5MeV = nGammaTotal;
    long long nNeutron = nNeutronTotal;
    int nEvents = nEventsForNorm;
    std::string physicsListName = config_.physics.physicsListName;
    int nNeutronExit = static_cast<int>(nNeutronExitTotal);
    int nNeutronModelExit = static_cast<int>(nNeutronModelExitTotal);
    std::vector<double> plate_edep_MeV = plateEdepMeV;
    std::vector<double> plate_neutron_track_len_mm = plateNeutronTrackLenMM;
    double beam_E0_MeV = config_.beam.energy_MeV;
    double beam_Esigma_rel = config_.beam.energy_sigma_rel_1sigma;
    std::string beam_spread_model = config_.beam.energy_spread_model;
    std::string beam_mode = config_.beam.mode;
    double beam_pulse_width_us = config_.beam.pulse_width_us;
    double beam_rep_rate_Hz = config_.beam.rep_rate_Hz;
    double beam_I_pulse_A = config_.beam.I_pulse_A;
    double beam_I_avg_A_input = config_.beam.I_avg_A;
    double beam_duty = (config_.beam.pulse_width_us * 1e-6) * config_.beam.rep_rate_Hz;
    double beam_I_avg_A_calc = beam_I_pulse_A * beam_duty;
    double beam_power_kW = config_.beam.beam_power_kW;

    runTree_->Branch("edep_substrate", &edep_substrate_MeV);
    runTree_->Branch("edep_coating", &edep_coating_MeV);
    runTree_->Branch("nGamma", &nGamma);
    runTree_->Branch("nGammaAbove5MeV", &nGammaAbove5MeV);
    runTree_->Branch("nNeutron", &nNeutron);
    runTree_->Branch("nNeutronExit", &nNeutronExit);
    runTree_->Branch("nNeutronModelExit", &nNeutronModelExit);
    runTree_->Branch("nEvents", &nEvents);
    runTree_->Branch("physicsListName", &physicsListName);
    runTree_->Branch("plate_edep_MeV", &plate_edep_MeV);
    runTree_->Branch("plate_neutron_track_len_mm", &plate_neutron_track_len_mm);
    runTree_->Branch("beam_E0_MeV", &beam_E0_MeV);
    runTree_->Branch("beam_Esigma_rel", &beam_Esigma_rel);
    runTree_->Branch("beam_spread_model", &beam_spread_model);
    runTree_->Branch("beam_mode", &beam_mode);
    runTree_->Branch("beam_pulse_width_us", &beam_pulse_width_us);
    runTree_->Branch("beam_rep_rate_Hz", &beam_rep_rate_Hz);
    runTree_->Branch("beam_I_pulse_A", &beam_I_pulse_A);
    runTree_->Branch("beam_I_avg_A_input", &beam_I_avg_A_input);
    runTree_->Branch("beam_duty", &beam_duty);
    runTree_->Branch("beam_I_avg_A_calc", &beam_I_avg_A_calc);
    runTree_->Branch("beam_power_kW", &beam_power_kW);
    runTree_->Fill();

    rootFile_->cd();
    {
      auto* metaTree = new TTree("RunMeta", "Run metadata");
      std::string target_type = config_.target.type;
      std::string geometry_json;
      {
        std::ostringstream os;
        os << "{";
        os << "\"plate_thicknesses_mm\":[";
        for (size_t i = 0; i < config_.target.plate_thicknesses_mm.size(); ++i) {
          os << config_.target.plate_thicknesses_mm[i];
          if (i + 1 < config_.target.plate_thicknesses_mm.size()) os << ",";
        }
        os << "],";
        os << "\"plate_xy_mm\":" << config_.target.plate_xy_mm << ",";
        os << "\"water_gap_mm\":" << config_.target.water_gap_mm << ",";
        os << "\"clad_ta_mm\":" << config_.target.clad_ta_mm << ",";
        os << "\"buffer_ti_mm\":" << config_.target.buffer_ti_mm << ",";
        os << "\"assembly_thickness_mm\":" << config_.target.assembly_thickness_mm << ",";
        os << "\"clad_thickness_front_mm\":" << config_.target.clad_thickness_front_mm << ",";
        os << "\"clad_thickness_rest_mm\":" << config_.target.clad_thickness_rest_mm << ",";
        os << "\"gap_inout_mm\":" << config_.target.gap_inout_mm << ",";
        os << "\"gap_mid_mm\":" << config_.target.gap_mid_mm << ",";
        os << "\"housing_inner_xy_mm\":" << config_.target.housing_inner_xy_mm << ",";
        os << "\"housing_wall_mm\":" << config_.target.housing_wall_mm << ",";
        os << "\"entrance_window_mm\":" << config_.target.entrance_window_mm << ",";
        os << "\"helium_chamber_len_mm\":" << config_.target.helium_chamber_len_mm << ",";
        os << "\"u7mo_density_g_cm3\":" << config_.target.u7mo_density_g_cm3 << ",";
        os << "\"fill_medium_in_target\":\"" << config_.target.fill_medium_in_target << "\"";
        os << "}";
        geometry_json = os.str();
      }
      double beam_energy_MeV = config_.beam.energy_MeV;
      double beam_energy_sigma_rel = config_.beam.energy_sigma_rel_1sigma;
      double beam_sigma_x_mm = config_.beam.sigmaX_mm;
      double beam_sigma_y_mm = config_.beam.sigmaY_mm;
      double beam_sigma_theta_x_mrad = config_.beam.sigma_theta_x_mrad;
      double beam_sigma_theta_y_mrad = config_.beam.sigma_theta_y_mrad;
      std::string pulse_mode = config_.beam.mode;
      double pulse_width_us = config_.beam.pulse_width_us;
      double rep_rate_Hz = config_.beam.rep_rate_Hz;
      double I_pulse_A = config_.beam.I_pulse_A;
      double I_avg_A = config_.beam.I_avg_A;
      double duty = (config_.beam.pulse_width_us * 1e-6) * config_.beam.rep_rate_Hz;
      double P_avg_kW = config_.beam.beam_power_kW;
      int nEvents = nEventsForNorm;
      int nThreads = config_.run.nThreads;
      double per_primary = 1.0;
      double N_e_per_s = electronsPerSecond;
      metaTree->Branch("target_type", &target_type);
      metaTree->Branch("geometry_json", &geometry_json);
      metaTree->Branch("beam_energy_MeV", &beam_energy_MeV);
      metaTree->Branch("beam_energy_sigma_rel", &beam_energy_sigma_rel);
      metaTree->Branch("beam_sigma_x_mm", &beam_sigma_x_mm);
      metaTree->Branch("beam_sigma_y_mm", &beam_sigma_y_mm);
      metaTree->Branch("beam_sigma_theta_x_mrad", &beam_sigma_theta_x_mrad);
      metaTree->Branch("beam_sigma_theta_y_mrad", &beam_sigma_theta_y_mrad);
      metaTree->Branch("pulse_mode", &pulse_mode);
      metaTree->Branch("pulse_width_us", &pulse_width_us);
      metaTree->Branch("rep_rate_Hz", &rep_rate_Hz);
      metaTree->Branch("I_pulse_A", &I_pulse_A);
      metaTree->Branch("I_avg_A", &I_avg_A);
      metaTree->Branch("duty", &duty);
      metaTree->Branch("P_avg_kW", &P_avg_kW);
      metaTree->Branch("nEvents", &nEvents);
      std::string normalization_mode = "per_primary";
      double beam_current_A = beamCurrentA;
      double irradiation_time_s = irradiationTimeS;
      double electrons_total = electronsTotal;
      std::string units_json = "{\"edep_MeV_per_primary\":\"MeV\",\"edep_J_per_primary\":\"J\",\"damage_energy_eV_per_primary\":\"eV\",\"flux_n_g\":\"1/cm2 per-primary\"}";
      std::ostringstream binsOs;
      binsOs << "[";
      for (size_t i = 0; i < neutronGroupEdgesMeV.size(); ++i) {
        binsOs << neutronGroupEdgesMeV[i];
        if (i + 1 < neutronGroupEdgesMeV.size()) binsOs << ",";
      }
      binsOs << "]";
      std::string neutron_energy_group_edges_MeV = binsOs.str();
      std::ostringstream binsLinearOs;
      binsLinearOs << "[";
      for (size_t i = 0; i < neutronGroupEdgesLinearMeV.size(); ++i) {
        binsLinearOs << neutronGroupEdgesLinearMeV[i];
        if (i + 1 < neutronGroupEdgesLinearMeV.size()) binsLinearOs << ",";
      }
      binsLinearOs << "]";
      std::string neutron_energy_group_edges_linear_MeV = binsLinearOs.str();
      std::string damage_source_status = "derived_from=damage_proxy";
      metaTree->Branch("nThreads", &nThreads);
      metaTree->Branch("per_primary", &per_primary);
      metaTree->Branch("N_e_per_s", &N_e_per_s);
      metaTree->Branch("normalization_mode", &normalization_mode);
      metaTree->Branch("beam_current_A", &beam_current_A);
      metaTree->Branch("irradiation_time_s", &irradiation_time_s);
      metaTree->Branch("electrons_total", &electrons_total);
      metaTree->Branch("units_json", &units_json);
      metaTree->Branch("neutron_energy_group_edges_MeV", &neutron_energy_group_edges_MeV);
      metaTree->Branch("neutron_energy_group_edges_linear_MeV", &neutron_energy_group_edges_linear_MeV);
      metaTree->Branch("damage_source_status", &damage_source_status);
      metaTree->Fill();
      metaTree->Write();
    }
    if (plateHeatmapBinsX_ > 0 && plateHeatmapBinsY_ > 0) {
      const double xMin = plateHeatmapBounds_.xMinMm;
      const double xMax = plateHeatmapBounds_.xMaxMm;
      const double yMin = plateHeatmapBounds_.yMinMm;
      const double yMax = plateHeatmapBounds_.yMaxMm;
      for (size_t plateIdx = 0; plateIdx < plateCount_; ++plateIdx) {
        const std::string name = "plate_neutron_heatmap_" + std::to_string(plateIdx + 1);
        const std::string title = "Plate neutron heatmap " + std::to_string(plateIdx + 1);
        auto* h2 = new TH2D(name.c_str(), title.c_str(), plateHeatmapBinsX_, xMin, xMax, plateHeatmapBinsY_, yMin, yMax);
        const size_t plateOffset = plateIdx * static_cast<size_t>(plateHeatmapBinsX_) *
                                   static_cast<size_t>(plateHeatmapBinsY_);
        for (int iy = 0; iy < plateHeatmapBinsY_; ++iy) {
          for (int ix = 0; ix < plateHeatmapBinsX_; ++ix) {
            const size_t idx = plateOffset + static_cast<size_t>(iy) * static_cast<size_t>(plateHeatmapBinsX_) +
                               static_cast<size_t>(ix);
            if (idx < plateNeutronHeatmapMM.size()) {
              h2->SetBinContent(ix + 1, iy + 1, plateNeutronHeatmapMM[idx]);
            }
          }
        }
        h2->Write();
      }
    }
    if (edepBinsX_ > 0 && edepBinsY_ > 0 && edepBinsZ_ > 0) {
      const double xMin = edepBounds_.xMinMm;
      const double xMax = edepBounds_.xMaxMm;
      const double yMin = edepBounds_.yMinMm;
      const double yMax = edepBounds_.yMaxMm;
      const auto plateStackBounds = DeterminePlateStackZBounds(config_);
      const double zMin = plateStackBounds.first;
      const double zMax = plateStackBounds.second;
      const double zOrigMin = edepBounds_.zMinMm;
      const double zOrigMax = edepBounds_.zMaxMm;
      const double zOrigSpan = zOrigMax - zOrigMin;
      const double zPlateSpan = zMax - zMin;
      auto* h3 = new TH3D("edep_3d", "Energy deposition 3D", edepBinsX_, xMin, xMax, edepBinsY_, yMin, yMax,
                          edepBinsZ_, zMin, zMax);
      auto* h2_edep_xy_mid = new TH2D("h2_edep_xy_mid", "Energy deposition XY (mid-plane)",
                                      edepBinsX_, xMin, xMax, edepBinsY_, yMin, yMax);
      const int midZ = edepBinsZ_ / 2;
      for (int iz = 0; iz < edepBinsZ_; ++iz) {
        const double zCenterOrig = zOrigMin + (static_cast<double>(iz) + 0.5) * zOrigSpan / edepBinsZ_;
        if (zCenterOrig < zMin || zCenterOrig > zMax || zPlateSpan <= 0.0) {
          continue;
        }
        int izPlate = static_cast<int>(std::floor((zCenterOrig - zMin) / zPlateSpan * edepBinsZ_));
        if (izPlate < 0) izPlate = 0;
        if (izPlate >= edepBinsZ_) izPlate = edepBinsZ_ - 1;
        for (int iy = 0; iy < edepBinsY_; ++iy) {
          for (int ix = 0; ix < edepBinsX_; ++ix) {
            const size_t idx = (static_cast<size_t>(iz) * static_cast<size_t>(edepBinsY_) + static_cast<size_t>(iy)) *
                                   static_cast<size_t>(edepBinsX_) +
                               static_cast<size_t>(ix);
            if (idx < edep3dMeV.size()) {
              const double old = h3->GetBinContent(ix + 1, iy + 1, izPlate + 1);
              h3->SetBinContent(ix + 1, iy + 1, izPlate + 1, old + edep3dMeV[idx]);
              if (izPlate == midZ) {
                h2_edep_xy_mid->SetBinContent(ix + 1, iy + 1, edep3dMeV[idx]);
              }
            }
          }
        }
      }
      h3->Write();
      h2_edep_xy_mid->Write();
    }

    if (config_.run.enableSwellingOutput && !meshRows.empty()) {
      auto* meshTree = new TTree("MeshData", "Voxelized mesh data for swelling postprocessing");
      int mesh_id = 0;
      long long voxel_id = 0;
      int ix = 0;
      int iy = 0;
      int iz = 0;
      double x_mm = 0.0;
      double y_mm = 0.0;
      double z_mm = 0.0;
      double volume_cm3 = 0.0;
      int material_id = 0;
      std::string material_name;
      int volume_id = 0;
      std::string volume_name;
      int layer_id = 0;
      std::string layer_name;
      double edep_MeV_per_primary = 0.0;
      double edep_J_per_primary = 0.0;
      double damage_energy_eV_per_primary = 0.0;
      int n_events_scored = 0;
      double flux_n_total_per_cm2_per_primary = 0.0;
      std::vector<double> flux_n;
      double h_prod_per_primary = 0.0;
      double he_prod_per_primary = 0.0;

      meshTree->Branch("mesh_id", &mesh_id);
      meshTree->Branch("voxel_id", &voxel_id);
      meshTree->Branch("ix", &ix);
      meshTree->Branch("iy", &iy);
      meshTree->Branch("iz", &iz);
      meshTree->Branch("x_mm", &x_mm);
      meshTree->Branch("y_mm", &y_mm);
      meshTree->Branch("z_mm", &z_mm);
      meshTree->Branch("volume_cm3", &volume_cm3);
      meshTree->Branch("material_id", &material_id);
      meshTree->Branch("material_name", &material_name);
      meshTree->Branch("volume_id", &volume_id);
      meshTree->Branch("volume_name", &volume_name);
      meshTree->Branch("layer_id", &layer_id);
      meshTree->Branch("layer_name", &layer_name);
      meshTree->Branch("edep_MeV_per_primary", &edep_MeV_per_primary);
      meshTree->Branch("edep_J_per_primary", &edep_J_per_primary);
      meshTree->Branch("damage_energy_eV_per_primary", &damage_energy_eV_per_primary);
      meshTree->Branch("n_events_scored", &n_events_scored);
      meshTree->Branch("flux_n_total_per_cm2_per_primary", &flux_n_total_per_cm2_per_primary);
      meshTree->Branch("flux_n_groups_per_cm2_per_primary", &flux_n);
      meshTree->Branch("h_prod_per_primary", &h_prod_per_primary);
      meshTree->Branch("he_prod_per_primary", &he_prod_per_primary);

      for (const auto& row : meshRows) {
        mesh_id = row.mesh_id;
        voxel_id = row.voxel_id;
        ix = row.ix;
        iy = row.iy;
        iz = row.iz;
        x_mm = row.x_mm;
        y_mm = row.y_mm;
        z_mm = row.z_mm;
        volume_cm3 = row.volume_cm3;
        material_id = row.label.material_id;
        material_name = row.label.material_name;
        volume_id = row.label.volume_id;
        volume_name = row.label.volume_name;
        layer_id = row.label.layer_id;
        layer_name = row.label.layer_name;
        edep_MeV_per_primary = row.edep_MeV_per_primary;
        edep_J_per_primary = row.edep_J_per_primary;
        damage_energy_eV_per_primary = row.damage_energy_eV_per_primary;
        n_events_scored = row.n_events_scored;
        flux_n_total_per_cm2_per_primary = row.flux_n_total_per_cm2_per_primary;
        flux_n = row.flux_n_groups_per_cm2_per_primary;
        h_prod_per_primary = row.h_prod_per_primary;
        he_prod_per_primary = row.he_prod_per_primary;
        meshTree->Fill();
      }
      meshTree->Write();
    }
    if (!plateNielMeV.empty()) {
      auto* h1_niel_plate = new TH1D("h1_niel_plate", "NIEL proxy per plate", static_cast<int>(plateNielMeV.size()), 0.5,
                                     static_cast<double>(plateNielMeV.size()) + 0.5);
      for (size_t i = 0; i < plateNielMeV.size(); ++i) {
        h1_niel_plate->SetBinContent(static_cast<int>(i) + 1, plateNielMeV[i]);
      }
      h1_niel_plate->Write();
    }
    if (!plateGasH.empty()) {
      auto* h1_gas_h_plate =
          new TH1D("h1_gas_h_plate", "Hydrogen production per plate", static_cast<int>(plateGasH.size()), 0.5,
                   static_cast<double>(plateGasH.size()) + 0.5);
      for (size_t i = 0; i < plateGasH.size(); ++i) {
        h1_gas_h_plate->SetBinContent(static_cast<int>(i) + 1, plateGasH[i]);
      }
      h1_gas_h_plate->Write();
    }
    if (!plateGasHe.empty()) {
      auto* h1_gas_he_plate =
          new TH1D("h1_gas_he_plate", "Helium production per plate", static_cast<int>(plateGasHe.size()), 0.5,
                   static_cast<double>(plateGasHe.size()) + 0.5);
      for (size_t i = 0; i < plateGasHe.size(); ++i) {
        h1_gas_he_plate->SetBinContent(static_cast<int>(i) + 1, plateGasHe[i]);
      }
      h1_gas_he_plate->Write();
    }
    if (!neutronSurfaceHits.empty()) {
      auto* neutronTree = new TTree("NeutronSurf", "Neutron surface crossings");
      int event_id = 0;
      double En_MeV = 0.0;
      double x_mm = 0.0;
      double y_mm = 0.0;
      double z_mm = 0.0;
      double cosTheta = 0.0;
      double weight = 1.0;
      double time_ns = 0.0;
      int surface_id = 0;
      std::string surface_name;
      neutronTree->Branch("event_id", &event_id);
      neutronTree->Branch("En_MeV", &En_MeV);
      neutronTree->Branch("x_mm", &x_mm);
      neutronTree->Branch("y_mm", &y_mm);
      neutronTree->Branch("z_mm", &z_mm);
      neutronTree->Branch("cosTheta", &cosTheta);
      neutronTree->Branch("weight", &weight);
      neutronTree->Branch("time_ns", &time_ns);
      neutronTree->Branch("surface_id", &surface_id);
      neutronTree->Branch("surface_name", &surface_name);

      const double xMin = edepBounds_.xMinMm;
      const double xMax = edepBounds_.xMaxMm;
      const double yMin = edepBounds_.yMinMm;
      const double yMax = edepBounds_.yMaxMm;
      const auto plateStackBounds = DeterminePlateStackZBounds(config_);
      const double zMin = plateStackBounds.first;
      const double zMax = plateStackBounds.second;
      auto* h2_exit_xy_downstream = new TH2D("h2_neutron_exit_xy_downstream", "Neutron exit (downstream)",
                                             kSurfBinsX, xMin, xMax, kSurfBinsY, yMin, yMax);
      auto* h2_exit_xy_upstream =
          new TH2D("h2_neutron_exit_xy_upstream", "Neutron exit (upstream)", kSurfBinsX, xMin, xMax, kSurfBinsY, yMin, yMax);
      auto* h2_exit_yz_side_x =
          new TH2D("h2_neutron_exit_yz_side_x", "Neutron exit (side x)", kSurfBinsY, yMin, yMax, kSurfBinsZ, zMin, zMax);
      auto* h2_exit_xz_side_y =
          new TH2D("h2_neutron_exit_xz_side_y", "Neutron exit (side y)", kSurfBinsX, xMin, xMax, kSurfBinsZ, zMin, zMax);
      const double xRange = xMax - xMin;
      const double yRange = yMax - yMin;
      const double sMax = 2.0 * yRange + 2.0 * xRange;
      auto* h2_exit_side_surface =
          new TH2D("h2_neutron_exit_side_surface", "Neutron exit (side surfaces)", kSurfBinsX, 0.0, sMax, kSurfBinsZ, zMin, zMax);
      for (const auto& hit : neutronSurfaceHits) {
        event_id = hit.event_id;
        En_MeV = hit.En_MeV;
        x_mm = hit.x_mm;
        y_mm = hit.y_mm;
        z_mm = hit.z_mm;
        cosTheta = hit.cosTheta;
        weight = hit.weight;
        time_ns = hit.time_ns;
        surface_id = hit.surface_id;
        switch (surface_id) {
          case 0:
            surface_name = "downstream";
            break;
          case 1:
            surface_name = "upstream";
            break;
          case 2:
            surface_name = "side_x";
            break;
          case 3:
            surface_name = "side_y";
            break;
          default:
            surface_name = "unknown";
            break;
        }
        if (surface_id == 0) {
          h2_exit_xy_downstream->Fill(x_mm, y_mm, weight);
        } else if (surface_id == 1) {
          h2_exit_xy_upstream->Fill(x_mm, y_mm, weight);
        } else if (surface_id == 2) {
          h2_exit_yz_side_x->Fill(y_mm, z_mm, weight);
          const double xMid = 0.5 * (xMin + xMax);
          double sCoord = 0.0;
          if (x_mm < xMid) {
            sCoord = y_mm - yMin;
          } else {
            sCoord = yRange + (yMax - y_mm);
          }
          h2_exit_side_surface->Fill(sCoord, z_mm, weight);
        } else if (surface_id == 3) {
          h2_exit_xz_side_y->Fill(x_mm, z_mm, weight);
          const double yMid = 0.5 * (yMin + yMax);
          double sCoord = 2.0 * yRange;
          if (y_mm < yMid) {
            sCoord += x_mm - xMin;
          } else {
            sCoord += xRange + (xMax - x_mm);
          }
          h2_exit_side_surface->Fill(sCoord, z_mm, weight);
        }
        neutronTree->Fill();
      }
      neutronTree->Write();
      h2_exit_xy_downstream->Write();
      h2_exit_xy_upstream->Write();
      h2_exit_yz_side_x->Write();
      h2_exit_xz_side_y->Write();
      h2_exit_side_surface->Write();
    }
    if (!photonSurfaceHits.empty()) {
      auto* photonTree = new TTree("PhotonSurf", "Photon surface crossings");
      int event_id = 0;
      double E_MeV = 0.0;
      double x_mm = 0.0;
      double y_mm = 0.0;
      double z_mm = 0.0;
      double cosTheta = 0.0;
      double weight = 1.0;
      double time_ns = 0.0;
      int surface_id = 0;
      std::string surface_name;
      photonTree->Branch("event_id", &event_id);
      photonTree->Branch("E_MeV", &E_MeV);
      photonTree->Branch("x_mm", &x_mm);
      photonTree->Branch("y_mm", &y_mm);
      photonTree->Branch("z_mm", &z_mm);
      photonTree->Branch("cosTheta", &cosTheta);
      photonTree->Branch("weight", &weight);
      photonTree->Branch("time_ns", &time_ns);
      photonTree->Branch("surface_id", &surface_id);
      photonTree->Branch("surface_name", &surface_name);
      for (const auto& hit : photonSurfaceHits) {
        event_id = hit.event_id;
        E_MeV = hit.E_MeV;
        x_mm = hit.x_mm;
        y_mm = hit.y_mm;
        z_mm = hit.z_mm;
        cosTheta = hit.cosTheta;
        weight = hit.weight;
        time_ns = hit.time_ns;
        surface_id = hit.surface_id;
        switch (surface_id) {
          case 0: surface_name = "downstream"; break;
          case 1: surface_name = "upstream"; break;
          case 2: surface_name = "side_x"; break;
          case 3: surface_name = "side_y"; break;
          default: surface_name = "unknown"; break;
        }
        photonTree->Fill();
      }
      photonTree->Write();
    }
    runTree_->Write();
    rootFile_->Close();
    delete rootFile_;
    rootFile_ = nullptr;
    runTree_ = nullptr;
  }
#endif
  const auto out = base / "logs" / "run_summary.json";
  std::ofstream os(out);
  os << "{\n";
  os << "  \"edep_substrate\": " << (edepSub / MeV) << ",\n";
  os << "  \"edep_coating\": " << (edepCoat / MeV) << ",\n";
  os << "  \"nGamma\": " << nGammaTotal << ",\n";
  os << "  \"nGammaAbove5MeV\": " << nGammaTotal << ",\n";
  os << "  \"nNeutron\": " << nNeutronTotal << ",\n";
  os << "  \"nNeutronExit\": " << nNeutronExitTotal << ",\n";
  os << "  \"nNeutronModelExit\": " << nNeutronModelExitTotal << ",\n";
  os << "  \"nEvents\": " << (run ? run->GetNumberOfEvent() : config_.run.nEvents) << ",\n";
  os << "  \"physicsListName\": \"" << config_.physics.physicsListName << "\",\n";
  os << "  \"heatmap_file\": \"logs/heatmaps.json\",\n";
  os << "  \"plate_summary\": {\n";
  os << "    \"count\": " << plateEdepMeV.size() << ",\n";
  os << "    \"edep_MeV\": [";
  for (size_t i = 0; i < plateEdepMeV.size(); ++i) {
    os << plateEdepMeV[i];
    if (i + 1 < plateEdepMeV.size()) os << ", ";
  }
  os << "],\n";
  os << "    \"neutron_track_len_mm\": [";
  for (size_t i = 0; i < plateNeutronTrackLenMM.size(); ++i) {
    os << plateNeutronTrackLenMM[i];
    if (i + 1 < plateNeutronTrackLenMM.size()) os << ", ";
  }
  os << "]\n";
  os << "  },\n";
  os << "  \"beam\": {\n";
  os << "    \"E0_MeV\": " << config_.beam.energy_MeV << ",\n";
  os << "    \"energy_spread_model\": \"" << config_.beam.energy_spread_model << "\",\n";
  os << "    \"energy_sigma_rel_1sigma\": " << config_.beam.energy_sigma_rel_1sigma << ",\n";
  os << "    \"mode\": \"" << config_.beam.mode << "\",\n";
  os << "    \"pulse_width_us\": " << config_.beam.pulse_width_us << ",\n";
  os << "    \"rep_rate_Hz\": " << config_.beam.rep_rate_Hz << ",\n";
  os << "    \"I_pulse_A\": " << config_.beam.I_pulse_A << ",\n";
  os << "    \"I_avg_A_input\": " << config_.beam.I_avg_A << ",\n";
  os << "    \"I_avg_A_calc\": " << (config_.beam.I_pulse_A * (config_.beam.pulse_width_us * 1e-6) * config_.beam.rep_rate_Hz) << ",\n";
  os << "    \"beam_power_kW\": " << config_.beam.beam_power_kW << "\n";
  os << "  }\n";
  os << "}\n";


  if (config_.run.enableSwellingOutput) {
    const auto runMetaOut = base / "logs" / "run_meta.json";
    std::ofstream runMetaOs(runMetaOut);
    runMetaOs << "{\n";
    runMetaOs << "  \"nEvents\": " << nEventsForNorm << ",\n";
    runMetaOs << "  \"normalization_mode\": \"per_primary\",\n";
    runMetaOs << "  \"beam_current_A\": " << beamCurrentA << ",\n";
    runMetaOs << "  \"irradiation_time_s\": " << irradiationTimeS << ",\n";
    runMetaOs << "  \"electrons_per_second\": " << electronsPerSecond << ",\n";
    runMetaOs << "  \"electrons_total\": " << electronsTotal << ",\n";
    runMetaOs << "  \"beam_energy_MeV\": " << config_.beam.energy_MeV << ",\n";
    runMetaOs << "  \"target_type\": \"" << config_.target.type << "\",\n";
    runMetaOs << "  \"damage_source_status\": \"derived_from=damage_proxy\",\n";
    runMetaOs << "  \"units\": {\n";
    runMetaOs << "    \"edep_MeV_per_primary\": \"MeV\",\n";
    runMetaOs << "    \"edep_J_per_primary\": \"J\",\n";
    runMetaOs << "    \"damage_energy_eV_per_primary\": \"eV\",\n";
    runMetaOs << "    \"flux_n_groups\": \"1/cm2 per-primary\"\n";
    runMetaOs << "  },\n";
    runMetaOs << "  \"neutron_energy_group_edges_MeV\": [";
    for (size_t i = 0; i < neutronGroupEdgesMeV.size(); ++i) {
      runMetaOs << neutronGroupEdgesMeV[i];
      if (i + 1 < neutronGroupEdgesMeV.size()) runMetaOs << ", ";
    }
    runMetaOs << "],\n";
    runMetaOs << "  \"neutron_energy_group_edges_linear_MeV\": [";
    for (size_t i = 0; i < neutronGroupEdgesLinearMeV.size(); ++i) {
      runMetaOs << neutronGroupEdgesLinearMeV[i];
      if (i + 1 < neutronGroupEdgesLinearMeV.size()) runMetaOs << ", ";
    }
    runMetaOs << "]\n";
    runMetaOs << "}\n";

    const auto zBoundsForMesh = DeterminePlateStackZBounds(config_);
    const auto meshDefOut = base / "logs" / "mesh_definition.json";
    std::ofstream meshDefOs(meshDefOut);
    meshDefOs << "{\n";
    meshDefOs << "  \"mesh_id\": 1,\n";
    meshDefOs << "  \"mesh_name\": \"target_plate_stack\",\n";
    meshDefOs << "  \"bins\": {\"x\": " << edepBinsX_ << ", \"y\": " << edepBinsY_ << ", \"z\": " << edepBinsZ_ << "},\n";
    meshDefOs << "  \"bounds_mm\": {\"x_min\": " << edepBounds_.xMinMm << ", \"x_max\": " << edepBounds_.xMaxMm
              << ", \"y_min\": " << edepBounds_.yMinMm << ", \"y_max\": " << edepBounds_.yMaxMm << ", \"z_min\": "
              << zBoundsForMesh.first << ", \"z_max\": " << zBoundsForMesh.second << "},\n";
    meshDefOs << "  \"voxel_count\": " << meshRows.size() << ",\n";
    meshDefOs << "  \"notes\": \"damage_energy is proxy from edep until dedicated displacement scoring is enabled\"\n";
    meshDefOs << "}\n";

    const auto meshDataOut = base / "logs" / "mesh_data.csv";
    std::ofstream meshCsv(meshDataOut);
    meshCsv << "mesh_id,voxel_id,ix,iy,iz,x_mm,y_mm,z_mm,volume_cm3,material_id,material_name,volume_id,volume_name,layer_id,layer_name,"
            << "edep_MeV_per_primary,edep_J_per_primary,damage_energy_eV_per_primary,n_events_scored,flux_n_total_per_cm2_per_primary";
    for (int g = 0; g < kNeutronFluxGroups; ++g) {
      meshCsv << ",flux_n_g" << g;
    }
    meshCsv << ",h_prod_per_primary,he_prod_per_primary\n";
    for (const auto& row : meshRows) {
      meshCsv << row.mesh_id << "," << row.voxel_id << "," << row.ix << "," << row.iy << "," << row.iz << ","
              << row.x_mm << "," << row.y_mm << "," << row.z_mm << "," << row.volume_cm3 << "," << row.label.material_id
              << "," << row.label.material_name << "," << row.label.volume_id << "," << row.label.volume_name << ","
              << row.label.layer_id << "," << row.label.layer_name << "," << row.edep_MeV_per_primary << ","
              << row.edep_J_per_primary << "," << row.damage_energy_eV_per_primary << "," << row.n_events_scored << ","
              << row.flux_n_total_per_cm2_per_primary;
      for (double v : row.flux_n_groups_per_cm2_per_primary) {
        meshCsv << "," << v;
      }
      meshCsv << "," << row.h_prod_per_primary << "," << row.he_prod_per_primary << "\n";
    }


  }

  const auto heatmapOut = base / "logs" / "heatmaps.json";
  std::ofstream heatmapOs(heatmapOut);
  heatmapOs << "{\n";
  heatmapOs << "  \"plate_neutron_heatmap\": {\n";
  heatmapOs << "    \"bins_x\": " << plateHeatmapBinsX_ << ",\n";
  heatmapOs << "    \"bins_y\": " << plateHeatmapBinsY_ << ",\n";
  heatmapOs << "    \"x_min_mm\": " << plateHeatmapBounds_.xMinMm << ",\n";
  heatmapOs << "    \"x_max_mm\": " << plateHeatmapBounds_.xMaxMm << ",\n";
  heatmapOs << "    \"y_min_mm\": " << plateHeatmapBounds_.yMinMm << ",\n";
  heatmapOs << "    \"y_max_mm\": " << plateHeatmapBounds_.yMaxMm << ",\n";
  heatmapOs << "    \"plates\": [\n";
  const size_t plateStride =
      static_cast<size_t>(plateHeatmapBinsX_) * static_cast<size_t>(plateHeatmapBinsY_);
  for (size_t plateIdx = 0; plateIdx < plateCount_; ++plateIdx) {
    heatmapOs << "      {\n";
    heatmapOs << "        \"index\": " << (plateIdx + 1) << ",\n";
    heatmapOs << "        \"values\": [";
    const size_t offset = plateIdx * plateStride;
    for (size_t i = 0; i < plateStride; ++i) {
      const size_t idx = offset + i;
      if (idx < plateNeutronHeatmapMM.size()) {
        heatmapOs << plateNeutronHeatmapMM[idx];
      } else {
        heatmapOs << 0.0;
      }
      if (i + 1 < plateStride) heatmapOs << ", ";
    }
    heatmapOs << "]\n";
    heatmapOs << "      }";
    if (plateIdx + 1 < plateCount_) heatmapOs << ",";
    heatmapOs << "\n";
  }
  heatmapOs << "    ]\n";
  heatmapOs << "  },\n";
  heatmapOs << "  \"edep_3d\": {\n";
  heatmapOs << "    \"bins_x\": " << edepBinsX_ << ",\n";
  heatmapOs << "    \"bins_y\": " << edepBinsY_ << ",\n";
  heatmapOs << "    \"bins_z\": " << edepBinsZ_ << ",\n";
  heatmapOs << "    \"x_min_mm\": " << edepBounds_.xMinMm << ",\n";
  heatmapOs << "    \"x_max_mm\": " << edepBounds_.xMaxMm << ",\n";
  heatmapOs << "    \"y_min_mm\": " << edepBounds_.yMinMm << ",\n";
  heatmapOs << "    \"y_max_mm\": " << edepBounds_.yMaxMm << ",\n";
  heatmapOs << "    \"z_min_mm\": " << edepBounds_.zMinMm << ",\n";
  heatmapOs << "    \"z_max_mm\": " << edepBounds_.zMaxMm << ",\n";
  heatmapOs << "    \"values\": [";
  for (size_t i = 0; i < edep3dMeV.size(); ++i) {
    heatmapOs << edep3dMeV[i];
    if (i + 1 < edep3dMeV.size()) heatmapOs << ", ";
  }
  heatmapOs << "]\n";
  heatmapOs << "  }\n";
  heatmapOs << "}\n";
}

void RunAction::AccumulateEvent(double edepSubstrate,
                                double edepCoating,
                                int nGamma,
                                int nNeutron,
                                int nNeutronExit,
                                int nNeutronModelExit,
                                int eventId,
                                const std::vector<double>& plateEdep,
                                const std::vector<double>& plateNeutronTrackLen,
                                const std::vector<double>& plateNeutronHeatmap,
                                const std::vector<double>& edep3d,
                                const std::vector<NeutronSurfaceHit>& neutronSurfaceHits,
                                const std::vector<PhotonSurfaceHit>& photonSurfaceHits,
                                const std::vector<double>& plateNiel,
                                const std::vector<double>& plateGasH,
                                const std::vector<double>& plateGasHe) {
  std::lock_guard<std::mutex> lock(gRunSummaryMutex);
  (void)eventId;
  gTotalEdepSubstrate += edepSubstrate;
  gTotalEdepCoating += edepCoating;
  gTotalGamma += nGamma;
  gTotalNeutron += nNeutron;
  gTotalNeutronExit += nNeutronExit;
  gTotalNeutronModelExit += nNeutronModelExit;
  if (gTotalPlateEdep.size() < plateEdep.size()) {
    gTotalPlateEdep.resize(plateEdep.size(), 0.0);
  }
  if (gTotalPlateNeutronTrackLen.size() < plateNeutronTrackLen.size()) {
    gTotalPlateNeutronTrackLen.resize(plateNeutronTrackLen.size(), 0.0);
  }
  for (size_t i = 0; i < plateEdep.size(); ++i) {
    gTotalPlateEdep[i] += plateEdep[i];
  }
  for (size_t i = 0; i < plateNeutronTrackLen.size(); ++i) {
    gTotalPlateNeutronTrackLen[i] += plateNeutronTrackLen[i];
  }
  if (gTotalPlateNeutronHeatmap.size() < plateNeutronHeatmap.size()) {
    gTotalPlateNeutronHeatmap.resize(plateNeutronHeatmap.size(), 0.0);
  }
  for (size_t i = 0; i < plateNeutronHeatmap.size(); ++i) {
    gTotalPlateNeutronHeatmap[i] += plateNeutronHeatmap[i];
  }
  if (gTotalEdep3d.size() < edep3d.size()) {
    gTotalEdep3d.resize(edep3d.size(), 0.0);
  }
  for (size_t i = 0; i < edep3d.size(); ++i) {
    gTotalEdep3d[i] += edep3d[i];
  }
  if (!neutronSurfaceHits.empty()) {
    gNeutronSurfaceHits.insert(gNeutronSurfaceHits.end(), neutronSurfaceHits.begin(), neutronSurfaceHits.end());
  }
  if (!photonSurfaceHits.empty()) {
    gPhotonSurfaceHits.insert(gPhotonSurfaceHits.end(), photonSurfaceHits.begin(), photonSurfaceHits.end());
  }
  if (gTotalPlateNiel.size() < plateNiel.size()) {
    gTotalPlateNiel.resize(plateNiel.size(), 0.0);
  }
  for (size_t i = 0; i < plateNiel.size(); ++i) {
    gTotalPlateNiel[i] += plateNiel[i];
  }
  if (gTotalPlateGasH.size() < plateGasH.size()) {
    gTotalPlateGasH.resize(plateGasH.size(), 0.0);
  }
  for (size_t i = 0; i < plateGasH.size(); ++i) {
    gTotalPlateGasH[i] += plateGasH[i];
  }
  if (gTotalPlateGasHe.size() < plateGasHe.size()) {
    gTotalPlateGasHe.resize(plateGasHe.size(), 0.0);
  }
  for (size_t i = 0; i < plateGasHe.size(); ++i) {
    gTotalPlateGasHe[i] += plateGasHe[i];
  }
  // TODO: migrate to G4Accumulable/G4AnalysisManager for richer MT-safe reporting.
}
