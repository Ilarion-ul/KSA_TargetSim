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
std::vector<double> gTotalPlateEdep;
std::vector<double> gTotalPlateNeutronTrackLen;
std::vector<double> gTotalPlateNeutronHeatmap;
std::vector<double> gTotalEdep3d;
std::vector<RunAction::NeutronSurfaceHit> gNeutronSurfaceHits;
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
    gTotalPlateEdep.assign(plateCount_, 0.0);
    gTotalPlateNeutronTrackLen.assign(plateCount_, 0.0);
    gTotalPlateNeutronHeatmap.assign(static_cast<size_t>(plateHeatmapBinsX_) * static_cast<size_t>(plateHeatmapBinsY_) *
                                         plateCount_,
                                     0.0);
    gTotalEdep3d.assign(static_cast<size_t>(edepBinsX_) * static_cast<size_t>(edepBinsY_) * static_cast<size_t>(edepBinsZ_),
                        0.0);
    gNeutronSurfaceHits.clear();
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
  std::vector<double> plateEdep;
  std::vector<double> plateNeutronTrackLen;
  std::vector<double> plateNeutronHeatmap;
  std::vector<double> edep3d;
  std::vector<NeutronSurfaceHit> neutronSurfaceHits;
  {
    std::lock_guard<std::mutex> lock(gRunSummaryMutex);
    edepSub = gTotalEdepSubstrate;
    edepCoat = gTotalEdepCoating;
    nGammaTotal = gTotalGamma;
    nNeutronTotal = gTotalNeutron;
    nNeutronExitTotal = gTotalNeutronExit;
    plateEdep = gTotalPlateEdep;
    plateNeutronTrackLen = gTotalPlateNeutronTrackLen;
    plateNeutronHeatmap = gTotalPlateNeutronHeatmap;
    edep3d = gTotalEdep3d;
    neutronSurfaceHits = gNeutronSurfaceHits;
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

#ifdef KSA_USE_ROOT
  if (rootFile_ && runTree_) {
    double edep_substrate_MeV = edepSub / MeV;
    double edep_coating_MeV = edepCoat / MeV;
    long long nGamma = nGammaTotal;
    long long nNeutron = nNeutronTotal;
    int nEvents = run ? run->GetNumberOfEvent() : config_.run.nEvents;
    std::string physicsListName = config_.physics.physicsListName;
    int nNeutronExit = static_cast<int>(nNeutronExitTotal);
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
    runTree_->Branch("nNeutron", &nNeutron);
    runTree_->Branch("nNeutronExit", &nNeutronExit);
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
      int nEvents = run ? run->GetNumberOfEvent() : config_.run.nEvents;
      int nThreads = config_.run.nThreads;
      double per_primary = 1.0;
      double N_e_per_s = I_avg_A / 1.602176634e-19;
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
      metaTree->Branch("nThreads", &nThreads);
      metaTree->Branch("per_primary", &per_primary);
      metaTree->Branch("N_e_per_s", &N_e_per_s);
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
      const double zMin = edepBounds_.zMinMm;
      const double zMax = edepBounds_.zMaxMm;
      auto* h3 = new TH3D("edep_3d", "Energy deposition 3D", edepBinsX_, xMin, xMax, edepBinsY_, yMin, yMax,
                          edepBinsZ_, zMin, zMax);
      for (int iz = 0; iz < edepBinsZ_; ++iz) {
        for (int iy = 0; iy < edepBinsY_; ++iy) {
          for (int ix = 0; ix < edepBinsX_; ++ix) {
            const size_t idx = (static_cast<size_t>(iz) * static_cast<size_t>(edepBinsY_) + static_cast<size_t>(iy)) *
                                   static_cast<size_t>(edepBinsX_) +
                               static_cast<size_t>(ix);
            if (idx < edep3dMeV.size()) {
              h3->SetBinContent(ix + 1, iy + 1, iz + 1, edep3dMeV[idx]);
            }
          }
        }
      }
      h3->Write();
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
      const double zMin = edepBounds_.zMinMm;
      const double zMax = edepBounds_.zMaxMm;
      auto* h2_exit_xy_downstream = new TH2D("h2_neutron_exit_xy_downstream", "Neutron exit (downstream)",
                                             kSurfBinsX, xMin, xMax, kSurfBinsY, yMin, yMax);
      auto* h2_exit_xy_upstream =
          new TH2D("h2_neutron_exit_xy_upstream", "Neutron exit (upstream)", kSurfBinsX, xMin, xMax, kSurfBinsY, yMin, yMax);
      auto* h2_exit_yz_side_x =
          new TH2D("h2_neutron_exit_yz_side_x", "Neutron exit (side x)", kSurfBinsY, yMin, yMax, kSurfBinsZ, zMin, zMax);
      auto* h2_exit_xz_side_y =
          new TH2D("h2_neutron_exit_xz_side_y", "Neutron exit (side y)", kSurfBinsX, xMin, xMax, kSurfBinsZ, zMin, zMax);
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
        } else if (surface_id == 3) {
          h2_exit_xz_side_y->Fill(x_mm, z_mm, weight);
        }
        neutronTree->Fill();
      }
      neutronTree->Write();
      h2_exit_xy_downstream->Write();
      h2_exit_xy_upstream->Write();
      h2_exit_yz_side_x->Write();
      h2_exit_xz_side_y->Write();
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
  os << "  \"nNeutron\": " << nNeutronTotal << ",\n";
  os << "  \"nNeutronExit\": " << nNeutronExitTotal << ",\n";
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
                                int eventId,
                                const std::vector<double>& plateEdep,
                                const std::vector<double>& plateNeutronTrackLen,
                                const std::vector<double>& plateNeutronHeatmap,
                                const std::vector<double>& edep3d,
                                const std::vector<NeutronSurfaceHit>& neutronSurfaceHits) {
  std::lock_guard<std::mutex> lock(gRunSummaryMutex);
  (void)eventId;
  gTotalEdepSubstrate += edepSubstrate;
  gTotalEdepCoating += edepCoating;
  gTotalGamma += nGamma;
  gTotalNeutron += nNeutron;
  gTotalNeutronExit += nNeutronExit;
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
  // TODO: migrate to G4Accumulable/G4AnalysisManager for richer MT-safe reporting.
}
