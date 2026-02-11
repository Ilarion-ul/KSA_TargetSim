#include "DetectorConstruction.h"

#include <G4Box.hh>
#include <G4Colour.hh>
#include <G4Element.hh>
#include <G4Exception.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4String.hh>
#include <G4SubtractionSolid.hh>
#include <G4SystemOfUnits.hh>
#include <G4Tubs.hh>
#include <G4VisAttributes.hh>

#include <algorithm>
#include <sstream>

DetectorConstruction::DetectorConstruction(AppConfig config) : config_(std::move(config)) {}

namespace {
void ConfigureVis(G4LogicalVolume* lv, const G4Colour& c) {
  auto* vis = new G4VisAttributes(c);
  vis->SetForceSolid(true);
  lv->SetVisAttributes(vis);
}

G4String IndexedName(const char* base, int idx1) {
  std::ostringstream os;
  os << base << "_" << idx1;
  return os.str();
}

std::vector<double> BuildUMoInterPlateGapsMm(const AppConfig& config) {
  const auto& t = config.target;
  const size_t n = t.plate_thicknesses_mm.size();
  if (n <= 1) return {};
  if (!t.inter_plate_gaps_mm.empty()) {
    return t.inter_plate_gaps_mm;
  }
  std::vector<double> gaps(n - 1, t.gap_rear_mm);
  const size_t split = static_cast<size_t>(std::min<int>(t.gap_split_index, static_cast<int>(n - 1)));
  for (size_t i = 0; i < split; ++i) {
    gaps[i] = t.gap_front_mm;
  }
  // backward-compatible fallback for legacy 12-plate setups
  if (n == 12 && t.gap_mid_mm > 0.0 && t.gap_front_mm == 3.0 && t.gap_rear_mm == 1.75) {
    std::fill(gaps.begin(), gaps.end(), t.gap_mid_mm);
  }
  return gaps;
}
} // namespace

G4VPhysicalVolume* DetectorConstruction::Construct() {
  auto* nist = G4NistManager::Instance();

  auto* worldMat = nist->FindOrBuildMaterial("G4_Galactic");
  auto* waterMat = nist->FindOrBuildMaterial("G4_WATER"); // TODO: replace with demineralized-water properties.

  // ---------------------------------------------------------------------------
  // U-Mo branch: realistic target with construction materials and beamline.
  // Coordinate convention: +Z is downstream (electron direction).
  // Origin (0,0,0) is the geometric center of TargetAssembly.
  // ---------------------------------------------------------------------------
  if (config_.target.type == "U-Mo") {
    // Materials.
    auto* elU = nist->FindOrBuildElement("U");
    auto* elMo = nist->FindOrBuildElement("Mo");
    auto* g4Al = nist->FindOrBuildMaterial("G4_Al");
    auto* heMat = nist->FindOrBuildMaterial("G4_He"); // TODO: pressure/density tuning for operational conditions.

    auto* u7mo = G4Material::GetMaterial("U7Mo", false);
    if (!u7mo) {
      u7mo = new G4Material("U7Mo", config_.target.u7mo_density_g_cm3 * g / cm3, 2);
      u7mo->AddElement(elU, 0.93);
      u7mo->AddElement(elMo, 0.07);
    }

    auto* sav1 = G4Material::GetMaterial("SAV1", false);
    if (!sav1) {
      // TODO: replace by real SAV-1 composition when available.
      sav1 = new G4Material("SAV1", g4Al->GetDensity(), 1);
      sav1->AddElement(nist->FindOrBuildElement("Al"), 1.0);
    }

    const double totalAssemblyLen = config_.geometry.total_assembly_len_mm * mm;
    const double beamlineVacLen = config_.geometry.beamline_vacuum_len_mm * mm;
    const double targetRegionLen = totalAssemblyLen - beamlineVacLen;
    const double clearance = config_.geometry.target_region_extra_clearance_mm * mm;

    if (targetRegionLen <= 0.0) {
      G4Exception("DetectorConstruction::Construct", "GeomUMo001", FatalException,
                  "TargetRegionLength must be > 0 (check total_assembly_len_mm and beamline_vacuum_len_mm)");
    }

    const auto& t = config_.target;
    const auto& plateTs = t.plate_thicknesses_mm;

    const double plateHalfXY = 0.5 * t.plate_xy_mm * mm;
    const double cladFront = t.clad_thickness_front_mm * mm;
    const double cladRest = t.clad_thickness_rest_mm * mm;
    const double gapInOut = t.gap_inout_mm * mm;
    const auto interPlateGapsMm = BuildUMoInterPlateGapsMm(config_);
    const double housingInnerHalfXY = 0.5 * t.housing_inner_xy_mm * mm;
    const double housingWall = t.housing_wall_mm * mm;
    const double housingOuterHalfXY = housingInnerHalfXY + housingWall;
    const double entranceWindowT = t.entrance_window_mm * mm;
    const double heliumLen = t.helium_chamber_len_mm * mm;

    const double maxCladXY = plateHalfXY + std::max(cladFront, cladRest);
    if (maxCladXY > housingInnerHalfXY) {
      G4Exception("DetectorConstruction::Construct", "GeomUMo002", FatalException,
                  "Plate with cladding does not fit housing inner aperture");
    }

    // Plate stack thickness (inlet gap + plates + inter-plate gaps + outlet gap).
    double stackLen = gapInOut;
    for (size_t i = 0; i < plateTs.size(); ++i) {
      const double coreT = plateTs[i] * mm;
      const double clad = (i < 4) ? cladFront : cladRest;
      stackLen += coreT + 2.0 * clad;
      if (i + 1 < plateTs.size()) {
        const double g = (i < interPlateGapsMm.size()) ? interPlateGapsMm[i] * mm : gapInOut;
        stackLen += g;
      }
    }
    stackLen += gapInOut;

    const double requiredTargetLen = entranceWindowT + stackLen + heliumLen + 2.0 * clearance;

    G4cout << "[geom U-Mo] Plate count = " << plateTs.size() << "\n"
           << "[geom U-Mo] BeamlineVacuumLength = " << beamlineVacLen / mm << " mm\n"
           << "[geom U-Mo] TotalAssemblyLength = " << totalAssemblyLen / mm << " mm\n"
           << "[geom U-Mo] TargetRegionLength = " << targetRegionLen / mm << " mm\n"
           << "[geom U-Mo] Plate stack length (with cladding+gaps) = " << stackLen / mm << " mm\n"
           << "[geom U-Mo] Required target length (window+stack+He+clearance) = " << requiredTargetLen / mm << " mm\n"
           << "[geom U-Mo] Housing inner XY = " << (2.0 * housingInnerHalfXY) / mm << " mm, plate+clad XY = "
           << (2.0 * maxCladXY) / mm << " mm" << G4endl;

    if (requiredTargetLen > targetRegionLen) {
      G4Exception("DetectorConstruction::Construct", "GeomUMo003", FatalException,
                  "Target region is too short for window + plate stack + helium chamber");
    }

    // World and top assembly.
    const double margin = config_.geometry.worldMargin_mm * mm;
    auto* worldS = new G4Box("World", housingOuterHalfXY + margin, housingOuterHalfXY + margin,
                             0.5 * totalAssemblyLen + margin);
    auto* worldLV = new G4LogicalVolume(worldS, worldMat, "World");
    auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr, false, 0, true);

    auto* targetAssemblyS = new G4Box("TargetAssemblyShape", housingOuterHalfXY, housingOuterHalfXY, 0.5 * totalAssemblyLen);
    auto* targetAssemblyLV = new G4LogicalVolume(targetAssemblyS, worldMat, "TargetAssemblyLV");
    new G4PVPlacement(nullptr, {}, targetAssemblyLV, "TargetAssembly", worldLV, false, 0, true);

    // Beamline vacuum section (upstream segment).
    const double beamlineCenterZ = -0.5 * totalAssemblyLen + 0.5 * beamlineVacLen;
    auto* beamlineS = new G4Box("BeamlineVacuumShape", housingInnerHalfXY, housingInnerHalfXY, 0.5 * beamlineVacLen);
    auto* beamlineLV = new G4LogicalVolume(beamlineS, worldMat, "BeamlineVacuumLV");
    new G4PVPlacement(nullptr, {0, 0, beamlineCenterZ}, beamlineLV, "BeamlineVacuum", targetAssemblyLV, false, 0, true);

    // Target region starts right after beamline section.
    const double targetStartZ = -0.5 * totalAssemblyLen + beamlineVacLen;
    const double targetCenterZ = targetStartZ + 0.5 * targetRegionLen;

    // Target housing: hollow square tube.
    auto* outerS = new G4Box("TargetHousingOuterShape", housingOuterHalfXY, housingOuterHalfXY, 0.5 * targetRegionLen);
    auto* innerS = new G4Box("TargetHousingInnerShape", housingInnerHalfXY, housingInnerHalfXY, 0.5 * targetRegionLen + 0.1 * mm);
    auto* housingS = new G4SubtractionSolid("TargetHousingShape", outerS, innerS);
    auto* housingLV = new G4LogicalVolume(housingS, sav1, "TargetHousingLV");
    new G4PVPlacement(nullptr, {0, 0, targetCenterZ}, housingLV, "TargetHousing", targetAssemblyLV, false, 0, true);

    // Inner target medium container in housing aperture.
    auto* fillMat = (t.fill_medium_in_target == "vacuum") ? worldMat : waterMat;
    auto* targetInnerFluidS = new G4Box("TargetInnerFluidShape", housingInnerHalfXY, housingInnerHalfXY, 0.5 * targetRegionLen);
    auto* targetInnerFluidLV = new G4LogicalVolume(targetInnerFluidS, fillMat, "TargetInnerFluidLV");
    new G4PVPlacement(nullptr, {0, 0, targetCenterZ}, targetInnerFluidLV, "TargetInnerFluid", targetAssemblyLV, false, 0, true);

    if (config_.geometry.enable_alignment_pin) {
      auto* pinMat = nist->FindOrBuildMaterial(config_.geometry.alignment_pin_material);
      if (!pinMat) pinMat = sav1;
      const auto& ps = config_.geometry.alignment_pin_size_mm;
      const auto& pp = config_.geometry.alignment_pin_pos_mm;
      const double hx = 0.5 * ps[0] * mm;
      const double hy = 0.5 * ps[1] * mm;
      const double hz = 0.5 * ps[2] * mm;
      if (std::abs(pp[0] * mm) + hx > housingInnerHalfXY || std::abs(pp[1] * mm) + hy > housingInnerHalfXY ||
          std::abs(pp[2] * mm) + hz > 0.5 * targetRegionLen) {
        G4Exception("DetectorConstruction::Construct", "GeomPIN001", FatalException,
                    "Alignment pin is outside U-Mo target inner volume");
      }
      auto* pinS = new G4Box("AlignmentPinShape", hx, hy, hz);
      auto* pinLV = new G4LogicalVolume(pinS, pinMat, "AlignmentPinLV");
      new G4PVPlacement(nullptr, {pp[0] * mm, pp[1] * mm, pp[2] * mm}, pinLV, "AlignmentPin", targetInnerFluidLV, false, 0,
                        true);
      ConfigureVis(pinLV, G4Colour(0.95, 0.1, 0.1));
    }

    // Build internal sequence in TargetInnerFluid local coordinates.
    double zCursor = -0.5 * targetRegionLen + clearance;

    // Entrance window at vacuum -> target interface.
    {
      const double hz = 0.5 * entranceWindowT;
      auto* winS = new G4Box("EntranceWindowShape", housingInnerHalfXY, housingInnerHalfXY, hz);
      auto* winLV = new G4LogicalVolume(winS, sav1, "EntranceWindowLV");
      new G4PVPlacement(nullptr, {0, 0, zCursor + hz}, winLV, "EntranceWindow", targetInnerFluidLV, false, 0, true);
      zCursor += entranceWindowT;
      ConfigureVis(winLV, G4Colour(0.85, 0.85, 0.85));
    }

    int gapIdx = 1;
    for (size_t i = 0; i < plateTs.size(); ++i) {
      // Gap before each plate: inlet for first, otherwise inter-plate scheme.
      {
        const double gapT = (i == 0) ? gapInOut : ((i - 1 < interPlateGapsMm.size()) ? interPlateGapsMm[i - 1] * mm : gapInOut);
        const double hz = 0.5 * gapT;
        auto* gapS = new G4Box("WaterGapShape", housingInnerHalfXY, housingInnerHalfXY, hz);
        auto* gapLV = new G4LogicalVolume(gapS, waterMat, "WaterGapLV");
        new G4PVPlacement(nullptr, {0, 0, zCursor + hz}, gapLV, IndexedName("WaterGap", gapIdx), targetInnerFluidLV,
                          false, gapIdx, true);
        zCursor += gapT;
        ++gapIdx;
        ConfigureVis(gapLV, G4Colour(0.6, 0.85, 1.0, 0.35));
      }

      const double coreT = plateTs[i] * mm;
      const double clad = (i < 4) ? cladFront : cladRest;
      const int i1 = static_cast<int>(i) + 1;

      // Cladding envelope.
      const double cladHalfXY = plateHalfXY + clad;
      const double cladTotT = coreT + 2.0 * clad;
      auto* cladS = new G4Box("PlateCladShape", cladHalfXY, cladHalfXY, 0.5 * cladTotT);
      auto* cladLV = new G4LogicalVolume(cladS, sav1, "PlateCladAlLV");
      new G4PVPlacement(nullptr, {0, 0, zCursor + 0.5 * cladTotT}, cladLV, IndexedName("PlateCladAl", i1),
                        targetInnerFluidLV, false, i1, true);
      if (i == 0) coatingLV_ = cladLV;

      // U-Mo core.
      auto* coreS = new G4Box("PlateUShape", plateHalfXY, plateHalfXY, 0.5 * coreT);
      auto* coreLV = new G4LogicalVolume(coreS, u7mo, "PlateULV");
      new G4PVPlacement(nullptr, {}, coreLV, IndexedName("PlateU", i1), cladLV, false, i1, true);
      if (i == 0) substrateLV_ = coreLV;

      ConfigureVis(cladLV, G4Colour(0.9, 0.9, 0.9));
      ConfigureVis(coreLV, G4Colour(0.3, 0.5, 0.2));

      zCursor += cladTotT;
    }

    // Final 13th (outlet) gap.
    {
      const double hz = 0.5 * gapInOut;
      auto* gapS = new G4Box("WaterGapShapeFinal", housingInnerHalfXY, housingInnerHalfXY, hz);
      auto* gapLV = new G4LogicalVolume(gapS, waterMat, "WaterGapLV");
      new G4PVPlacement(nullptr, {0, 0, zCursor + hz}, gapLV, IndexedName("WaterGap", gapIdx), targetInnerFluidLV,
                        false, gapIdx, true);
      zCursor += gapInOut;
      ConfigureVis(gapLV, G4Colour(0.6, 0.85, 1.0, 0.35));
    }

    // Helium chamber downstream (replaces water in this segment).
    {
      const double hz = 0.5 * heliumLen;
      auto* heS = new G4Box("HeliumChamberShape", housingInnerHalfXY, housingInnerHalfXY, hz);
      auto* heLV = new G4LogicalVolume(heS, heMat, "HeliumChamberLV");
      new G4PVPlacement(nullptr, {0, 0, zCursor + hz}, heLV, "HeliumChamber", targetInnerFluidLV, false, 0, true);
      ConfigureVis(heLV, G4Colour(1.0, 0.95, 0.6, 0.25));
    }

    ConfigureVis(targetAssemblyLV, G4Colour(0.7, 0.7, 0.7, 0.1));
    ConfigureVis(beamlineLV, G4Colour(0.5, 0.5, 0.8, 0.1));
    ConfigureVis(housingLV, G4Colour(0.8, 0.8, 0.85));
    ConfigureVis(targetInnerFluidLV, G4Colour(0.6, 0.85, 1.0, 0.1));
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    return worldPV;
  }

  // ---------------------------------------------------------------------------
  // W-Ta sectional branch (existing realistic model).
  // ---------------------------------------------------------------------------
  if (config_.target.type == "W-Ta" && config_.geometry.simpleCylinder) {
    auto* wMat = nist->FindOrBuildMaterial("G4_W");
    if (config_.target.w_substrate_material == "W_Fe_Ni") {
      auto* elW = nist->FindOrBuildElement("W");
      auto* elFe = nist->FindOrBuildElement("Fe");
      auto* elNi = nist->FindOrBuildElement("Ni");
      auto* wFeNi = G4Material::GetMaterial("W_Fe_Ni", false);
      if (!wFeNi) {
        wFeNi = new G4Material("W_Fe_Ni", 18.5 * g / cm3, 3);
        wFeNi->AddElement(elW, 0.955);
        wFeNi->AddElement(elFe, 0.015);
        wFeNi->AddElement(elNi, 0.030);
      }
      wMat = wFeNi;
    }
    auto* taMat = nist->FindOrBuildMaterial("G4_Ta");
    auto* tiMat = nist->FindOrBuildMaterial("G4_Ti");

    const double plateXY = config_.target.plate_xy_mm * mm;
    const double plateHalfXY = 0.5 * plateXY;
    const double ta = config_.target.clad_ta_mm * mm;
    const double ti = config_.target.buffer_ti_mm * mm;
    const double waterGap = config_.target.water_gap_mm * mm;
    const double assemblyT = config_.target.assembly_thickness_mm * mm;
    const auto& wPlateTsMM = config_.target.plate_thicknesses_mm;

    const double taHalfX = plateHalfXY + ti + ta;
    const double taHalfY = plateHalfXY + ti + ta;
    const double margin = config_.geometry.worldMargin_mm * mm;

    auto* worldS = new G4Box("World", taHalfX + margin, taHalfY + margin, 0.5 * assemblyT + margin);
    auto* worldLV = new G4LogicalVolume(worldS, worldMat, "World");
    auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr, false, 0, true);

    auto* assemblyS = new G4Box("TargetAssembly", taHalfX, taHalfY, 0.5 * assemblyT);
    auto* assemblyLV = new G4LogicalVolume(assemblyS, waterMat, "TargetAssemblyLV");
    new G4PVPlacement(nullptr, {}, assemblyLV, "TargetAssembly", worldLV, false, 0, true);

    if (config_.geometry.enable_alignment_pin) {
      auto* pinMat = nist->FindOrBuildMaterial(config_.geometry.alignment_pin_material);
      if (!pinMat) pinMat = taMat;
      const auto& ps = config_.geometry.alignment_pin_size_mm;
      const auto& pp = config_.geometry.alignment_pin_pos_mm;
      const double hx = 0.5 * ps[0] * mm;
      const double hy = 0.5 * ps[1] * mm;
      const double hz = 0.5 * ps[2] * mm;
      if (std::abs(pp[0] * mm) + hx > taHalfX || std::abs(pp[1] * mm) + hy > taHalfY ||
          std::abs(pp[2] * mm) + hz > 0.5 * assemblyT) {
        G4Exception("DetectorConstruction::Construct", "GeomPIN002", FatalException,
                    "Alignment pin is outside W-Ta assembly volume");
      }
      auto* pinS = new G4Box("AlignmentPinShape", hx, hy, hz);
      auto* pinLV = new G4LogicalVolume(pinS, pinMat, "AlignmentPinLV");
      new G4PVPlacement(nullptr, {pp[0] * mm, pp[1] * mm, pp[2] * mm}, pinLV, "AlignmentPin", assemblyLV, false, 0, true);
      ConfigureVis(pinLV, G4Colour(0.95, 0.1, 0.1));
    }

    double stackT = 8.0 * waterGap;
    for (double tmm : wPlateTsMM) {
      stackT += tmm * mm + 2.0 * ti + 2.0 * ta;
    }
    G4cout << "[geom W-Ta] sectional stack total thickness = " << stackT / mm
           << " mm; assembly thickness = " << assemblyT / mm << " mm" << G4endl;
    if (stackT > assemblyT) {
      G4Exception("DetectorConstruction::Construct", "GeomWTa001", FatalException,
                  "Sectional stack thickness exceeds assembly_thickness_mm");
    }

    double zCursor = -0.5 * stackT;
    for (size_t i = 0; i < wPlateTsMM.size(); ++i) {
      {
        const double hz = 0.5 * waterGap;
        auto* s = new G4Box("TargetWaterGapShape", taHalfX, taHalfY, hz);
        auto* lv = new G4LogicalVolume(s, waterMat, "TargetWaterGap");
        new G4PVPlacement(nullptr, {0, 0, zCursor + hz}, lv, IndexedName("TargetWaterGap", static_cast<int>(i) + 1),
                          assemblyLV, false, static_cast<int>(i) + 1, true);
        zCursor += waterGap;
        ConfigureVis(lv, G4Colour(0.6, 0.85, 1.0, 0.35));
      }

      const double wT = wPlateTsMM[i] * mm;
      const double taTot = wT + 2.0 * ti + 2.0 * ta;
      const double tiTot = wT + 2.0 * ti;

      auto* taS = new G4Box("TargetCoatingShape", taHalfX, taHalfY, 0.5 * taTot);
      auto* taLV = new G4LogicalVolume(taS, taMat, "TargetCoating");
      new G4PVPlacement(nullptr, {0, 0, zCursor + 0.5 * taTot}, taLV,
                        IndexedName("TargetCoating", static_cast<int>(i) + 1), assemblyLV, false,
                        static_cast<int>(i) + 1, true);
      if (i == 0) coatingLV_ = taLV;

      auto* tiS = new G4Box("TargetBufferTiShape", plateHalfXY + ti, plateHalfXY + ti, 0.5 * tiTot);
      auto* tiLV = new G4LogicalVolume(tiS, tiMat, "TargetBufferTi");
      new G4PVPlacement(nullptr, {}, tiLV, IndexedName("TargetBufferTi", static_cast<int>(i) + 1), taLV, false,
                        static_cast<int>(i) + 1, true);

      auto* wS = new G4Box("TargetSubstrateShape", plateHalfXY, plateHalfXY, 0.5 * wT);
      auto* wLV = new G4LogicalVolume(wS, wMat, "TargetSubstrate");
      new G4PVPlacement(nullptr, {}, wLV, IndexedName("TargetSubstrate", static_cast<int>(i) + 1), tiLV, false,
                        static_cast<int>(i) + 1, true);
      if (i == 0) substrateLV_ = wLV;

      ConfigureVis(taLV, G4Colour(0.95, 0.5, 0.2));
      ConfigureVis(tiLV, G4Colour(0.7, 0.7, 0.75));
      ConfigureVis(wLV, G4Colour(0.2, 0.4, 0.9));

      zCursor += taTot;
    }

    {
      const double hz = 0.5 * waterGap;
      auto* s = new G4Box("TargetWaterGapShapeFinal", taHalfX, taHalfY, hz);
      auto* lv = new G4LogicalVolume(s, waterMat, "TargetWaterGap");
      new G4PVPlacement(nullptr, {0, 0, zCursor + hz}, lv, IndexedName("TargetWaterGap", 8), assemblyLV, false, 8,
                        true);
      ConfigureVis(lv, G4Colour(0.6, 0.85, 1.0, 0.35));
    }

    ConfigureVis(assemblyLV, G4Colour(0.7, 0.9, 0.95, 0.15));
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    return worldPV;
  }

  // ---------------------------------------------------------------------------
  // Legacy/simple target branch.
  // ---------------------------------------------------------------------------
  auto* substrateMat = nist->FindOrBuildMaterial(config_.target.type == "U-Al" ? "G4_U" : "G4_W");
  auto* coatingMat = nist->FindOrBuildMaterial(config_.target.type == "U-Al" ? "G4_Al" : "G4_Ta");

  const double r = config_.target.radius_mm * mm;
  const double subT = config_.target.substrate_thickness_mm * mm;
  const double coatT = config_.target.coating_thickness_mm * mm;
  const double totalTargetT = subT + coatT;
  const double margin = config_.geometry.worldMargin_mm * mm;

  const double worldXY = 2.0 * (r + margin);
  const double worldZ = 2.0 * (0.5 * totalTargetT + margin);

  auto* worldS = new G4Box("World", 0.5 * worldXY, 0.5 * worldXY, 0.5 * worldZ);
  auto* worldLV = new G4LogicalVolume(worldS, worldMat, "World");
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr, false, 0, true);

  auto* subS = new G4Tubs("TargetSubstrate", 0.0, r, 0.5 * subT, 0.0, 360.0 * deg);
  substrateLV_ = new G4LogicalVolume(subS, substrateMat, "TargetSubstrate");

  auto* coatS = new G4Tubs("TargetCoating", 0.0, r, 0.5 * coatT, 0.0, 360.0 * deg);
  coatingLV_ = new G4LogicalVolume(coatS, coatingMat, "TargetCoating");

  const double zSub = -0.5 * coatT;
  const double zCoat = +0.5 * subT;

  new G4PVPlacement(nullptr, {0, 0, zSub}, substrateLV_, "TargetSubstrate", worldLV, false, 0, true);
  new G4PVPlacement(nullptr, {0, 0, zCoat}, coatingLV_, "TargetCoating", worldLV, false, 0, true);

  ConfigureVis(substrateLV_, G4Colour(0.2, 0.4, 0.9));
  ConfigureVis(coatingLV_, G4Colour(0.9, 0.4, 0.2));
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  return worldPV;
}
