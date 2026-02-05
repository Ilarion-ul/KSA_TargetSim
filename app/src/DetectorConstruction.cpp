#include "DetectorConstruction.h"

#include <G4Box.hh>
#include <G4Colour.hh>
#include <G4Exception.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4String.hh>
#include <G4SystemOfUnits.hh>
#include <G4Tubs.hh>
#include <G4VisAttributes.hh>

#include <string>

DetectorConstruction::DetectorConstruction(AppConfig config) : config_(std::move(config)) {}

namespace {
void ConfigureVis(G4LogicalVolume* lv, const G4Colour& c) {
  auto* vis = new G4VisAttributes(c);
  vis->SetForceSolid(true);
  lv->SetVisAttributes(vis);
}

} // namespace

G4VPhysicalVolume* DetectorConstruction::Construct() {
  auto* nist = G4NistManager::Instance();

  auto* worldMat = nist->FindOrBuildMaterial("G4_Galactic");
  auto* waterMat = nist->FindOrBuildMaterial("G4_WATER"); // TODO: replace by demineralized-water composition.

  // Keep legacy simple cylinder for non W-Ta presets or if explicitly requested.
  const bool useSectionalWTa = (config_.target.type == "W-Ta" && config_.geometry.simpleCylinder);

  if (!useSectionalWTa) {
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

  // Realistic W-Ta sectional assembly geometry.
  auto* wMat = nist->FindOrBuildMaterial("G4_W");
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

  // Integrated 12 cm target node container, filled with water (explicit choice).
  auto* assemblyS = new G4Box("TargetAssembly", taHalfX, taHalfY, 0.5 * assemblyT);
  auto* assemblyLV = new G4LogicalVolume(assemblyS, waterMat, "TargetAssemblyLV");
  new G4PVPlacement(nullptr, {}, assemblyLV, "TargetAssembly", worldLV, false, 0, true);
  ConfigureVis(assemblyLV, G4Colour(0.7, 0.9, 0.95, 0.15));

  double stackT = 8.0 * waterGap;
  for (double tmm : wPlateTsMM) {
    stackT += tmm * mm + 2.0 * ti + 2.0 * ta;
  }

  G4cout << "[geom] W-Ta sectional stack total thickness = " << stackT / mm
         << " mm; assembly thickness = " << assemblyT / mm << " mm" << G4endl;

  if (stackT > assemblyT) {
    G4Exception("DetectorConstruction::Construct", "Geom001", FatalException,
                "Sectional stack thickness exceeds assembly_thickness_mm");
  }

  double zCursor = -0.5 * stackT;

  for (size_t i = 0; i < wPlateTsMM.size(); ++i) {
    // Water layer before plate i.
    {
      const double hz = 0.5 * waterGap;
      auto* s = new G4Box("TargetWaterGapShape", taHalfX, taHalfY, hz);
      auto* lv = new G4LogicalVolume(s, waterMat, "TargetWaterGap");
      new G4PVPlacement(nullptr, {0, 0, zCursor + hz}, lv,
                        G4String("TargetWaterGap_") + std::to_string(i), assemblyLV, false, static_cast<G4int>(i), true);
      ConfigureVis(lv, G4Colour(0.6, 0.85, 1.0, 0.35));
      zCursor += waterGap;
    }

    const double wT = wPlateTsMM[i] * mm;
    const double taTot = wT + 2.0 * ti + 2.0 * ta;
    const double tiTot = wT + 2.0 * ti;

    // Ta cladding (outer envelope).
    auto* taS = new G4Box("TargetCoatingShape", taHalfX, taHalfY, 0.5 * taTot);
    auto* taLV = new G4LogicalVolume(taS, taMat, "TargetCoating");
    new G4PVPlacement(nullptr, {0, 0, zCursor + 0.5 * taTot}, taLV,
                      G4String("TargetCoating_") + std::to_string(i), assemblyLV, false, static_cast<G4int>(i), true);
    ConfigureVis(taLV, G4Colour(0.95, 0.5, 0.2));
    if (i == 0) coatingLV_ = taLV;

    // Ti buffer nested inside Ta.
    auto* tiS = new G4Box("TargetBufferTiShape", plateHalfXY + ti, plateHalfXY + ti, 0.5 * tiTot);
    auto* tiLV = new G4LogicalVolume(tiS, tiMat, "TargetBufferTi");
    new G4PVPlacement(nullptr, {}, tiLV,
                      G4String("TargetBufferTi_") + std::to_string(i), taLV, false, static_cast<G4int>(i), true);
    ConfigureVis(tiLV, G4Colour(0.7, 0.7, 0.75));

    // W core nested inside Ti.
    auto* wS = new G4Box("TargetSubstrateShape", plateHalfXY, plateHalfXY, 0.5 * wT);
    auto* wLV = new G4LogicalVolume(wS, wMat, "TargetSubstrate");
    new G4PVPlacement(nullptr, {}, wLV,
                      G4String("TargetSubstrate_") + std::to_string(i), tiLV, false, static_cast<G4int>(i), true);
    ConfigureVis(wLV, G4Colour(0.2, 0.4, 0.9));
    if (i == 0) substrateLV_ = wLV;

    zCursor += taTot;
  }

  // Final 8th water layer after plate 7.
  {
    const double hz = 0.5 * waterGap;
    auto* s = new G4Box("TargetWaterGapShapeFinal", taHalfX, taHalfY, hz);
    auto* lv = new G4LogicalVolume(s, waterMat, "TargetWaterGap");
    new G4PVPlacement(nullptr, {0, 0, zCursor + hz}, lv,
                      "TargetWaterGap_7", assemblyLV, false, 7, true);
    ConfigureVis(lv, G4Colour(0.6, 0.85, 1.0, 0.35));
  }

  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  return worldPV;
}
