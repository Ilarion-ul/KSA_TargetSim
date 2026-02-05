#include "DetectorConstruction.h"

#include <G4Box.hh>
#include <G4Colour.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4SystemOfUnits.hh>
#include <G4Tubs.hh>
#include <G4VisAttributes.hh>

DetectorConstruction::DetectorConstruction(AppConfig config) : config_(std::move(config)) {}

G4VPhysicalVolume* DetectorConstruction::Construct() {
  auto* nist = G4NistManager::Instance();

  auto* worldMat = nist->FindOrBuildMaterial("G4_Galactic");
  auto* substrateMat = nist->FindOrBuildMaterial(config_.target.type == "U-Al" ? "G4_U" : "G4_W");
  auto* coatingMat = nist->FindOrBuildMaterial(config_.target.type == "U-Al" ? "G4_Al" : "G4_Ta");

  // Coordinate convention: beam travels generally along +Z.
  // Target is two cylinders placed one after another along Z:
  // [substrate] then [coating] on the +Z end-face of substrate.
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

  auto* subVis = new G4VisAttributes(G4Colour(0.2, 0.4, 0.9));
  subVis->SetForceSolid(true);
  substrateLV_->SetVisAttributes(subVis);

  auto* coatVis = new G4VisAttributes(G4Colour(0.9, 0.4, 0.2));
  coatVis->SetForceSolid(true);
  coatingLV_->SetVisAttributes(coatVis);

  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  return worldPV;
}
