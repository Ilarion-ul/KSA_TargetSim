#pragma once

#include "Config.h"

#include <G4LogicalVolume.hh>
#include <G4VUserDetectorConstruction.hh>

class G4VPhysicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction {
 public:
  explicit DetectorConstruction(AppConfig config);
  ~DetectorConstruction() override = default;

  G4VPhysicalVolume* Construct() override;

  G4LogicalVolume* GetSubstrateLV() const { return substrateLV_; }
  G4LogicalVolume* GetCoatingLV() const { return coatingLV_; }

 private:
  AppConfig config_;
  G4LogicalVolume* substrateLV_{nullptr};
  G4LogicalVolume* coatingLV_{nullptr};

  // TODO: add realistic target shell/cooling/mechanical assemblies.
};
