#pragma once

#include "Config.h"

#include <G4VUserPrimaryGeneratorAction.hh>

class G4Event;
class G4ParticleGun;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  explicit PrimaryGeneratorAction(const AppConfig& config);
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* event) override;

 private:
  AppConfig config_;
  G4ParticleGun* gun_{nullptr};
};
