#include "PrimaryGeneratorAction.h"

#include <G4Event.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>

#include <cmath>

PrimaryGeneratorAction::PrimaryGeneratorAction(const AppConfig& config) : config_(config) {
  gun_ = new G4ParticleGun(1);
  G4ParticleDefinition* electron = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  gun_->SetParticleDefinition(electron);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete gun_;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  // TODO: expose deterministic seed control via JSON/CLI.
  const auto& b = config_.beam;
  const double div = b.divergence_mrad * 1e-3;

  double sigmaX = b.sigmaX_mm;
  double sigmaY = b.sigmaY_mm;

  if (b.defect_mode == "halo" && G4UniformRand() < b.halo_fraction) {
    sigmaX = b.halo_sigma_mm;
    sigmaY = b.halo_sigma_mm;
  }

  double x = G4RandGauss::shoot(0.0, sigmaX) + b.position_mm[0];
  double y = G4RandGauss::shoot(0.0, sigmaY) + b.position_mm[1];
  double z = b.position_mm[2];

  if (b.defect_mode == "offset") {
    x += b.offset_mm[0];
    y += b.offset_mm[1];
    z += b.offset_mm[2];
  }

  double dx = b.direction[0];
  double dy = b.direction[1];
  double dz = b.direction[2];

  dx += G4RandGauss::shoot(0.0, div);
  dy += G4RandGauss::shoot(0.0, div);

  if (b.defect_mode == "tilt") {
    dx += b.tilt_mrad[0] * 1e-3;
    dy += b.tilt_mrad[1] * 1e-3;
  }

  const double norm = std::sqrt(dx * dx + dy * dy + dz * dz);
  dx /= norm;
  dy /= norm;
  dz /= norm;

  gun_->SetParticleEnergy(b.energy_MeV * MeV);
  gun_->SetParticlePosition({x * mm, y * mm, z * mm});
  gun_->SetParticleMomentumDirection({dx, dy, dz});
  gun_->GeneratePrimaryVertex(event);
}
