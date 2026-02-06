#include "PrimaryGeneratorAction.h"

#include <G4Event.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4ios.hh>
#include <Randomize.hh>

#include <algorithm>
#include <cmath>

namespace {
constexpr double kElectronMassMeV = 0.51099895;
}

PrimaryGeneratorAction::PrimaryGeneratorAction(const AppConfig& config) : config_(config) {
  gun_ = new G4ParticleGun(1);
  G4ParticleDefinition* electron = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  gun_->SetParticleDefinition(electron);

  const auto& b = config_.beam;
  const double gamma = 1.0 + b.energy_MeV / kElectronMassMeV;
  const double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
  const double duty = (b.pulse_width_us * 1e-6) * b.rep_rate_Hz;
  const double iAvgCalc = b.I_pulse_A * duty;

  G4cout << "[beam] E0 = " << b.energy_MeV << " MeV, spread_model=" << b.energy_spread_model
         << ", sigma_rel=" << b.energy_sigma_rel_1sigma
         << ", uniform_halfspan=" << b.energy_uniform_rel_halfspan << G4endl;
  if (b.use_emit_model) {
    const double epsGeom = b.norm_emit_m_rad / (beta * gamma);
    const double sigmaX = std::max(1e-12, b.sigmaX_mm * 1e-3);
    const double sigmaY = std::max(1e-12, b.sigmaY_mm * 1e-3);
    G4cout << "[beam] emit_model=ON, eps_n=" << b.norm_emit_m_rad << " m*rad, beta*gamma=" << beta * gamma
           << ", approx sigma_theta_x=" << (epsGeom / sigmaX) * 1e3 << " mrad"
           << ", sigma_theta_y=" << (epsGeom / sigmaY) * 1e3 << " mrad" << G4endl;
  } else {
    G4cout << "[beam] sigma_x/y(mm)=(" << b.sigmaX_mm << "," << b.sigmaY_mm << ")"
           << ", sigma_theta_x/y(mrad)=(" << b.sigma_theta_x_mrad << "," << b.sigma_theta_y_mrad << ")"
           << G4endl;
  }

  G4cout << "[beam] mode=" << b.mode << ", pulse_width_us=" << b.pulse_width_us
         << ", rep_rate_Hz=" << b.rep_rate_Hz << ", I_pulse_A=" << b.I_pulse_A
         << ", I_avg_A(input)=" << b.I_avg_A << ", I_avg_A(calc)=" << iAvgCalc
         << ", duty=" << duty << ", P_avg_kW(meta)=" << b.beam_power_kW
         << "  [per-event results remain normalized per primary electron]" << G4endl;

  // TODO: if time structure is required, set G4PrimaryVertex time in GeneratePrimaries.
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete gun_;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  const auto& b = config_.beam;

  // Energy spread model.
  double rel = 0.0;
  if (b.energy_spread_model == "uniform") {
    rel = (2.0 * G4UniformRand() - 1.0) * b.energy_uniform_rel_halfspan;
  } else {
    rel = G4RandGauss::shoot(0.0, b.energy_sigma_rel_1sigma);
  }
  double energyMeV = b.energy_MeV * (1.0 + rel);
  energyMeV = std::max(1e-6, energyMeV);

  double sigmaX = b.sigmaX_mm;
  double sigmaY = b.sigmaY_mm;

  if (b.defect_mode == "halo" && G4UniformRand() < b.halo_fraction) {
    const double scale = std::max(1.0, b.halo_sigma_scale);
    sigmaX *= scale;
    sigmaY *= scale;
  }

  double x = G4RandGauss::shoot(0.0, sigmaX) + b.position_mm[0];
  double y = G4RandGauss::shoot(0.0, sigmaY) + b.position_mm[1];
  double z = b.position_mm[2];

  if (b.defect_mode == "offset") {
    x += b.offset_mm[0];
    y += b.offset_mm[1];
  }

  double dx = b.direction[0];
  double dy = b.direction[1];
  double dz = b.direction[2];

  double sigmaThetaX = b.sigma_theta_x_mrad * 1e-3;
  double sigmaThetaY = b.sigma_theta_y_mrad * 1e-3;

  if (b.use_emit_model && b.emit_sigma_theta_from == "eps_over_sigma") {
    const double gamma = 1.0 + b.energy_MeV / kElectronMassMeV;
    const double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
    const double epsGeom = b.norm_emit_m_rad / (beta * gamma);
    const double sigmaX_m = std::max(1e-12, b.sigmaX_mm * 1e-3);
    const double sigmaY_m = std::max(1e-12, b.sigmaY_mm * 1e-3);
    // Approximation without full Twiss transport.
    sigmaThetaX = epsGeom / sigmaX_m;
    sigmaThetaY = epsGeom / sigmaY_m;
  }

  // Backward compatibility: if explicit theta values are absent, use divergence.
  if (sigmaThetaX <= 0.0 && sigmaThetaY <= 0.0) {
    sigmaThetaX = b.divergence_mrad * 1e-3;
    sigmaThetaY = b.divergence_mrad * 1e-3;
  }

  dx += G4RandGauss::shoot(0.0, sigmaThetaX);
  dy += G4RandGauss::shoot(0.0, sigmaThetaY);

  if (b.defect_mode == "tilt") {
    dx += b.tilt_mrad[0] * 1e-3;
    dy += b.tilt_mrad[1] * 1e-3;
  }

  const double norm = std::sqrt(dx * dx + dy * dy + dz * dz);
  dx /= norm;
  dy /= norm;
  dz /= norm;

  gun_->SetParticleEnergy(energyMeV * MeV);
  gun_->SetParticlePosition({x * mm, y * mm, z * mm});
  gun_->SetParticleMomentumDirection({dx, dy, dz});
  gun_->GeneratePrimaryVertex(event);
}
