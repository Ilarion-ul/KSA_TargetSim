#pragma once

#include <array>
#include <optional>
#include <string>

// NOTE: Units are stored in human-friendly values (MeV/mm/mrad).
// Conversion to Geant4 units is done at use-sites in simulation code.

struct BeamConfig {
  double energy_MeV{100.0};
  double sigmaX_mm{1.0};
  double sigmaY_mm{1.0};
  double divergence_mrad{1.0};
  std::array<double, 3> position_mm{0.0, 0.0, -5.0};
  std::array<double, 3> direction{0.0, 0.0, 1.0};

  std::string defect_mode{"ideal"}; // ideal / offset / tilt / halo
  std::array<double, 3> offset_mm{0.0, 0.0, 0.0};
  std::array<double, 2> tilt_mrad{0.0, 0.0};
  double halo_fraction{0.0};
  double halo_sigma_mm{3.0};
};

struct TargetConfig {
  std::string type{"W-Ta"};
  double substrate_thickness_mm{10.0};
  double coating_thickness_mm{0.5};
  double radius_mm{10.0};
  double temperature_K{300.0};
};

struct RunConfig {
  int nEvents{1000};
  int nThreads{0};
  std::string outputRootFile{"run.root"};
  std::string outputDir{"results"};
  bool enableVis{false};
  bool enableEventTree{false}; // TODO: optional per-event ROOT tree switch.
};

struct PhysicsConfig {
  std::string physicsListName{"QGSP_BIC_HPT"};
  bool enablePhotonuclear{true}; // placeholder for future toggles.
  std::optional<double> cut_mm{0.1};
};

struct GeometryConfig {
  bool simpleCylinder{true};
  double worldMargin_mm{50.0};
};

struct DefectConfig {
  // TODO: extended defects block for future detailed models.
  std::string mode{"ideal"};
};

struct AppConfig {
  BeamConfig beam;
  TargetConfig target;
  RunConfig run;
  PhysicsConfig physics;
  GeometryConfig geometry;
  DefectConfig defects;
};

AppConfig LoadConfig(const std::string& path);
