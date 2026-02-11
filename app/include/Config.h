#pragma once

#include <array>
#include <optional>
#include <string>
#include <vector>

// NOTE: Units are stored in human-friendly values (MeV/mm/mrad).
// Conversion to Geant4 units is done at use-sites in simulation code.

struct BeamConfig {
  // Central energy and energy spread.
  double energy_MeV{100.0};
  double energy_sigma_rel_1sigma{0.01};
  double energy_uniform_rel_halfspan{0.04};
  std::string energy_spread_model{"gauss"}; // gauss | uniform

  // Position and nominal direction.
  std::array<double, 3> position_mm{0.0, 0.0, -5.0}; // alias: pos_mm
  std::array<double, 3> direction{0.0, 0.0, 1.0};    // alias: dir

  // Transverse beam spot.
  double sigmaX_mm{1.0}; // alias: sigma_x_mm
  double sigmaY_mm{1.0}; // alias: sigma_y_mm

  // Angular spread (direct mode).
  double sigma_theta_x_mrad{1.0};
  double sigma_theta_y_mrad{1.0};

  // Legacy angle key (kept for compatibility).
  double divergence_mrad{1.0};

  // Optional emit-based approximation.
  double norm_emit_m_rad{5e-7};
  bool use_emit_model{false};
  std::string emit_sigma_theta_from{"eps_over_sigma"};

  // Beam operation mode metadata.
  std::string mode{"cw"}; // cw | pulsed
  double pulse_width_us{2.7};
  double rep_rate_Hz{625.0};
  double I_pulse_A{0.6};
  double I_avg_A{0.001};
  double beam_power_kW{100.0};

  // Defect model.
  std::string defect_mode{"ideal"}; // ideal / offset / tilt / halo
  std::array<double, 3> offset_mm{0.0, 0.0, 0.0}; // defect.offset_mm
  std::array<double, 2> tilt_mrad{0.0, 0.0};      // defect.tilt_mrad
  double halo_fraction{0.0};
  double halo_sigma_mm{3.0};
  double halo_sigma_scale{5.0};
};

struct TargetConfig {
  std::string type{"W-Ta"};

  // Legacy/simple geometry fields.
  double substrate_thickness_mm{10.0};
  double coating_thickness_mm{0.5};
  double radius_mm{10.0};

  // Realistic sectional W-Ta geometry fields.
  double plate_xy_mm{65.8};
  std::vector<double> plate_thicknesses_mm{2.5, 2.5, 2.5, 3.5, 3.5, 5.5, 9.5};
  double water_gap_mm{2.0};
  double clad_ta_mm{0.25};
  double buffer_ti_mm{0.04};
  double assembly_thickness_mm{120.0};

  // U-Mo target and construction materials.
  double clad_thickness_front_mm{0.95}; // plates 1..4
  double clad_thickness_rest_mm{0.70};  // plates 5..12
  double gap_inout_mm{1.0};
  double gap_mid_mm{1.75};
  double housing_inner_xy_mm{66.0};
  double housing_wall_mm{2.0};
  double entrance_window_mm{2.0};
  double helium_chamber_len_mm{237.0};
  double u7mo_density_g_cm3{17.0}; // TODO: validate exact tech-specific density.
  std::string fill_medium_in_target{"water"};
  std::string w_substrate_material{"pure_W"}; // pure_W | W_Fe_Ni

  // Optional explicit U-Mo inter-plate gaps (size = plate_count - 1).
  // If empty, a rule-based fallback is used (front/rear split or legacy gap_mid).
  std::vector<double> inter_plate_gaps_mm{};
  double gap_front_mm{3.0};
  double gap_rear_mm{1.75};
  int gap_split_index{3}; // number of first inter-plate gaps using gap_front_mm

  double temperature_K{300.0};
};

struct RunConfig {
  int nEvents{1000};
  int nThreads{0};
  std::string outputRootFile{"run.root"};
  std::string outputDir{"results"};
  bool enableVis{false};
  bool enableEventTree{false}; // TODO: optional per-event ROOT tree switch.
  double irradiation_time_s{1.0};
  bool enableSwellingOutput{true};
};

struct PhysicsConfig {
  std::string physicsListName{"QGSP_BIC_HPT"};
  bool enablePhotonuclear{true}; // placeholder for future toggles.
  std::optional<double> cut_mm{0.1};
};

struct GeometryConfig {
  bool simpleCylinder{true};
  double worldMargin_mm{50.0};

  // Assembly/beamline split (used by U-Mo geometry branch).
  double total_assembly_len_mm{2620.0};
  double beamline_vacuum_len_mm{2210.0};
  double target_region_extra_clearance_mm{0.0};

  // Optional centering/alignment pin ("finger").
  bool enable_alignment_pin{false};
  std::array<double, 3> alignment_pin_pos_mm{0.0, 0.0, 0.0};
  std::array<double, 3> alignment_pin_size_mm{2.0, 2.0, 10.0};
  std::string alignment_pin_material{"G4_Al"};
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
