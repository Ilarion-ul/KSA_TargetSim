#include "Config.h"

#include <algorithm>
#include <fstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

namespace {
using json = nlohmann::json;

template <typename T>
void get_if_exists(const json& j, const char* key, T& out) {
  if (j.contains(key)) {
    out = j.at(key).get<T>();
  }
}

void validate_required_sections(const json& j) {
  if (!j.contains("beam") || !j.contains("target") || !j.contains("run")) {
    throw std::runtime_error("Config validation failed: required sections beam/target/run are missing");
  }
}
} // namespace

AppConfig LoadConfig(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Failed to open config file: " + path);
  }

  json j;
  in >> j;
  validate_required_sections(j);

  AppConfig cfg;

  const auto& jb = j.at("beam");
  get_if_exists(jb, "energy_MeV", cfg.beam.energy_MeV);
  get_if_exists(jb, "energy_sigma_rel_1sigma", cfg.beam.energy_sigma_rel_1sigma);
  get_if_exists(jb, "energy_uniform_rel_halfspan", cfg.beam.energy_uniform_rel_halfspan);
  get_if_exists(jb, "energy_spread_model", cfg.beam.energy_spread_model);

  get_if_exists(jb, "position_mm", cfg.beam.position_mm);
  get_if_exists(jb, "pos_mm", cfg.beam.position_mm);
  get_if_exists(jb, "direction", cfg.beam.direction);
  get_if_exists(jb, "dir", cfg.beam.direction);

  get_if_exists(jb, "sigmaX_mm", cfg.beam.sigmaX_mm);
  get_if_exists(jb, "sigmaY_mm", cfg.beam.sigmaY_mm);
  get_if_exists(jb, "sigma_x_mm", cfg.beam.sigmaX_mm);
  get_if_exists(jb, "sigma_y_mm", cfg.beam.sigmaY_mm);

  get_if_exists(jb, "sigma_theta_x_mrad", cfg.beam.sigma_theta_x_mrad);
  get_if_exists(jb, "sigma_theta_y_mrad", cfg.beam.sigma_theta_y_mrad);
  get_if_exists(jb, "divergence_mrad", cfg.beam.divergence_mrad);

  get_if_exists(jb, "norm_emit_m_rad", cfg.beam.norm_emit_m_rad);
  get_if_exists(jb, "use_emit_model", cfg.beam.use_emit_model);
  get_if_exists(jb, "emit_sigma_theta_from", cfg.beam.emit_sigma_theta_from);

  get_if_exists(jb, "mode", cfg.beam.mode);
  get_if_exists(jb, "pulse_width_us", cfg.beam.pulse_width_us);
  get_if_exists(jb, "rep_rate_Hz", cfg.beam.rep_rate_Hz);
  get_if_exists(jb, "I_pulse_A", cfg.beam.I_pulse_A);
  get_if_exists(jb, "I_avg_A", cfg.beam.I_avg_A);
  get_if_exists(jb, "beam_power_kW", cfg.beam.beam_power_kW);

  get_if_exists(jb, "defect_mode", cfg.beam.defect_mode);
  get_if_exists(jb, "offset_mm", cfg.beam.offset_mm);
  get_if_exists(jb, "tilt_mrad", cfg.beam.tilt_mrad);
  get_if_exists(jb, "halo_fraction", cfg.beam.halo_fraction);
  get_if_exists(jb, "halo_sigma_mm", cfg.beam.halo_sigma_mm);
  get_if_exists(jb, "halo_sigma_scale", cfg.beam.halo_sigma_scale);

  if (jb.contains("defect")) {
    const auto& jd = jb.at("defect");
    if (jd.contains("offset_mm")) {
      std::vector<double> tmp = jd.at("offset_mm").get<std::vector<double>>();
      if (tmp.size() == 2) {
        cfg.beam.offset_mm = {tmp[0], tmp[1], 0.0};
      } else if (tmp.size() == 3) {
        cfg.beam.offset_mm = {tmp[0], tmp[1], tmp[2]};
      } else {
        throw std::runtime_error("Config validation failed: beam.defect.offset_mm must have 2 or 3 values");
      }
    }
    get_if_exists(jd, "tilt_mrad", cfg.beam.tilt_mrad);
    get_if_exists(jd, "halo_fraction", cfg.beam.halo_fraction);
    get_if_exists(jd, "halo_sigma_scale", cfg.beam.halo_sigma_scale);
  }

  if (cfg.beam.mode != "cw" && cfg.beam.mode != "pulsed") {
    throw std::runtime_error("Config validation failed: beam.mode must be 'cw' or 'pulsed'");
  }
  if (cfg.beam.energy_spread_model != "gauss" && cfg.beam.energy_spread_model != "uniform") {
    throw std::runtime_error("Config validation failed: beam.energy_spread_model must be 'gauss' or 'uniform'");
  }

  const auto& jt = j.at("target");
  get_if_exists(jt, "type", cfg.target.type);
  if (cfg.target.type != "W-Ta" && cfg.target.type != "U-Al" && cfg.target.type != "U-Mo") {
    throw std::runtime_error("Config validation failed: target.type must be 'W-Ta', 'U-Al', or 'U-Mo'");
  }

  get_if_exists(jt, "substrate_thickness_mm", cfg.target.substrate_thickness_mm);
  get_if_exists(jt, "coating_thickness_mm", cfg.target.coating_thickness_mm);
  get_if_exists(jt, "radius_mm", cfg.target.radius_mm);
  get_if_exists(jt, "temperature_K", cfg.target.temperature_K);

  get_if_exists(jt, "plate_xy_mm", cfg.target.plate_xy_mm);
  get_if_exists(jt, "plate_thicknesses_mm", cfg.target.plate_thicknesses_mm);
  get_if_exists(jt, "water_gap_mm", cfg.target.water_gap_mm);
  get_if_exists(jt, "clad_ta_mm", cfg.target.clad_ta_mm);
  get_if_exists(jt, "buffer_ti_mm", cfg.target.buffer_ti_mm);
  get_if_exists(jt, "assembly_thickness_mm", cfg.target.assembly_thickness_mm);

  get_if_exists(jt, "clad_thickness_front_mm", cfg.target.clad_thickness_front_mm);
  get_if_exists(jt, "clad_thickness_rest_mm", cfg.target.clad_thickness_rest_mm);
  get_if_exists(jt, "gap_inout_mm", cfg.target.gap_inout_mm);
  get_if_exists(jt, "gap_mid_mm", cfg.target.gap_mid_mm);
  get_if_exists(jt, "housing_inner_xy_mm", cfg.target.housing_inner_xy_mm);
  get_if_exists(jt, "housing_wall_mm", cfg.target.housing_wall_mm);
  get_if_exists(jt, "entrance_window_mm", cfg.target.entrance_window_mm);
  get_if_exists(jt, "helium_chamber_len_mm", cfg.target.helium_chamber_len_mm);
  get_if_exists(jt, "u7mo_density_g_cm3", cfg.target.u7mo_density_g_cm3);
  get_if_exists(jt, "fill_medium_in_target", cfg.target.fill_medium_in_target);
  get_if_exists(jt, "w_substrate_material", cfg.target.w_substrate_material);
  get_if_exists(jt, "inter_plate_gaps_mm", cfg.target.inter_plate_gaps_mm);
  get_if_exists(jt, "gap_front_mm", cfg.target.gap_front_mm);
  get_if_exists(jt, "gap_rear_mm", cfg.target.gap_rear_mm);
  get_if_exists(jt, "gap_split_index", cfg.target.gap_split_index);

  for (double v : cfg.target.plate_thicknesses_mm) {
    if (v <= 0.0) {
      throw std::runtime_error("Config validation failed: each target.plate_thicknesses_mm element must be > 0");
    }
  }

  if (cfg.target.type == "W-Ta") {
    if (cfg.target.plate_thicknesses_mm.size() != 7) {
      throw std::runtime_error("Config validation failed: target.plate_thicknesses_mm must contain 7 values for W-Ta");
    }
    if (cfg.target.buffer_ti_mm < 0.03 || cfg.target.buffer_ti_mm > 0.06) {
      throw std::runtime_error("Config validation failed: target.buffer_ti_mm must be in [0.03, 0.06] mm");
    }
    if (cfg.target.clad_ta_mm < 0.25 || cfg.target.clad_ta_mm > 0.27) {
      throw std::runtime_error("Config validation failed: target.clad_ta_mm must be in [0.25, 0.27] mm");
    }
    if (cfg.target.water_gap_mm != 2.0 && cfg.target.water_gap_mm != 1.75) {
      throw std::runtime_error("Config validation failed: target.water_gap_mm must be 2.0 or 1.75 mm for W-Ta presets");
    }
    if (cfg.target.w_substrate_material != "pure_W" && cfg.target.w_substrate_material != "W_Fe_Ni") {
      throw std::runtime_error("Config validation failed: target.w_substrate_material must be 'pure_W' or 'W_Fe_Ni'");
    }
  }

  if (cfg.target.type == "U-Mo") {
    const auto plateCount = cfg.target.plate_thicknesses_mm.size();
    if (plateCount != 6 && plateCount != 12) {
      throw std::runtime_error("Config validation failed: target.plate_thicknesses_mm must contain 6 (new spec) or 12 (legacy) values for U-Mo");
    }
    if (cfg.target.fill_medium_in_target != "water" && cfg.target.fill_medium_in_target != "vacuum") {
      throw std::runtime_error("Config validation failed: target.fill_medium_in_target must be 'water' or 'vacuum'");
    }
    if (!cfg.target.inter_plate_gaps_mm.empty() && cfg.target.inter_plate_gaps_mm.size() != plateCount - 1) {
      throw std::runtime_error("Config validation failed: target.inter_plate_gaps_mm size must be plate_count-1");
    }
    for (double g : cfg.target.inter_plate_gaps_mm) {
      if (g <= 0.0) {
        throw std::runtime_error("Config validation failed: target.inter_plate_gaps_mm values must be > 0");
      }
    }
    if (cfg.target.gap_front_mm <= 0.0 || cfg.target.gap_rear_mm <= 0.0) {
      throw std::runtime_error("Config validation failed: target.gap_front_mm and target.gap_rear_mm must be > 0");
    }
    if (cfg.target.gap_split_index < 0 || static_cast<size_t>(cfg.target.gap_split_index) > plateCount - 1) {
      throw std::runtime_error("Config validation failed: target.gap_split_index must be in [0, plate_count-1]");
    }
    const double maxClad = std::max(cfg.target.clad_thickness_front_mm, cfg.target.clad_thickness_rest_mm);
    if (cfg.target.plate_xy_mm + 2.0 * maxClad > cfg.target.housing_inner_xy_mm) {
      throw std::runtime_error("Config validation failed: U-Mo plate+cladding does not fit housing_inner_xy_mm");
    }
  }

  const auto& jr = j.at("run");
  get_if_exists(jr, "nEvents", cfg.run.nEvents);
  get_if_exists(jr, "nThreads", cfg.run.nThreads);
  get_if_exists(jr, "outputRootFile", cfg.run.outputRootFile);
  get_if_exists(jr, "outputDir", cfg.run.outputDir);
  get_if_exists(jr, "enableVis", cfg.run.enableVis);
  get_if_exists(jr, "enableEventTree", cfg.run.enableEventTree);
  get_if_exists(jr, "irradiation_time_s", cfg.run.irradiation_time_s);
  get_if_exists(jr, "enableSwellingOutput", cfg.run.enableSwellingOutput);
  if (cfg.run.nThreads < 0) {
    cfg.run.nThreads = 0;
  }
  if (cfg.run.outputDir.empty()) {
    cfg.run.outputDir = "results";
  }
  if (cfg.run.irradiation_time_s <= 0.0) {
    cfg.run.irradiation_time_s = 1.0;
  }

  if (j.contains("physics")) {
    const auto& jp = j.at("physics");
    get_if_exists(jp, "physicsListName", cfg.physics.physicsListName);
    get_if_exists(jp, "enablePhotonuclear", cfg.physics.enablePhotonuclear);
    if (jp.contains("cut_mm") && !jp.at("cut_mm").is_null()) {
      cfg.physics.cut_mm = jp.at("cut_mm").get<double>();
    }
  }

  if (j.contains("geometry")) {
    const auto& jg = j.at("geometry");
    get_if_exists(jg, "simpleCylinder", cfg.geometry.simpleCylinder);
    get_if_exists(jg, "worldMargin_mm", cfg.geometry.worldMargin_mm);
    get_if_exists(jg, "total_assembly_len_mm", cfg.geometry.total_assembly_len_mm);
    get_if_exists(jg, "beamline_vacuum_len_mm", cfg.geometry.beamline_vacuum_len_mm);
    get_if_exists(jg, "target_region_extra_clearance_mm", cfg.geometry.target_region_extra_clearance_mm);
    get_if_exists(jg, "enable_alignment_pin", cfg.geometry.enable_alignment_pin);
    get_if_exists(jg, "alignment_pin_pos_mm", cfg.geometry.alignment_pin_pos_mm);
    get_if_exists(jg, "alignment_pin_size_mm", cfg.geometry.alignment_pin_size_mm);
    get_if_exists(jg, "alignment_pin_material", cfg.geometry.alignment_pin_material);
  }

  if (cfg.geometry.alignment_pin_size_mm[0] <= 0.0 || cfg.geometry.alignment_pin_size_mm[1] <= 0.0 ||
      cfg.geometry.alignment_pin_size_mm[2] <= 0.0) {
    throw std::runtime_error("Config validation failed: geometry.alignment_pin_size_mm must be positive in all components");
  }

  if (j.contains("defects")) {
    const auto& jd = j.at("defects");
    get_if_exists(jd, "mode", cfg.defects.mode);
  } else {
    cfg.defects.mode = cfg.beam.defect_mode;
  }

  return cfg;
}
