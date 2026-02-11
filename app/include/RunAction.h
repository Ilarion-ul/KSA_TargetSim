#pragma once

#include "Config.h"

#include <G4UserRunAction.hh>

#include <mutex>
#include <string>
#include <vector>

class G4Run;

class RunAction : public G4UserRunAction {
 public:
  struct NeutronSurfaceHit {
    int event_id{0};
    double En_MeV{0.0};
    double x_mm{0.0};
    double y_mm{0.0};
    double z_mm{0.0};
    double cosTheta{0.0};
    double weight{1.0};
    double time_ns{0.0};
    int surface_id{0};
  };

  struct PhotonSurfaceHit {
    int event_id{0};
    double E_MeV{0.0};
    double x_mm{0.0};
    double y_mm{0.0};
    double z_mm{0.0};
    double cosTheta{0.0};
    double weight{1.0};
    double time_ns{0.0};
    int surface_id{0};
  };

  struct HeatmapBounds {
    double xMinMm;
    double xMaxMm;
    double yMinMm;
    double yMaxMm;
    double zMinMm;
    double zMaxMm;
  };

  explicit RunAction(const AppConfig& config);
  ~RunAction() override;

  void BeginOfRunAction(const G4Run* run) override;
  void EndOfRunAction(const G4Run* run) override;

  void AccumulateEvent(double edepSubstrate,
                       double edepCoating,
                       int nGamma,
                       int nNeutron,
                       int nNeutronExit,
                       int nNeutronModelExit,
                       int eventId,
                       const std::vector<double>& plateEdep,
                       const std::vector<double>& plateNeutronTrackLen,
                       const std::vector<double>& plateNeutronHeatmap,
                       const std::vector<double>& edep3d,
                       const std::vector<NeutronSurfaceHit>& neutronSurfaceHits,
                       const std::vector<PhotonSurfaceHit>& photonSurfaceHits,
                       const std::vector<double>& plateNiel,
                       const std::vector<double>& plateGasH,
                       const std::vector<double>& plateGasHe);
  size_t PlateCount() const { return plateCount_; }
  int PlateHeatmapBinsX() const { return plateHeatmapBinsX_; }
  int PlateHeatmapBinsY() const { return plateHeatmapBinsY_; }
  int EdepBinsX() const { return edepBinsX_; }
  int EdepBinsY() const { return edepBinsY_; }
  int EdepBinsZ() const { return edepBinsZ_; }
  HeatmapBounds PlateHeatmapBounds() const { return plateHeatmapBounds_; }
  HeatmapBounds EdepBounds() const { return edepBounds_; }

 private:
  AppConfig config_;
  std::mutex mutex_;
  size_t plateCount_{0};
  int plateHeatmapBinsX_{0};
  int plateHeatmapBinsY_{0};
  int edepBinsX_{0};
  int edepBinsY_{0};
  int edepBinsZ_{0};
  HeatmapBounds plateHeatmapBounds_{};
  HeatmapBounds edepBounds_{};

  double totalEdepSubstrate_{0.0};
  double totalEdepCoating_{0.0};
  long long totalGamma_{0};
  long long totalNeutron_{0};
  long long totalNeutronExit_{0};
  std::vector<double> totalPlateEdep_;
  std::vector<double> totalPlateNeutronTrackLen_;
  std::vector<double> totalPlateNeutronHeatmap_;
  std::vector<double> totalEdep3d_;
  std::vector<double> totalPlateNiel_;
  std::vector<double> totalPlateGasH_;
  std::vector<double> totalPlateGasHe_;

#ifdef KSA_USE_ROOT
  class TFile* rootFile_{nullptr};
  class TTree* runTree_{nullptr};
#endif
};
