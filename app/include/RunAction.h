#pragma once

#include "Config.h"

#include <G4UserRunAction.hh>

#include <mutex>
#include <string>
#include <vector>

class G4Run;

class RunAction : public G4UserRunAction {
 public:
  explicit RunAction(const AppConfig& config);
  ~RunAction() override;

  void BeginOfRunAction(const G4Run* run) override;
  void EndOfRunAction(const G4Run* run) override;

  void AccumulateEvent(double edepSubstrate,
                       double edepCoating,
                       int nGamma,
                       int nNeutron,
                       int nNeutronExit,
                       const std::vector<double>& plateEdep,
                       const std::vector<double>& plateNeutronTrackLen);
  size_t PlateCount() const { return plateCount_; }

 private:
  AppConfig config_;
  std::mutex mutex_;
  size_t plateCount_{0};

  double totalEdepSubstrate_{0.0};
  double totalEdepCoating_{0.0};
  long long totalGamma_{0};
  long long totalNeutron_{0};
  long long totalNeutronExit_{0};
  std::vector<double> totalPlateEdep_;
  std::vector<double> totalPlateNeutronTrackLen_;

#ifdef KSA_USE_ROOT
  class TFile* rootFile_{nullptr};
  class TTree* runTree_{nullptr};
#endif
};
