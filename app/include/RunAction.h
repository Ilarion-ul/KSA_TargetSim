#pragma once

#include "Config.h"

#include <G4UserRunAction.hh>

#include <mutex>
#include <string>

class G4Run;

class RunAction : public G4UserRunAction {
 public:
  explicit RunAction(const AppConfig& config);
  ~RunAction() override;

  void BeginOfRunAction(const G4Run* run) override;
  void EndOfRunAction(const G4Run* run) override;

  void AccumulateEvent(double edepSubstrate, double edepCoating, int nGamma, int nNeutron);

 private:
  AppConfig config_;
  std::mutex mutex_;

  double totalEdepSubstrate_{0.0};
  double totalEdepCoating_{0.0};
  long long totalGamma_{0};
  long long totalNeutron_{0};

#ifdef KSA_USE_ROOT
  class TFile* rootFile_{nullptr};
  class TTree* runTree_{nullptr};
#endif
};
