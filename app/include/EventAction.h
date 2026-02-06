#pragma once

#include <G4UserEventAction.hh>

#include <vector>
#include <unordered_set>

class G4Event;
class RunAction;

class EventAction : public G4UserEventAction {
 public:
  explicit EventAction(RunAction* runAction);
  ~EventAction() override = default;

  void BeginOfEventAction(const G4Event* event) override;
  void EndOfEventAction(const G4Event* event) override;

  void AddEdepSubstrate(double edep) { edepSubstrate_ += edep; }
  void AddEdepCoating(double edep) { edepCoating_ += edep; }
  void AddPlateEdep(int plateIndex, double edep);
  void AddPlateNeutronTrackLen(int plateIndex, double stepLen);
  void AddPlateNeutronHeatmap(int plateIndex, double xMm, double yMm, double stepLen);
  void AddEdep3d(double xMm, double yMm, double zMm, double edep);
  void CountGamma(int trackId);
  void CountNeutron(int trackId);
  void CountNeutronExit(int trackId);

 private:
  int plateHeatmapBinsX_{0};
  int plateHeatmapBinsY_{0};
  int edepBinsX_{0};
  int edepBinsY_{0};
  int edepBinsZ_{0};
  double plateXMinMm_{0.0};
  double plateXMaxMm_{0.0};
  double plateYMinMm_{0.0};
  double plateYMaxMm_{0.0};
  double edepXMinMm_{0.0};
  double edepXMaxMm_{0.0};
  double edepYMinMm_{0.0};
  double edepYMaxMm_{0.0};
  double edepZMinMm_{0.0};
  double edepZMaxMm_{0.0};

  RunAction* runAction_{nullptr};

  double edepSubstrate_{0.0};
  double edepCoating_{0.0};
  int nGamma_{0};
  int nNeutron_{0};
  int nNeutronExit_{0};
  std::vector<double> plateEdep_;
  std::vector<double> plateNeutronTrackLen_;
  std::vector<double> plateNeutronHeatmap_;
  std::vector<double> edep3d_;

  std::unordered_set<int> gammaTrackIds_;
  std::unordered_set<int> neutronTrackIds_;
  std::unordered_set<int> neutronExitTrackIds_;
};
