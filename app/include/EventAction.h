#pragma once

#include <G4UserEventAction.hh>

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
  void CountGamma(int trackId);
  void CountNeutron(int trackId);

 private:
  RunAction* runAction_{nullptr};

  double edepSubstrate_{0.0};
  double edepCoating_{0.0};
  int nGamma_{0};
  int nNeutron_{0};

  std::unordered_set<int> gammaTrackIds_;
  std::unordered_set<int> neutronTrackIds_;
};
