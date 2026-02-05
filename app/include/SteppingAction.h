#pragma once

#include <G4UserSteppingAction.hh>

class G4Step;
class EventAction;

class SteppingAction : public G4UserSteppingAction {
 public:
  explicit SteppingAction(EventAction* eventAction);
  ~SteppingAction() override = default;

  void UserSteppingAction(const G4Step* step) override;

 private:
  EventAction* eventAction_{nullptr};

  // TODO: replace this minimal metric with mesh scoring / fluence estimators.
};
