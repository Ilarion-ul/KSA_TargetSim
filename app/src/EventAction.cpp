#include "EventAction.h"

#include "RunAction.h"

EventAction::EventAction(RunAction* runAction) : runAction_(runAction) {}

void EventAction::BeginOfEventAction(const G4Event*) {
  edepSubstrate_ = 0.0;
  edepCoating_ = 0.0;
  nGamma_ = 0;
  nNeutron_ = 0;
  gammaTrackIds_.clear();
  neutronTrackIds_.clear();
}

void EventAction::EndOfEventAction(const G4Event*) {
  if (runAction_) {
    runAction_->AccumulateEvent(edepSubstrate_, edepCoating_, nGamma_, nNeutron_);
  }
}

void EventAction::CountGamma(int trackId) {
  if (gammaTrackIds_.insert(trackId).second) {
    ++nGamma_;
  }
}

void EventAction::CountNeutron(int trackId) {
  if (neutronTrackIds_.insert(trackId).second) {
    ++nNeutron_;
  }
}
