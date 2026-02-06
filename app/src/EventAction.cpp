#include "EventAction.h"

#include "RunAction.h"

EventAction::EventAction(RunAction* runAction) : runAction_(runAction) {
  if (runAction_) {
    plateEdep_.assign(runAction_->PlateCount(), 0.0);
    plateNeutronTrackLen_.assign(runAction_->PlateCount(), 0.0);
  }
}

void EventAction::BeginOfEventAction(const G4Event*) {
  edepSubstrate_ = 0.0;
  edepCoating_ = 0.0;
  nGamma_ = 0;
  nNeutron_ = 0;
  nNeutronExit_ = 0;
  gammaTrackIds_.clear();
  neutronTrackIds_.clear();
  neutronExitTrackIds_.clear();
  std::fill(plateEdep_.begin(), plateEdep_.end(), 0.0);
  std::fill(plateNeutronTrackLen_.begin(), plateNeutronTrackLen_.end(), 0.0);
}

void EventAction::EndOfEventAction(const G4Event*) {
  if (runAction_) {
    runAction_->AccumulateEvent(edepSubstrate_, edepCoating_, nGamma_, nNeutron_, nNeutronExit_, plateEdep_,
                                plateNeutronTrackLen_);
  }
}

void EventAction::AddPlateEdep(int plateIndex, double edep) {
  if (plateIndex < 0 || static_cast<size_t>(plateIndex) >= plateEdep_.size()) {
    return;
  }
  plateEdep_[static_cast<size_t>(plateIndex)] += edep;
}

void EventAction::AddPlateNeutronTrackLen(int plateIndex, double stepLen) {
  if (plateIndex < 0 || static_cast<size_t>(plateIndex) >= plateNeutronTrackLen_.size()) {
    return;
  }
  plateNeutronTrackLen_[static_cast<size_t>(plateIndex)] += stepLen;
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

void EventAction::CountNeutronExit(int trackId) {
  if (neutronExitTrackIds_.insert(trackId).second) {
    ++nNeutronExit_;
  }
}
