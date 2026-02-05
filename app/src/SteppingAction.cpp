#include "SteppingAction.h"

#include "EventAction.h"

#include <G4ParticleDefinition.hh>
#include <G4Step.hh>
#include <G4Track.hh>

namespace {
bool StartsWith(const G4String& name, const char* prefix) {
  const G4String p(prefix);
  return name.size() >= p.size() && name.substr(0, p.size()) == p;
}
}

SteppingAction::SteppingAction(EventAction* eventAction) : eventAction_(eventAction) {}

void SteppingAction::UserSteppingAction(const G4Step* step) {
  if (!eventAction_ || !step || !step->GetPreStepPoint() || !step->GetPreStepPoint()->GetTouchableHandle()) {
    return;
  }

  const auto* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  if (!volume) return;

  const auto name = volume->GetName();
  const double edep = step->GetTotalEnergyDeposit();

  const bool inSub = StartsWith(name, "TargetSubstrate");
  const bool inCoat = StartsWith(name, "TargetCoating");
  if (!(inSub || inCoat)) {
    return;
  }

  if (inSub) {
    eventAction_->AddEdepSubstrate(edep);
  } else if (inCoat) {
    eventAction_->AddEdepCoating(edep);
  }

  // Minimal particle counting metric:
  // count unique gamma/neutron track IDs that appear in target volume.
  const auto* track = step->GetTrack();
  if (!track || !track->GetDefinition()) return;

  const auto* p = track->GetDefinition();
  if (p->GetParticleName() == "gamma") {
    eventAction_->CountGamma(track->GetTrackID());
  } else if (p->GetParticleName() == "neutron") {
    eventAction_->CountNeutron(track->GetTrackID());
  }
}
