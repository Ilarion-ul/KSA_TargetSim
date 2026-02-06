#include "SteppingAction.h"

#include "EventAction.h"

#include <G4ParticleDefinition.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>

namespace {
bool StartsWith(const G4String& name, const char* prefix) {
  const G4String p(prefix);
  return name.size() >= p.size() && name.substr(0, p.size()) == p;
}

bool IsPlateVolume(const G4String& name) {
  return StartsWith(name, "TargetSubstrate") || StartsWith(name, "TargetCoating") || StartsWith(name, "TargetBufferTi") ||
         StartsWith(name, "PlateU_") || StartsWith(name, "PlateCladAl_");
}

bool IsSubstrateVolume(const G4String& name) {
  return StartsWith(name, "TargetSubstrate") || StartsWith(name, "PlateU_");
}

bool IsCoatingVolume(const G4String& name) {
  return StartsWith(name, "TargetCoating") || StartsWith(name, "TargetBufferTi") || StartsWith(name, "PlateCladAl_");
}

int PlateIndexFromCopyNo(int copyNo) {
  if (copyNo > 0) {
    return copyNo - 1;
  }
  return 0;
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

  if (!IsPlateVolume(name)) {
    return;
  }

  const auto* prePoint = step->GetPreStepPoint();
  const auto globalPos = prePoint->GetPosition();
  const auto localPos = prePoint->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(globalPos);

  if (IsSubstrateVolume(name)) {
    eventAction_->AddEdepSubstrate(edep);
  } else if (IsCoatingVolume(name)) {
    eventAction_->AddEdepCoating(edep);
  }
  const int plateIndex = PlateIndexFromCopyNo(step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber());
  eventAction_->AddPlateEdep(plateIndex, edep);
  if (edep > 0.0) {
    eventAction_->AddEdep3d(globalPos.x() / mm, globalPos.y() / mm, globalPos.z() / mm, edep);
  }

  // Minimal particle counting metric:
  // count unique gamma/neutron track IDs that appear in target volume.
  const auto* track = step->GetTrack();
  if (!track || !track->GetDefinition()) return;

  const auto* p = track->GetDefinition();
  const auto particleName = p->GetParticleName();
  if (particleName == "gamma") {
    eventAction_->CountGamma(track->GetTrackID());
  } else if (particleName == "neutron") {
    eventAction_->CountNeutron(track->GetTrackID());
    const double stepLen = step->GetStepLength();
    eventAction_->AddPlateNeutronTrackLen(plateIndex, stepLen);
    eventAction_->AddPlateNeutronHeatmap(plateIndex, localPos.x() / mm, localPos.y() / mm, stepLen);
    if (step->GetPostStepPoint()) {
      const auto* postVolume = step->GetPostStepPoint()->GetTouchableHandle()
                                   ? step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()
                                   : nullptr;
      const auto postName = postVolume ? postVolume->GetName() : G4String();
      if (!postVolume || !IsPlateVolume(postName)) {
        eventAction_->CountNeutronExit(track->GetTrackID());
      }
    }
  }
}
