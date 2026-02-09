#include "EventAction.h"

#include "RunAction.h"

#include <G4Event.hh>

namespace {
int ClampBin(double value, double minValue, double maxValue, int bins) {
  if (bins <= 0 || value < minValue || value >= maxValue) {
    return -1;
  }
  const double span = maxValue - minValue;
  if (span <= 0.0) {
    return -1;
  }
  const double scaled = (value - minValue) / span;
  int index = static_cast<int>(scaled * bins);
  if (index < 0) index = 0;
  if (index >= bins) index = bins - 1;
  return index;
}
} // namespace

EventAction::EventAction(RunAction* runAction) : runAction_(runAction) {
  if (runAction_) {
    plateEdep_.assign(runAction_->PlateCount(), 0.0);
    plateNeutronTrackLen_.assign(runAction_->PlateCount(), 0.0);

    plateHeatmapBinsX_ = runAction_->PlateHeatmapBinsX();
    plateHeatmapBinsY_ = runAction_->PlateHeatmapBinsY();
    edepBinsX_ = runAction_->EdepBinsX();
    edepBinsY_ = runAction_->EdepBinsY();
    edepBinsZ_ = runAction_->EdepBinsZ();
    const auto plateBounds = runAction_->PlateHeatmapBounds();
    plateXMinMm_ = plateBounds.xMinMm;
    plateXMaxMm_ = plateBounds.xMaxMm;
    plateYMinMm_ = plateBounds.yMinMm;
    plateYMaxMm_ = plateBounds.yMaxMm;
    const auto edepBounds = runAction_->EdepBounds();
    edepXMinMm_ = edepBounds.xMinMm;
    edepXMaxMm_ = edepBounds.xMaxMm;
    edepYMinMm_ = edepBounds.yMinMm;
    edepYMaxMm_ = edepBounds.yMaxMm;
    edepZMinMm_ = edepBounds.zMinMm;
    edepZMaxMm_ = edepBounds.zMaxMm;
    targetXMinMm_ = edepBounds.xMinMm;
    targetXMaxMm_ = edepBounds.xMaxMm;
    targetYMinMm_ = edepBounds.yMinMm;
    targetYMaxMm_ = edepBounds.yMaxMm;
    targetZMinMm_ = edepBounds.zMinMm;
    targetZMaxMm_ = edepBounds.zMaxMm;

    const size_t plateHeatmapSize = runAction_->PlateCount() * static_cast<size_t>(plateHeatmapBinsX_) *
                                    static_cast<size_t>(plateHeatmapBinsY_);
    plateNeutronHeatmap_.assign(plateHeatmapSize, 0.0);
    const size_t edep3dSize = static_cast<size_t>(edepBinsX_) * static_cast<size_t>(edepBinsY_) *
                              static_cast<size_t>(edepBinsZ_);
    edep3d_.assign(edep3dSize, 0.0);
  }
}

void EventAction::BeginOfEventAction(const G4Event* event) {
  eventId_ = event ? event->GetEventID() : 0;
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
  std::fill(plateNeutronHeatmap_.begin(), plateNeutronHeatmap_.end(), 0.0);
  std::fill(edep3d_.begin(), edep3d_.end(), 0.0);
  neutronSurfaceHits_.clear();
}

void EventAction::EndOfEventAction(const G4Event*) {
  if (runAction_) {
    runAction_->AccumulateEvent(edepSubstrate_, edepCoating_, nGamma_, nNeutron_, nNeutronExit_, eventId_, plateEdep_,
                                plateNeutronTrackLen_, plateNeutronHeatmap_, edep3d_, neutronSurfaceHits_);
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

void EventAction::AddPlateNeutronHeatmap(int plateIndex, double xMm, double yMm, double stepLen) {
  if (!runAction_ || plateIndex < 0 || static_cast<size_t>(plateIndex) >= plateEdep_.size() || plateHeatmapBinsX_ <= 0 ||
      plateHeatmapBinsY_ <= 0) {
    return;
  }
  const int ix = ClampBin(xMm, plateXMinMm_, plateXMaxMm_, plateHeatmapBinsX_);
  const int iy = ClampBin(yMm, plateYMinMm_, plateYMaxMm_, plateHeatmapBinsY_);
  if (ix < 0 || iy < 0) return;
  const size_t plateOffset = static_cast<size_t>(plateIndex) * static_cast<size_t>(plateHeatmapBinsX_) *
                             static_cast<size_t>(plateHeatmapBinsY_);
  const size_t idx = plateOffset + static_cast<size_t>(iy) * static_cast<size_t>(plateHeatmapBinsX_) +
                     static_cast<size_t>(ix);
  if (idx >= plateNeutronHeatmap_.size()) return;
  plateNeutronHeatmap_[idx] += stepLen;
}

void EventAction::AddEdep3d(double xMm, double yMm, double zMm, double edep) {
  if (edepBinsX_ <= 0 || edepBinsY_ <= 0 || edepBinsZ_ <= 0) {
    return;
  }
  const int ix = ClampBin(xMm, edepXMinMm_, edepXMaxMm_, edepBinsX_);
  const int iy = ClampBin(yMm, edepYMinMm_, edepYMaxMm_, edepBinsY_);
  const int iz = ClampBin(zMm, edepZMinMm_, edepZMaxMm_, edepBinsZ_);
  if (ix < 0 || iy < 0 || iz < 0) return;
  const size_t idx = (static_cast<size_t>(iz) * static_cast<size_t>(edepBinsY_) + static_cast<size_t>(iy)) *
                         static_cast<size_t>(edepBinsX_) +
                     static_cast<size_t>(ix);
  if (idx >= edep3d_.size()) return;
  edep3d_[idx] += edep;
}

void EventAction::AddNeutronSurfaceHit(double EnMeV,
                                       double xMm,
                                       double yMm,
                                       double zMm,
                                       double cosTheta,
                                       double weight,
                                       double timeNs,
                                       int surfaceId) {
  RunAction::NeutronSurfaceHit hit;
  hit.event_id = eventId_;
  hit.En_MeV = EnMeV;
  hit.x_mm = xMm;
  hit.y_mm = yMm;
  hit.z_mm = zMm;
  hit.cosTheta = cosTheta;
  hit.weight = weight;
  hit.time_ns = timeNs;
  hit.surface_id = surfaceId;
  neutronSurfaceHits_.push_back(hit);
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
