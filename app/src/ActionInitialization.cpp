#include "ActionInitialization.h"

#include "EventAction.h"
#include "PrimaryGeneratorAction.h"
#include "RunAction.h"
#include "SteppingAction.h"

ActionInitialization::ActionInitialization(AppConfig config) : config_(std::move(config)) {}

void ActionInitialization::BuildForMaster() const {
  SetUserAction(new RunAction(config_));
}

void ActionInitialization::Build() const {
  auto* runAction = new RunAction(config_);
  SetUserAction(runAction);

  auto* eventAction = new EventAction(runAction);
  SetUserAction(eventAction);

  SetUserAction(new PrimaryGeneratorAction(config_));
  SetUserAction(new SteppingAction(eventAction));
}
