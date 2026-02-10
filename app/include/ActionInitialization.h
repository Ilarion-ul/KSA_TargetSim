#pragma once

#include "Config.h"

#include <G4VUserActionInitialization.hh>

class ActionInitialization : public G4VUserActionInitialization {
 public:
  explicit ActionInitialization(AppConfig config);
  ~ActionInitialization() override = default;

  void BuildForMaster() const override;
  void Build() const override;

 private:
  AppConfig config_;
};
