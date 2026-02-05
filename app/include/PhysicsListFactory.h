#pragma once

#include "Config.h"

class G4VUserPhysicsList;

G4VUserPhysicsList* CreatePhysicsList(const PhysicsConfig& config);

// TODO: add fine-grained EM controls (bremsstrahlung knobs) and explicit
// photonuclear toggles when integrating custom physics constructors.
