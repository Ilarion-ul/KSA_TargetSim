#include "PhysicsListFactory.h"

#include <G4PhysListFactory.hh>
#include <G4SystemOfUnits.hh>
#include <G4VModularPhysicsList.hh>
#include <G4VUserPhysicsList.hh>

#include <iostream>

G4VUserPhysicsList* CreatePhysicsList(const PhysicsConfig& config) {
  G4PhysListFactory factory;
  const std::string requested = config.physicsListName.empty() ? "QGSP_BIC_HPT" : config.physicsListName;

  G4VModularPhysicsList* list = nullptr;
  if (factory.IsReferencePhysList(requested)) {
    list = factory.GetReferencePhysList(requested);
  } else {
    std::cerr << "[physics] Unknown physics list '" << requested
              << "', fallback to QGSP_BIC_HPT\n";
    list = factory.GetReferencePhysList("QGSP_BIC_HPT");
  }

  if (config.cut_mm.has_value()) {
    list->SetDefaultCutValue(config.cut_mm.value() * mm);
  }

  return list;
}
