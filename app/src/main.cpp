#include "ActionInitialization.h"
#include "Cli.h"
#include "Config.h"
#include "DetectorConstruction.h"
#include "PhysicsListFactory.h"

#include <G4RunManagerFactory.hh>
#include <G4UImanager.hh>
#include <G4UIExecutive.hh>
#include <G4VisExecutive.hh>
#include <G4VisManager.hh>

#include <iostream>
#include <memory>

int main(int argc, char** argv) {
  try {
    auto cli = ParseCli(argc, argv);
    if (cli.showHelp) {
      PrintHelp();
      return 0;
    }

    AppConfig cfg = LoadConfig(cli.configPath.empty() ? "app/config/default_WTa.json" : cli.configPath);

    if (cli.nEventsOverride > 0) cfg.run.nEvents = cli.nEventsOverride;
    if (cli.nThreadsOverride >= 0) cfg.run.nThreads = cli.nThreadsOverride;
    if (!cli.physicsOverride.empty()) cfg.physics.physicsListName = cli.physicsOverride;
    if (cli.enableVis) cfg.run.enableVis = true;

    std::cout << "[ksasim] physics=" << cfg.physics.physicsListName
              << ", target=" << cfg.target.type
              << ", Ebeam(MeV)=" << cfg.beam.energy_MeV
              << ", nEvents=" << cfg.run.nEvents << "\n";

    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    if (cfg.run.nThreads > 0) {
      runManager->SetNumberOfThreads(cfg.run.nThreads);
    }

    runManager->SetUserInitialization(new DetectorConstruction(cfg));
    runManager->SetUserInitialization(CreatePhysicsList(cfg.physics));
    runManager->SetUserInitialization(new ActionInitialization(cfg));

    auto* ui = G4UImanager::GetUIpointer();

    if (cfg.run.enableVis) {
      G4UIExecutive uiExec(argc, argv);
      G4VisManager* visManager = new G4VisExecutive();
      visManager->Initialize();

      const std::string visMacro = cli.macroPath.empty() ? "app/macros/vis_ogl.mac" : cli.macroPath;
      ui->ApplyCommand("/control/execute " + visMacro);
      uiExec.SessionStart();

      delete visManager;
    } else {
      runManager->Initialize();
      if (!cli.macroPath.empty()) {
        ui->ApplyCommand("/control/execute " + cli.macroPath);
      } else {
        ui->ApplyCommand("/control/execute app/macros/run_batch.mac");
        if (cfg.run.nEvents > 0) {
          ui->ApplyCommand("/run/beamOn " + std::to_string(cfg.run.nEvents));
        }
      }
    }

    delete runManager;
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "[ksasim] Fatal error: " << ex.what() << "\n";
    return 2;
  }
}
