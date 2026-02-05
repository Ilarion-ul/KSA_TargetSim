#include "Cli.h"

#include <iostream>
#include <stdexcept>

CliOptions ParseCli(int argc, char** argv) {
  CliOptions opt;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      opt.showHelp = true;
    } else if (arg == "--vis") {
      opt.enableVis = true; // --vis overrides batch behavior in main.
    } else if (arg == "-m") {
      if (i + 1 >= argc) throw std::runtime_error("Option -m requires a macro path");
      opt.macroPath = argv[++i];
    } else if (arg == "-c") {
      if (i + 1 >= argc) throw std::runtime_error("Option -c requires a config path");
      opt.configPath = argv[++i];
    } else if (arg == "-t") {
      if (i + 1 >= argc) throw std::runtime_error("Option -t requires a thread count");
      opt.nThreadsOverride = std::stoi(argv[++i]);
    } else if (arg == "-n") {
      if (i + 1 >= argc) throw std::runtime_error("Option -n requires number of events");
      opt.nEventsOverride = std::stoi(argv[++i]);
    } else if (arg == "--physics") {
      if (i + 1 >= argc) throw std::runtime_error("Option --physics requires a list name");
      opt.physicsOverride = argv[++i];
    } else {
      throw std::runtime_error("Unknown option: " + arg + " (use --help)");
    }
  }

  return opt;
}

void PrintHelp() {
  std::cout
      << "KSA_TargetSim (ksasim)\n"
      << "Usage:\n"
      << "  ksasim [options]\n\n"
      << "Options:\n"
      << "  --help, -h              Show this help\n"
      << "  --vis                   Run in visualization mode (overrides batch macro)\n"
      << "  -m <macro.mac>          Batch macro path\n"
      << "  -c <config.json>        Configuration JSON path (default: app/config/default_WTa.json)\n"
      << "  -t <nThreads>           Override thread count\n"
      << "  -n <nEvents>            Override number of events\n"
      << "  --physics <name>        Override Geant4 physics list name\n\n"
      << "Examples:\n"
      << "  ksasim --help\n"
      << "  ksasim -c app/config/default_WTa.json -m app/macros/run_batch.mac -n 100\n"
      << "  ksasim --vis -c app/config/quick_vis.json\n";
}
