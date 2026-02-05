#pragma once

#include <string>

struct CliOptions {
  bool showHelp{false};
  bool enableVis{false};
  std::string macroPath;
  std::string configPath{"app/config/default_WTa.json"};
  int nThreadsOverride{-1};
  int nEventsOverride{-1};
  std::string physicsOverride;
};

CliOptions ParseCli(int argc, char** argv);
void PrintHelp();
