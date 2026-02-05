#include "RunAction.h"

#include <G4Run.hh>
#include <G4SystemOfUnits.hh>
#include <G4Threading.hh>

#if __has_include(<filesystem>)
#include <filesystem>
#define KSA_HAVE_STD_FILESYSTEM 1
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
#define KSA_HAVE_STD_FILESYSTEM 0
#else
#error "No filesystem support available"
#endif

#include <fstream>

#ifdef KSA_USE_ROOT
#include <TFile.h>
#include <TTree.h>
#endif

RunAction::RunAction(const AppConfig& config) : config_(config) {}
RunAction::~RunAction() = default;

namespace {
#if KSA_HAVE_STD_FILESYSTEM
namespace fs = std::filesystem;
#else
namespace fs = std::experimental::filesystem;
#endif

// Shared accumulators across worker/master RunAction instances.
// In MT Geant4 creates separate action objects per thread, so member counters
// are not automatically merged.
std::mutex gRunSummaryMutex;
double gTotalEdepSubstrate = 0.0;
double gTotalEdepCoating = 0.0;
long long gTotalGamma = 0;
long long gTotalNeutron = 0;
}

void RunAction::BeginOfRunAction(const G4Run*) {
  // Reset global counters once at run start by master only.
  if (G4Threading::IsMasterThread()) {
    std::lock_guard<std::mutex> lock(gRunSummaryMutex);
    gTotalEdepSubstrate = 0.0;
    gTotalEdepCoating = 0.0;
    gTotalGamma = 0;
    gTotalNeutron = 0;
  }

  const auto base = fs::path(config_.run.outputDir.empty() ? "results" : config_.run.outputDir);
  fs::create_directories(base / "logs");
  fs::create_directories(base / "root");
  fs::create_directories(base / "vis");

#ifdef KSA_USE_ROOT
  // ROOT global state is not thread-safe for opening/writing one file from
  // all workers; keep file creation/writing on master only.
  if (G4Threading::IsMasterThread()) {
    const auto rootPath = base / "root" / config_.run.outputRootFile;
    rootFile_ = TFile::Open(rootPath.string().c_str(), "RECREATE");
    runTree_ = new TTree("run_summary", "KSA run summary");
  }
#endif
}

void RunAction::EndOfRunAction(const G4Run* run) {
  // Master writes final summary once. Workers only contribute counters.
  if (!G4Threading::IsMasterThread()) {
    return;
  }

  const auto base = fs::path(config_.run.outputDir.empty() ? "results" : config_.run.outputDir);

  double edepSub = 0.0;
  double edepCoat = 0.0;
  long long nGammaTotal = 0;
  long long nNeutronTotal = 0;
  {
    std::lock_guard<std::mutex> lock(gRunSummaryMutex);
    edepSub = gTotalEdepSubstrate;
    edepCoat = gTotalEdepCoating;
    nGammaTotal = gTotalGamma;
    nNeutronTotal = gTotalNeutron;
  }

#ifdef KSA_USE_ROOT
  if (rootFile_ && runTree_) {
    double edep_substrate_MeV = edepSub / MeV;
    double edep_coating_MeV = edepCoat / MeV;
    long long nGamma = nGammaTotal;
    long long nNeutron = nNeutronTotal;
    int nEvents = run ? run->GetNumberOfEvent() : config_.run.nEvents;
    std::string physicsListName = config_.physics.physicsListName;

    runTree_->Branch("edep_substrate", &edep_substrate_MeV);
    runTree_->Branch("edep_coating", &edep_coating_MeV);
    runTree_->Branch("nGamma", &nGamma);
    runTree_->Branch("nNeutron", &nNeutron);
    runTree_->Branch("nEvents", &nEvents);
    runTree_->Branch("physicsListName", &physicsListName);
    runTree_->Fill();

    rootFile_->cd();
    runTree_->Write();
    rootFile_->Close();
    delete rootFile_;
    rootFile_ = nullptr;
    runTree_ = nullptr;
  }
#else
  const auto out = base / "logs" / "run_summary.json";
  std::ofstream os(out);
  os << "{\n";
  os << "  \"edep_substrate\": " << (edepSub / MeV) << ",\n";
  os << "  \"edep_coating\": " << (edepCoat / MeV) << ",\n";
  os << "  \"nGamma\": " << nGammaTotal << ",\n";
  os << "  \"nNeutron\": " << nNeutronTotal << ",\n";
  os << "  \"nEvents\": " << (run ? run->GetNumberOfEvent() : config_.run.nEvents) << ",\n";
  os << "  \"physicsListName\": \"" << config_.physics.physicsListName << "\"\n";
  os << "}\n";
#endif
}

void RunAction::AccumulateEvent(double edepSubstrate, double edepCoating, int nGamma, int nNeutron) {
  std::lock_guard<std::mutex> lock(gRunSummaryMutex);
  gTotalEdepSubstrate += edepSubstrate;
  gTotalEdepCoating += edepCoating;
  gTotalGamma += nGamma;
  gTotalNeutron += nNeutron;
  // TODO: migrate to G4Accumulable/G4AnalysisManager for richer MT-safe reporting.
}
