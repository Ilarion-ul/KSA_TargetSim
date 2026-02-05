#include "RunAction.h"

#include <G4Run.hh>
#include <G4SystemOfUnits.hh>

#include <filesystem>
#include <fstream>

#ifdef KSA_USE_ROOT
#include <TFile.h>
#include <TTree.h>
#endif

RunAction::RunAction(const AppConfig& config) : config_(config) {}
RunAction::~RunAction() = default;

void RunAction::BeginOfRunAction(const G4Run*) {
  totalEdepSubstrate_ = 0.0;
  totalEdepCoating_ = 0.0;
  totalGamma_ = 0;
  totalNeutron_ = 0;

  const auto base = std::filesystem::path(config_.run.outputDir.empty() ? "results" : config_.run.outputDir);
  std::filesystem::create_directories(base / "logs");
  std::filesystem::create_directories(base / "root");
  std::filesystem::create_directories(base / "vis");

#ifdef KSA_USE_ROOT
  const auto rootPath = base / "root" / config_.run.outputRootFile;
  rootFile_ = TFile::Open(rootPath.string().c_str(), "RECREATE");
  runTree_ = new TTree("run_summary", "KSA run summary");
#endif
}

void RunAction::EndOfRunAction(const G4Run* run) {
  const auto base = std::filesystem::path(config_.run.outputDir.empty() ? "results" : config_.run.outputDir);

#ifdef KSA_USE_ROOT
  if (rootFile_ && runTree_) {
    double edep_substrate_MeV = totalEdepSubstrate_ / MeV;
    double edep_coating_MeV = totalEdepCoating_ / MeV;
    long long nGamma = totalGamma_;
    long long nNeutron = totalNeutron_;
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
  os << "  \"edep_substrate\": " << (totalEdepSubstrate_ / MeV) << ",\n";
  os << "  \"edep_coating\": " << (totalEdepCoating_ / MeV) << ",\n";
  os << "  \"nGamma\": " << totalGamma_ << ",\n";
  os << "  \"nNeutron\": " << totalNeutron_ << ",\n";
  os << "  \"nEvents\": " << (run ? run->GetNumberOfEvent() : config_.run.nEvents) << ",\n";
  os << "  \"physicsListName\": \"" << config_.physics.physicsListName << "\"\n";
  os << "}\n";
#endif
}

void RunAction::AccumulateEvent(double edepSubstrate, double edepCoating, int nGamma, int nNeutron) {
  std::lock_guard<std::mutex> lock(mutex_);
  totalEdepSubstrate_ += edepSubstrate;
  totalEdepCoating_ += edepCoating;
  totalGamma_ += nGamma;
  totalNeutron_ += nNeutron;
  // TODO: migrate to G4Accumulable/G4AnalysisManager for richer MT-safe reporting.
}
