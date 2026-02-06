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
    double beam_E0_MeV = config_.beam.energy_MeV;
    double beam_Esigma_rel = config_.beam.energy_sigma_rel_1sigma;
    std::string beam_spread_model = config_.beam.energy_spread_model;
    std::string beam_mode = config_.beam.mode;
    double beam_pulse_width_us = config_.beam.pulse_width_us;
    double beam_rep_rate_Hz = config_.beam.rep_rate_Hz;
    double beam_I_pulse_A = config_.beam.I_pulse_A;
    double beam_I_avg_A_input = config_.beam.I_avg_A;
    double beam_duty = (config_.beam.pulse_width_us * 1e-6) * config_.beam.rep_rate_Hz;
    double beam_I_avg_A_calc = beam_I_pulse_A * beam_duty;
    double beam_power_kW = config_.beam.beam_power_kW;

    runTree_->Branch("edep_substrate", &edep_substrate_MeV);
    runTree_->Branch("edep_coating", &edep_coating_MeV);
    runTree_->Branch("nGamma", &nGamma);
    runTree_->Branch("nNeutron", &nNeutron);
    runTree_->Branch("nEvents", &nEvents);
    runTree_->Branch("physicsListName", &physicsListName);
    runTree_->Branch("beam_E0_MeV", &beam_E0_MeV);
    runTree_->Branch("beam_Esigma_rel", &beam_Esigma_rel);
    runTree_->Branch("beam_spread_model", &beam_spread_model);
    runTree_->Branch("beam_mode", &beam_mode);
    runTree_->Branch("beam_pulse_width_us", &beam_pulse_width_us);
    runTree_->Branch("beam_rep_rate_Hz", &beam_rep_rate_Hz);
    runTree_->Branch("beam_I_pulse_A", &beam_I_pulse_A);
    runTree_->Branch("beam_I_avg_A_input", &beam_I_avg_A_input);
    runTree_->Branch("beam_duty", &beam_duty);
    runTree_->Branch("beam_I_avg_A_calc", &beam_I_avg_A_calc);
    runTree_->Branch("beam_power_kW", &beam_power_kW);
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
  os << "  \"physicsListName\": \"" << config_.physics.physicsListName << "\",\n";
  os << "  \"beam\": {\n";
  os << "    \"E0_MeV\": " << config_.beam.energy_MeV << ",\n";
  os << "    \"energy_spread_model\": \"" << config_.beam.energy_spread_model << "\",\n";
  os << "    \"energy_sigma_rel_1sigma\": " << config_.beam.energy_sigma_rel_1sigma << ",\n";
  os << "    \"mode\": \"" << config_.beam.mode << "\",\n";
  os << "    \"pulse_width_us\": " << config_.beam.pulse_width_us << ",\n";
  os << "    \"rep_rate_Hz\": " << config_.beam.rep_rate_Hz << ",\n";
  os << "    \"I_pulse_A\": " << config_.beam.I_pulse_A << ",\n";
  os << "    \"I_avg_A_input\": " << config_.beam.I_avg_A << ",\n";
  os << "    \"I_avg_A_calc\": " << (config_.beam.I_pulse_A * (config_.beam.pulse_width_us * 1e-6) * config_.beam.rep_rate_Hz) << ",\n";
  os << "    \"beam_power_kW\": " << config_.beam.beam_power_kW << "\n";
  os << "  }\n";
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
