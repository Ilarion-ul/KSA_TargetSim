#include <TCanvas.h>
#include <TDatime.h>
#include <TFile.h>
#include <TKey.h>
#include <TObjString.h>
#include <TObject.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TH1D.h>
#include <TStyle.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace {

void PrepareTreeForSafeRead(TTree* tree) {
  if (!tree) return;
  tree->SetCacheSize(0);
  tree->SetBranchStatus("*", 0);
}

void EnableBranches(TTree* tree, const std::vector<std::string>& names) {
  if (!tree) return;
  for (const auto& n : names) {
    if (tree->GetBranch(n.c_str())) {
      tree->SetBranchStatus(n.c_str(), 1);
    }
  }
}

bool SafeGetEntry(TTree* tree, Long64_t entry, const char* label) {
  if (!tree) return false;
  if (entry < 2) {
    ::Info("export_run_artifacts", "Reading %s entry %lld", label, entry);
  }
  const Long64_t bytes = tree->GetEntry(entry);
  if (bytes < 0) {
    ::Error("export_run_artifacts", "Failed to read %s entry %lld", label, entry);
    return false;
  }
  return true;
}

bool ReadNumericRunSummaryValue(TTree* tree, const char* branchName, long long& outValue) {
  if (!tree || !branchName || !tree->GetBranch(branchName)) {
    return false;
  }
  const Long64_t n = tree->Draw(branchName, "", "goff", 1, 0);
  if (n <= 0 || !tree->GetV1()) {
    return false;
  }
  outValue = static_cast<long long>(tree->GetV1()[0]);
  return true;
}

std::string Sanitize(const std::string& value) {
  std::string out = value;
  std::replace(out.begin(), out.end(), ' ', '_');
  std::replace(out.begin(), out.end(), '/', '_');
  return out;
}

std::string TimestampFromFile(TFile* file, const std::string& path) {
  TDatime dt = file->GetCreationDate();
  if (dt.GetYear() > 0) {
    std::ostringstream os;
    os << std::setw(4) << std::setfill('0') << dt.GetYear()
       << std::setw(2) << std::setfill('0') << dt.GetMonth()
       << std::setw(2) << std::setfill('0') << dt.GetDay()
       << "_"
       << std::setw(2) << std::setfill('0') << dt.GetHour()
       << std::setw(2) << std::setfill('0') << dt.GetMinute()
       << std::setw(2) << std::setfill('0') << dt.GetSecond();
    return os.str();
  }
  Long_t id, size, flags, modtime;
  if (gSystem->GetPathInfo(path.c_str(), &id, &size, &flags, &modtime) == 0) {
    TDatime mdt(modtime);
    std::ostringstream os;
    os << std::setw(4) << std::setfill('0') << mdt.GetYear()
       << std::setw(2) << std::setfill('0') << mdt.GetMonth()
       << std::setw(2) << std::setfill('0') << mdt.GetDay()
       << "_"
       << std::setw(2) << std::setfill('0') << mdt.GetHour()
       << std::setw(2) << std::setfill('0') << mdt.GetMinute()
       << std::setw(2) << std::setfill('0') << mdt.GetSecond();
    return os.str();
  }
  return "unknown_time";
}

std::string ReadTargetType(TFile* file) {
  auto* tree = dynamic_cast<TTree*>(file->Get("RunMeta"));
  if (!tree) return "unknown";
  PrepareTreeForSafeRead(tree);
  EnableBranches(tree, {"target_type"});
  std::string* target_type = nullptr;
  if (!tree->GetBranch("target_type")) return "unknown";
  tree->SetBranchAddress("target_type", &target_type);
  tree->Print();
  if (!SafeGetEntry(tree, 0, "RunMeta(target_type)")) return "unknown";
  const std::string value = target_type ? *target_type : "unknown";
  return value.empty() ? "unknown" : value;
}

void WriteRunMetaJson(TFile* file, const std::string& outPath) {
  auto* tree = dynamic_cast<TTree*>(file->Get("RunMeta"));
  if (!tree) return;
  PrepareTreeForSafeRead(tree);
  EnableBranches(tree, {"target_type", "geometry_json", "beam_energy_MeV", "beam_energy_sigma_rel", "beam_sigma_x_mm",
                        "beam_sigma_y_mm", "beam_sigma_theta_x_mrad", "beam_sigma_theta_y_mrad", "pulse_mode",
                        "pulse_width_us", "rep_rate_Hz", "I_pulse_A", "I_avg_A", "duty", "P_avg_kW", "nEvents",
                        "nThreads", "per_primary", "N_e_per_s"});
  std::string* target_type = nullptr;
  std::string* geometry_json = nullptr;
  double beam_energy_MeV = 0.0;
  double beam_energy_sigma_rel = 0.0;
  double beam_sigma_x_mm = 0.0;
  double beam_sigma_y_mm = 0.0;
  double beam_sigma_theta_x_mrad = 0.0;
  double beam_sigma_theta_y_mrad = 0.0;
  std::string* pulse_mode = nullptr;
  double pulse_width_us = 0.0;
  double rep_rate_Hz = 0.0;
  double I_pulse_A = 0.0;
  double I_avg_A = 0.0;
  double duty = 0.0;
  double P_avg_kW = 0.0;
  int nEvents = 0;
  int nThreads = 0;
  double per_primary = 1.0;
  double N_e_per_s = 0.0;
  if (!tree->GetBranch("target_type") || !tree->GetBranch("geometry_json") || !tree->GetBranch("pulse_mode")) {
    return;
  }
  tree->SetBranchAddress("target_type", &target_type);
  tree->SetBranchAddress("geometry_json", &geometry_json);
  tree->SetBranchAddress("beam_energy_MeV", &beam_energy_MeV);
  tree->SetBranchAddress("beam_energy_sigma_rel", &beam_energy_sigma_rel);
  tree->SetBranchAddress("beam_sigma_x_mm", &beam_sigma_x_mm);
  tree->SetBranchAddress("beam_sigma_y_mm", &beam_sigma_y_mm);
  tree->SetBranchAddress("beam_sigma_theta_x_mrad", &beam_sigma_theta_x_mrad);
  tree->SetBranchAddress("beam_sigma_theta_y_mrad", &beam_sigma_theta_y_mrad);
  tree->SetBranchAddress("pulse_mode", &pulse_mode);
  tree->SetBranchAddress("pulse_width_us", &pulse_width_us);
  tree->SetBranchAddress("rep_rate_Hz", &rep_rate_Hz);
  tree->SetBranchAddress("I_pulse_A", &I_pulse_A);
  tree->SetBranchAddress("I_avg_A", &I_avg_A);
  tree->SetBranchAddress("duty", &duty);
  tree->SetBranchAddress("P_avg_kW", &P_avg_kW);
  tree->SetBranchAddress("nEvents", &nEvents);
  tree->SetBranchAddress("nThreads", &nThreads);
  tree->SetBranchAddress("per_primary", &per_primary);
  tree->SetBranchAddress("N_e_per_s", &N_e_per_s);
  tree->Print();
  if (!SafeGetEntry(tree, 0, "RunMeta")) return;

  const std::string tt = target_type ? *target_type : "unknown";
  const std::string gj = geometry_json ? *geometry_json : "{}";
  const std::string pm = pulse_mode ? *pulse_mode : "unknown";


  std::ofstream os(outPath);
  os << "{\n";
  os << "  \"target_type\": \"" << tt << "\",\n";
  os << "  \"geometry\": " << gj << ",\n";
  os << "  \"beam\": {\n";
  os << "    \"energy_MeV\": " << beam_energy_MeV << ",\n";
  os << "    \"energy_sigma_rel\": " << beam_energy_sigma_rel << ",\n";
  os << "    \"sigma_x_mm\": " << beam_sigma_x_mm << ",\n";
  os << "    \"sigma_y_mm\": " << beam_sigma_y_mm << ",\n";
  os << "    \"sigma_theta_x_mrad\": " << beam_sigma_theta_x_mrad << ",\n";
  os << "    \"sigma_theta_y_mrad\": " << beam_sigma_theta_y_mrad << "\n";
  os << "  },\n";
  os << "  \"pulse\": {\n";
  os << "    \"mode\": \"" << pm << "\",\n";
  os << "    \"pulse_width_us\": " << pulse_width_us << ",\n";
  os << "    \"rep_rate_Hz\": " << rep_rate_Hz << ",\n";
  os << "    \"I_pulse_A\": " << I_pulse_A << ",\n";
  os << "    \"I_avg_A\": " << I_avg_A << ",\n";
  os << "    \"duty\": " << duty << ",\n";
  os << "    \"P_avg_kW\": " << P_avg_kW << "\n";
  os << "  },\n";
  os << "  \"run\": {\n";
  os << "    \"nEvents\": " << nEvents << ",\n";
  os << "    \"nThreads\": " << nThreads << "\n";
  os << "  },\n";
  os << "  \"normalization\": {\n";
  os << "    \"per_primary\": " << per_primary << ",\n";
  os << "    \"N_e_per_s\": " << N_e_per_s << "\n";
  os << "  }\n";
  os << "}\n";
}

void WriteNeutronSurfCsv(TFile* file, const std::string& outPath) {
  auto* tree = dynamic_cast<TTree*>(file->Get("NeutronSurf"));
  if (!tree) return;
  PrepareTreeForSafeRead(tree);
  EnableBranches(tree, {"event_id", "En_MeV", "x_mm", "y_mm", "z_mm", "cosTheta", "weight", "time_ns", "surface_id", "surface_name"});
  int event_id = 0;
  double En_MeV = 0.0;
  double x_mm = 0.0;
  double y_mm = 0.0;
  double z_mm = 0.0;
  double cosTheta = 0.0;
  double weight = 1.0;
  double time_ns = 0.0;
  int surface_id = 0;
  std::string* surface_name = nullptr;
  tree->SetBranchAddress("event_id", &event_id);
  tree->SetBranchAddress("En_MeV", &En_MeV);
  tree->SetBranchAddress("x_mm", &x_mm);
  tree->SetBranchAddress("y_mm", &y_mm);
  tree->SetBranchAddress("z_mm", &z_mm);
  tree->SetBranchAddress("cosTheta", &cosTheta);
  tree->SetBranchAddress("weight", &weight);
  tree->SetBranchAddress("time_ns", &time_ns);
  tree->SetBranchAddress("surface_id", &surface_id);
  if (tree->GetBranch("surface_name")) {
    tree->SetBranchAddress("surface_name", &surface_name);
  }

  std::ofstream os(outPath);
  os << "event_id,En_MeV,x_mm,y_mm,z_mm,cosTheta,weight,time_ns,surface_id,surface_name\n";
  tree->Print();
  const Long64_t n = tree->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    if (!SafeGetEntry(tree, i, "NeutronSurf(csv)")) break;
    os << event_id << "," << En_MeV << "," << x_mm << "," << y_mm << "," << z_mm << "," << cosTheta << "," << weight << ","
       << time_ns << "," << surface_id << "," << (surface_name ? *surface_name : "") << "\n";
  }
}



double ReadPrimaryElectronCount(TFile* file) {
  auto* tree = dynamic_cast<TTree*>(file->Get("RunMeta"));
  if (!tree) return 0.0;
  PrepareTreeForSafeRead(tree);
  EnableBranches(tree, {"nEvents"});
  int nEvents = 0;
  if (!tree->GetBranch("nEvents")) return 0.0;
  tree->SetBranchAddress("nEvents", &nEvents);
  if (!SafeGetEntry(tree, 0, "RunMeta(nEvents)")) return 0.0;
  return static_cast<double>(nEvents);
}

void ExportPhotonSourceDataAndSpectra(TFile* file, const std::string& outDir) {
  auto* tree = dynamic_cast<TTree*>(file->Get("PhotonSurf"));
  if (!tree) return;
  PrepareTreeForSafeRead(tree);
  EnableBranches(tree, {"event_id", "E_MeV", "x_mm", "y_mm", "z_mm", "cosTheta", "weight", "time_ns", "surface_id"});

  int event_id = 0;
  double E_MeV = 0.0;
  double x_mm = 0.0;
  double y_mm = 0.0;
  double z_mm = 0.0;
  double cosTheta = 0.0;
  double weight = 1.0;
  double time_ns = 0.0;
  int surface_id = 0;
  tree->SetBranchAddress("event_id", &event_id);
  tree->SetBranchAddress("E_MeV", &E_MeV);
  tree->SetBranchAddress("x_mm", &x_mm);
  tree->SetBranchAddress("y_mm", &y_mm);
  tree->SetBranchAddress("z_mm", &z_mm);
  tree->SetBranchAddress("cosTheta", &cosTheta);
  tree->SetBranchAddress("weight", &weight);
  tree->SetBranchAddress("time_ns", &time_ns);
  tree->SetBranchAddress("surface_id", &surface_id);

  auto* h_linear = new TH1D("h_photon_source_spectrum_linear", "Photon source spectrum (linear)", 200, 0.0, 100.0);
  auto* h_log = new TH1D("h_photon_source_spectrum_log", "Photon source spectrum (log)", 200, 1e-6, 100.0);
  auto* h_focus = new TH1D("h_photon_source_spectrum_4p5_30", "Photon source spectrum (4.5-30 MeV)", 255, 4.5, 30.0);
  h_log->GetXaxis()->SetMoreLogLabels(true);

  tree->Print();
  const Long64_t n = tree->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    if (!SafeGetEntry(tree, i, "PhotonSurf(source_export)")) break;
    if (E_MeV >= 0.0) {
      h_linear->Fill(E_MeV, weight);
    }
    if (E_MeV > 0.0) {
      h_log->Fill(E_MeV, weight);
    }
    if (E_MeV >= 4.5 && E_MeV <= 30.0) {
      h_focus->Fill(E_MeV, weight);
    }
  }

  TCanvas c1("c_ph_lin", "c_ph_lin", 1000, 800);
  h_linear->GetXaxis()->SetTitle("Photon energy E_{#gamma} (MeV)");
  h_linear->GetYaxis()->SetTitle("Counts (weighted)");
  h_linear->Draw("HIST");
  c1.SaveAs((outDir + "/photon_source_spectrum_linear.png").c_str());

  TCanvas c2("c_ph_log", "c_ph_log", 1000, 800);
  c2.SetLogx();
  h_log->GetXaxis()->SetTitle("Photon energy E_{#gamma} (MeV)");
  h_log->GetYaxis()->SetTitle("Counts (weighted)");
  h_log->Draw("HIST");
  c2.SaveAs((outDir + "/photon_source_spectrum_log.png").c_str());

  TCanvas c3("c_ph_focus", "c_ph_focus", 1000, 800);
  h_focus->GetXaxis()->SetTitle("Photon energy E_{#gamma} (MeV)");
  h_focus->GetYaxis()->SetTitle("Counts (weighted)");
  h_focus->Draw("HIST");
  c3.SaveAs((outDir + "/photon_source_spectrum_4p5_30.png").c_str());
}

void WriteParticleYieldsPerElectron(TFile* file, const std::string& outPath) {
  const double nElectrons = ReadPrimaryElectronCount(file);
  if (nElectrons <= 0.0) return;

  double neutronWeighted = 0.0;
  double photonWeighted = 0.0;
  Long64_t neutronEntries = 0;
  Long64_t photonEntries = 0;
  long long nNeutronSummary = -1;
  long long nGammaAbove5MeVSummary = -1;
  long long nNeutronModelExitSummary = -1;

  if (auto* ntree = dynamic_cast<TTree*>(file->Get("NeutronSurf"))) {
    PrepareTreeForSafeRead(ntree);
    EnableBranches(ntree, {"weight"});
    double weight = 1.0;
    if (ntree->GetBranch("weight")) {
      ntree->SetBranchAddress("weight", &weight);
      neutronEntries = ntree->GetEntries();
      for (Long64_t i = 0; i < neutronEntries; ++i) {
        if (!SafeGetEntry(ntree, i, "NeutronSurf(yield)")) break;
        neutronWeighted += weight;
      }
    }
  }

  if (auto* ptree = dynamic_cast<TTree*>(file->Get("PhotonSurf"))) {
    PrepareTreeForSafeRead(ptree);
    EnableBranches(ptree, {"weight"});
    double weight = 1.0;
    if (ptree->GetBranch("weight")) {
      ptree->SetBranchAddress("weight", &weight);
      photonEntries = ptree->GetEntries();
      for (Long64_t i = 0; i < photonEntries; ++i) {
        if (!SafeGetEntry(ptree, i, "PhotonSurf(yield)")) break;
        photonWeighted += weight;
      }
    }
  }


  if (auto* sumTree = dynamic_cast<TTree*>(file->Get("run_summary"))) {
    PrepareTreeForSafeRead(sumTree);
    EnableBranches(sumTree, {"nNeutron", "nGammaAbove5MeV", "nNeutronModelExit"});
    ReadNumericRunSummaryValue(sumTree, "nNeutron", nNeutronSummary);
    ReadNumericRunSummaryValue(sumTree, "nGammaAbove5MeV", nGammaAbove5MeVSummary);
    ReadNumericRunSummaryValue(sumTree, "nNeutronModelExit", nNeutronModelExitSummary);
  }

  std::ofstream os(outPath);
  os << "{\n";
  os << "  \"n_primary_electrons\": " << nElectrons << ",\n";
  os << "  \"neutron_entries\": " << neutronEntries << ",\n";
  os << "  \"photon_entries\": " << photonEntries << ",\n";
  os << "  \"neutrons_per_electron\": " << (static_cast<double>(neutronEntries) / nElectrons) << ",\n";
  os << "  \"photons_per_electron\": " << (static_cast<double>(photonEntries) / nElectrons) << ",\n";
  if (nNeutronSummary >= 0) {
    os << "  \"neutrons_per_electron_from_run_summary\": " << (static_cast<double>(nNeutronSummary) / nElectrons) << ",\n";
  }
  if (nGammaAbove5MeVSummary >= 0) {
    os << "  \"photons_above5MeV_per_electron_from_run_summary\": " << (static_cast<double>(nGammaAbove5MeVSummary) / nElectrons) << ",\n";
  }
  if (nNeutronModelExitSummary >= 0) {
    os << "  \"neutrons_model_exit_per_electron_from_run_summary\": " << (static_cast<double>(nNeutronModelExitSummary) / nElectrons) << ",\n";
  }
  os << "  \"neutrons_weighted_per_electron\": " << (neutronWeighted / nElectrons) << ",\n";
  os << "  \"photons_weighted_per_electron\": " << (photonWeighted / nElectrons) << "\n";
  os << "}\n";
}
void ExportNeutronSourceDataAndSpectra(TFile* file, const std::string& outDir) {
  auto* tree = dynamic_cast<TTree*>(file->Get("NeutronSurf"));
  if (!tree) return;
  PrepareTreeForSafeRead(tree);
  EnableBranches(tree, {"event_id", "En_MeV", "x_mm", "y_mm", "z_mm", "cosTheta", "weight", "time_ns", "surface_id"});
  int event_id = 0;
  double En_MeV = 0.0;
  double x_mm = 0.0;
  double y_mm = 0.0;
  double z_mm = 0.0;
  double cosTheta = 0.0;
  double weight = 1.0;
  double time_ns = 0.0;
  int surface_id = 0;
  tree->SetBranchAddress("event_id", &event_id);
  tree->SetBranchAddress("En_MeV", &En_MeV);
  tree->SetBranchAddress("x_mm", &x_mm);
  tree->SetBranchAddress("y_mm", &y_mm);
  tree->SetBranchAddress("z_mm", &z_mm);
  tree->SetBranchAddress("cosTheta", &cosTheta);
  tree->SetBranchAddress("weight", &weight);
  tree->SetBranchAddress("time_ns", &time_ns);
  tree->SetBranchAddress("surface_id", &surface_id);

  auto* h_linear = new TH1D("h_neutron_source_spectrum_linear", "Neutron source spectrum (linear)", 250, 0.0, 5.0);
  auto* h_log = new TH1D("h_neutron_source_spectrum_log", "Neutron source spectrum (log)", 250, 2.5e-9, 5.0);
  h_log->GetXaxis()->SetMoreLogLabels(true);

  std::ofstream os(outDir + "/neutron_source.csv");
  os << "event_id,En_MeV,x_mm,y_mm,z_mm,cosTheta,weight,time_ns,surface_id\n";
  tree->Print();
  const Long64_t n = tree->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    if (!SafeGetEntry(tree, i, "NeutronSurf(source_export)")) break;
    h_linear->Fill(En_MeV, weight);
    if (En_MeV > 0.0) {
      h_log->Fill(En_MeV, weight);
    }
    os << event_id << "," << En_MeV << "," << x_mm << "," << y_mm << "," << z_mm << "," << cosTheta << "," << weight << ","
       << time_ns << "," << surface_id << "\n";
  }

  TCanvas c1("c1", "c1", 1000, 800);
  h_linear->GetXaxis()->SetTitle("Neutron energy En (MeV)");
  h_linear->GetYaxis()->SetTitle("Counts (weighted)");
  h_linear->Draw("HIST");
  c1.SaveAs((outDir + "/neutron_source_spectrum_linear.png").c_str());

  TCanvas c2("c2", "c2", 1000, 800);
  c2.SetLogx();
  h_log->GetXaxis()->SetTitle("Neutron energy En (MeV)");
  h_log->GetYaxis()->SetTitle("Counts (weighted)");
  h_log->Draw("HIST");
  c2.SaveAs((outDir + "/neutron_source_spectrum_log.png").c_str());
}
} // namespace

void export_run_artifacts(const char* rootPath, const char* outBase = "results/root/png") {
  std::unique_ptr<TFile> file(TFile::Open(rootPath, "READ"));
  if (!file || file->IsZombie()) {
    ::Error("export_run_artifacts", "Failed to open ROOT file: %s", rootPath);
    return;
  }

  const std::string targetType = ReadTargetType(file.get());
  const std::string timestamp = TimestampFromFile(file.get(), rootPath);
  const std::string runLabel = Sanitize(timestamp + "_" + targetType);
  const std::string outDir = std::string(outBase) + "/" + runLabel;
  gSystem->mkdir(outDir.c_str(), true);

  WriteRunMetaJson(file.get(), outDir + "/runmeta.json");
  WriteNeutronSurfCsv(file.get(), outDir + "/neutron_surface.csv");
  ExportNeutronSourceDataAndSpectra(file.get(), outDir);
  ExportPhotonSourceDataAndSpectra(file.get(), outDir);
  WriteParticleYieldsPerElectron(file.get(), outDir + "/particle_yields_per_electron.json");

  gStyle->SetOptStat(0);

  TIter next(file->GetListOfKeys());
  while (auto* key = dynamic_cast<TKey*>(next())) {
    TObject* obj = key->ReadObj();
    if (!obj) continue;
    if (obj->InheritsFrom("TH3")) {
      auto* h3 = dynamic_cast<TH3*>(obj);
      if (!h3) continue;
      TCanvas c("c", "c", 1200, 900);
      c.SetLeftMargin(0.12);
      c.SetRightMargin(0.20);
      c.SetBottomMargin(0.12);
      c.SetTopMargin(0.08);
      auto* projXY = h3->Project3D("xy");
      projXY->SetStats(0);
      projXY->Draw("COLZ");
      c.SaveAs((outDir + "/" + std::string(h3->GetName()) + "_xy.png").c_str());
      auto* projXZ = h3->Project3D("xz");
      projXZ->SetStats(0);
      projXZ->Draw("COLZ");
      c.SaveAs((outDir + "/" + std::string(h3->GetName()) + "_xz.png").c_str());
      auto* projYZ = h3->Project3D("yz");
      projYZ->SetStats(0);
      projYZ->Draw("COLZ");
      c.SaveAs((outDir + "/" + std::string(h3->GetName()) + "_yz.png").c_str());
    } else if (obj->InheritsFrom("TH2")) {
      auto* h2 = dynamic_cast<TH2*>(obj);
      if (!h2) continue;
      TCanvas c("c", "c", 1200, 900);
      c.SetLeftMargin(0.12);
      c.SetRightMargin(0.20);
      c.SetBottomMargin(0.12);
      c.SetTopMargin(0.08);
      h2->SetStats(0);
      h2->Draw("COLZ");
      c.SaveAs((outDir + "/" + std::string(h2->GetName()) + ".png").c_str());
    } else if (obj->InheritsFrom("TH1")) {
      auto* h1 = dynamic_cast<TH1*>(obj);
      if (!h1) continue;
      TCanvas c("c", "c", 1200, 900);
      c.SetLeftMargin(0.12);
      c.SetRightMargin(0.08);
      c.SetBottomMargin(0.12);
      c.SetTopMargin(0.08);
      h1->SetStats(0);
      h1->Draw("HIST");
      c.SaveAs((outDir + "/" + std::string(h1->GetName()) + ".png").c_str());
    }
  }
  ::Info("export_run_artifacts", "Artifacts saved to %s", outDir.c_str());
}
