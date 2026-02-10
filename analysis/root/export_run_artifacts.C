#include <TCanvas.h>
#include <TDatime.h>
#include <TFile.h>
#include <TKey.h>
#include <TObjString.h>
#include <TObject.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace {
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
  std::string target_type;
  tree->SetBranchAddress("target_type", &target_type);
  tree->GetEntry(0);
  return target_type.empty() ? "unknown" : target_type;
}

void WriteRunMetaJson(TFile* file, const std::string& outPath) {
  auto* tree = dynamic_cast<TTree*>(file->Get("RunMeta"));
  if (!tree) return;
  std::string target_type;
  std::string geometry_json;
  double beam_energy_MeV = 0.0;
  double beam_energy_sigma_rel = 0.0;
  double beam_sigma_x_mm = 0.0;
  double beam_sigma_y_mm = 0.0;
  double beam_sigma_theta_x_mrad = 0.0;
  double beam_sigma_theta_y_mrad = 0.0;
  std::string pulse_mode;
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
  tree->GetEntry(0);

  std::ofstream os(outPath);
  os << "{\n";
  os << "  \"target_type\": \"" << target_type << "\",\n";
  os << "  \"geometry\": " << geometry_json << ",\n";
  os << "  \"beam\": {\n";
  os << "    \"energy_MeV\": " << beam_energy_MeV << ",\n";
  os << "    \"energy_sigma_rel\": " << beam_energy_sigma_rel << ",\n";
  os << "    \"sigma_x_mm\": " << beam_sigma_x_mm << ",\n";
  os << "    \"sigma_y_mm\": " << beam_sigma_y_mm << ",\n";
  os << "    \"sigma_theta_x_mrad\": " << beam_sigma_theta_x_mrad << ",\n";
  os << "    \"sigma_theta_y_mrad\": " << beam_sigma_theta_y_mrad << "\n";
  os << "  },\n";
  os << "  \"pulse\": {\n";
  os << "    \"mode\": \"" << pulse_mode << "\",\n";
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
  int event_id = 0;
  double En_MeV = 0.0;
  double x_mm = 0.0;
  double y_mm = 0.0;
  double z_mm = 0.0;
  double cosTheta = 0.0;
  double weight = 1.0;
  double time_ns = 0.0;
  int surface_id = 0;
  std::string surface_name;
  tree->SetBranchAddress("event_id", &event_id);
  tree->SetBranchAddress("En_MeV", &En_MeV);
  tree->SetBranchAddress("x_mm", &x_mm);
  tree->SetBranchAddress("y_mm", &y_mm);
  tree->SetBranchAddress("z_mm", &z_mm);
  tree->SetBranchAddress("cosTheta", &cosTheta);
  tree->SetBranchAddress("weight", &weight);
  tree->SetBranchAddress("time_ns", &time_ns);
  tree->SetBranchAddress("surface_id", &surface_id);
  tree->SetBranchAddress("surface_name", &surface_name);
  std::ofstream os(outPath);
  os << "event_id,En_MeV,x_mm,y_mm,z_mm,cosTheta,weight,time_ns,surface_id,surface_name\n";
  const Long64_t n = tree->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    tree->GetEntry(i);
    os << event_id << "," << En_MeV << "," << x_mm << "," << y_mm << "," << z_mm << "," << cosTheta << "," << weight << ","
       << time_ns << "," << surface_id << "," << surface_name << "\n";
  }
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

  TIter next(file->GetListOfKeys());
  while (auto* key = dynamic_cast<TKey*>(next())) {
    TObject* obj = key->ReadObj();
    if (!obj) continue;
    if (obj->InheritsFrom("TH3")) {
      auto* h3 = dynamic_cast<TH3*>(obj);
      if (!h3) continue;
      TCanvas c("c", "c", 1000, 800);
      auto* projXY = h3->Project3D("xy");
      projXY->Draw("COLZ");
      c.SaveAs((outDir + "/" + std::string(h3->GetName()) + "_xy.png").c_str());
      auto* projXZ = h3->Project3D("xz");
      projXZ->Draw("COLZ");
      c.SaveAs((outDir + "/" + std::string(h3->GetName()) + "_xz.png").c_str());
      auto* projYZ = h3->Project3D("yz");
      projYZ->Draw("COLZ");
      c.SaveAs((outDir + "/" + std::string(h3->GetName()) + "_yz.png").c_str());
    } else if (obj->InheritsFrom("TH2")) {
      auto* h2 = dynamic_cast<TH2*>(obj);
      if (!h2) continue;
      TCanvas c("c", "c", 1000, 800);
      h2->Draw("COLZ");
      c.SaveAs((outDir + "/" + std::string(h2->GetName()) + ".png").c_str());
    } else if (obj->InheritsFrom("TH1")) {
      auto* h1 = dynamic_cast<TH1*>(obj);
      if (!h1) continue;
      TCanvas c("c", "c", 1000, 800);
      h1->Draw("HIST");
      c.SaveAs((outDir + "/" + std::string(h1->GetName()) + ".png").c_str());
    }
  }
  ::Info("export_run_artifacts", "Artifacts saved to %s", outDir.c_str());
}
