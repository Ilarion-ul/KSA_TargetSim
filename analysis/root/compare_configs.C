// compare_configs.C
// Usage:
//   root -l -q 'analysis/root/compare_configs.C("results/root/run_WTa.root","results/root/run_UAl.root")'
// Note: compare normalized metrics per primary e- where possible.

#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include <iostream>

void compare_configs(const char* fileA = "results/root/run_WTa.root", const char* fileB = "results/root/run_UAl.root") {
  gSystem->mkdir("results/vis", true);

  auto readSummary = [](const char* fn, double& sub, double& coat) -> bool {
    TFile f(fn, "READ");
    if (f.IsZombie()) return false;
    auto* t = dynamic_cast<TTree*>(f.Get("run_summary"));
    if (!t) return false;
    t->SetBranchAddress("edep_substrate", &sub);
    t->SetBranchAddress("edep_coating", &coat);
    t->GetEntry(0);
    return true;
  };

  double subA = 0.0, coatA = 0.0, subB = 0.0, coatB = 0.0;
  if (!readSummary(fileA, subA, coatA) || !readSummary(fileB, subB, coatB)) {
    std::cerr << "Failed to read one or both run_summary trees." << std::endl;
    return;
  }

  auto* c = new TCanvas("c_cmp", "WTa vs UAl", 700, 500);
  auto* h = new TH1D("hCmp", "Edep comparison (per run);Case/Region;MeV", 4, 0.5, 4.5);
  h->GetXaxis()->SetBinLabel(1, "WTa Sub");
  h->GetXaxis()->SetBinLabel(2, "WTa Coat");
  h->GetXaxis()->SetBinLabel(3, "UAl Sub");
  h->GetXaxis()->SetBinLabel(4, "UAl Coat");
  h->SetBinContent(1, subA);
  h->SetBinContent(2, coatA);
  h->SetBinContent(3, subB);
  h->SetBinContent(4, coatB);
  h->Draw("BAR");

  c->SaveAs("results/vis/compare.png");
}
