// plot_run_summary.C
// Usage:
//   root -l -q 'analysis/root/plot_run_summary.C("results/root/run_WTa.root")'

#include <TCanvas.h>
#include <TH1D.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include <iostream>

void plot_run_summary(const char* fileName = "results/root/run_WTa.root") {
  gSystem->mkdir("results/vis", true);

  TFile f(fileName, "READ");
  if (f.IsZombie()) {
    std::cerr << "Cannot open file: " << fileName << std::endl;
    return;
  }

  auto* t = dynamic_cast<TTree*>(f.Get("run_summary"));
  if (!t) {
    std::cerr << "Tree run_summary not found. Per-event tree is optional; macro reads run summary only." << std::endl;
    return;
  }

  double edepSub = 0.0, edepCoat = 0.0;
  long long nGamma = 0, nNeutron = 0;
  t->SetBranchAddress("edep_substrate", &edepSub);
  t->SetBranchAddress("edep_coating", &edepCoat);
  t->SetBranchAddress("nGamma", &nGamma);
  t->SetBranchAddress("nNeutron", &nNeutron);
  t->GetEntry(0);

  std::cout << "Edep substrate [MeV]: " << edepSub << "\n"
            << "Edep coating   [MeV]: " << edepCoat << "\n"
            << "nGamma: " << nGamma << ", nNeutron: " << nNeutron << std::endl;

  auto* c = new TCanvas("c_summary", "Run Summary", 900, 450);
  c->Divide(2, 1);

  c->cd(1);
  auto* hEdep = new TH1D("hEdep", "Energy deposition;Region;MeV", 2, 0.5, 2.5);
  hEdep->GetXaxis()->SetBinLabel(1, "Substrate");
  hEdep->GetXaxis()->SetBinLabel(2, "Coating");
  hEdep->SetBinContent(1, edepSub);
  hEdep->SetBinContent(2, edepCoat);
  hEdep->Draw("BAR");

  c->cd(2);
  auto* hN = new TH1D("hN", "Particle counts;Type;Count", 2, 0.5, 2.5);
  hN->GetXaxis()->SetBinLabel(1, "#gamma");
  hN->GetXaxis()->SetBinLabel(2, "n");
  hN->SetBinContent(1, nGamma);
  hN->SetBinContent(2, nNeutron);
  hN->Draw("BAR");

  c->SaveAs("results/vis/summary.png");
}
