#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <ROOT/RVec.hxx>
#include <cmath>

void z_boson_resonance() {

    // Abrir archivo ROOT
   TFile *file = TFile::Open(
    "/mnt/c/Users/jonat/Downloads/ODEO_FEB2025_v0_2to4lep_data15_periodD.2to4lep.root"
);


    if (!file || file->IsZombie()) {
        Error("z_boson_resonance", "No se pudo abrir el archivo");
        return;
    }

    // Obtener el TTree
    TTree *tree = (TTree*) file->Get("analysis");

    if (!tree) {
        Error("z_boson_resonance", "No se encontró el TTree 'analysis'");
        return;
    }

    // Histograma de masa invariante
    TH1F *h_mll = new TH1F(
        "h_mll",
        "Z #rightarrow l^{+}l^{-}; m_{ll} [GeV]; Events",
        100, 60, 120
    );

    // Variables
    int lep_n;
    ROOT::VecOps::RVec<float> *lep_pt = nullptr;
    ROOT::VecOps::RVec<float> *lep_eta = nullptr;
    ROOT::VecOps::RVec<float> *lep_phi = nullptr;
    ROOT::VecOps::RVec<int>   *lep_type = nullptr;
    ROOT::VecOps::RVec<bool>  *lep_isTightIso = nullptr;
    ROOT::VecOps::RVec<bool>  *lep_isTrigMatched = nullptr;

    // Branches
    tree->SetBranchAddress("lep_n", &lep_n);
    tree->SetBranchAddress("lep_pt", &lep_pt);
    tree->SetBranchAddress("lep_eta", &lep_eta);
    tree->SetBranchAddress("lep_phi", &lep_phi);
    tree->SetBranchAddress("lep_type", &lep_type);
    tree->SetBranchAddress("lep_isTightIso", &lep_isTightIso);
    tree->SetBranchAddress("lep_isTrigMatched", &lep_isTrigMatched);

    // Loop de eventos
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        if (lep_n != 2) continue;

        if (!lep_isTightIso->at(0)) continue;
        if (!lep_isTightIso->at(1)) continue;

        if (!(lep_isTrigMatched->at(0) || lep_isTrigMatched->at(1))) continue;

        // mismo sabor (ee o μμ)
        if (lep_type->at(0) != lep_type->at(1)) continue;

        float deta = lep_eta->at(0) - lep_eta->at(1);
        float dphi = lep_phi->at(0) - lep_phi->at(1);

        if (dphi > M_PI)  dphi -= 2*M_PI;
        if (dphi < -M_PI) dphi += 2*M_PI;

        float mll = std::sqrt(
            2 * lep_pt->at(0) * lep_pt->at(1) *
            (std::cosh(deta) - std::cos(dphi))
        );

        h_mll->Fill(mll);
    }

    // Dibujar histograma
    h_mll->Draw();

    // Ajuste Voigt (Breit–Wigner convolucionada con Gauss)
    TF1 *f_voigt = new TF1(
        "f_voigt",
        "[0]*TMath::Voigt(x-[1],[2],[3])",
        80, 100
    );

    f_voigt->SetParameters(5000, 91.2, 2.0, 2.5);
    f_voigt->SetParLimits(1, 90.0, 92.0);
    f_voigt->SetParLimits(2, 0.5, 5.0);
    f_voigt->SetParLimits(3, 1.0, 5.0);

    h_mll->Fit(f_voigt, "R");
    f_voigt->Draw("same");
}
