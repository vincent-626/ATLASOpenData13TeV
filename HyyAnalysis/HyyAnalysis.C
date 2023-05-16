// H to yy Analysis using RDataFrame
// Run the code with root -b -l -q HyyAnalysis.C

#include <iostream>
#include <TROOT.h>

using namespace ROOT;

float ComputeInvMass(RVecF, RVecF, RVecF, RVecF);
void Analyze(TString, TString);

void HyyAnalysis() {
    TString inDir = "/Users/vincent626/Code/CERN/opendata/ATLASOpenData13TeV/Samples/GamGam/";

    // Data Samples
    Analyze(inDir + "Data/data_A.GamGam.root", "data_A.GamGam");
    Analyze(inDir + "Data/data_B.GamGam.root", "data_B.GamGam");
    Analyze(inDir + "Data/data_C.GamGam.root", "data_C.GamGam");
    Analyze(inDir + "Data/data_D.GamGam.root", "data_D.GamGam");

    // MC Samples
    Analyze(inDir + "MC/mc_343981.ggH125_gamgam.GamGam.root", "ggH125_gamgam");
    Analyze(inDir + "MC/mc_345041.VBFH125_gamgam.GamGam.root", "VBFH125_gamgam");
    Analyze(inDir + "MC/mc_345318.WpH125J_Wincl_gamgam.GamGam.root", "WpH125J_Wincl_gamgam");
    Analyze(inDir + "MC/mc_345319.ZH125J_Zincl_gamgam.GamGam.root", "ZH125J_Zincl_gamgam");
    Analyze(inDir + "MC/mc_341081.ttH125_gamgam.GamGam.root", "ttH125_gamgam");
}

void Analyze(TString inFile, TString sampleName) {
    // Parallel Processing
    EnableImplicitMT();

    // Printing Status
    std::cout << "Analyzing " << sampleName << "..." << std::endl;

    // Construct RDataFrame
    TFile *f = TFile::Open(inFile);
    RDataFrame df("mini", f);

    // Define weight
    auto weighted = df.Define("weight", "1.0");

    if (!sampleName.Contains("data")) {
        weighted = df.Define("weight", "scaleFactor_PHOTON*scaleFactor_PhotonTRIGGER*scaleFactor_PILEUP * mcWeight");
    }

    // Good Photon Selection + 2 Photons Filtering
    auto good_gammas = weighted.Filter("trigP")
                            .Define("good_photon", "photon_isTightID && photon_pt > 25000 && abs(photon_eta) < 2.37 && (abs(photon_eta) < 1.37 || abs(photon_eta) > 1.52)")
                            .Filter("Sum(good_photon) == 2");
    
    // Isolated Photons
    auto iso = good_gammas.Filter("Sum(photon_ptcone30[good_photon] / photon_pt[good_photon] < 0.065) == 2")
                            .Filter("Sum(photon_etcone20[good_photon] / photon_pt[good_photon] < 0.065) == 2");

    // Kinematic Cuts
    auto selected = iso.Define("m_yy", "ComputeInvMass(photon_pt[good_photon], photon_eta[good_photon], photon_phi[good_photon], photon_E[good_photon])")
                        .Filter("photon_pt[good_photon][0] / 1000. / m_yy > 0.35")
                        .Filter("photon_pt[good_photon][1] / 1000. / m_yy > 0.25")
                        .Filter("m_yy > 105 && m_yy < 160");

    // Filling Histogram
    auto histMass = selected.Histo1D({"mass", "m_yy", 55u, 105., 160.}, "m_yy", "weight");

    // Save Histogram
    TString output_name = sampleName + ".root";
    TFile physicsOutput(output_name, "recreate");
    histMass->Write();
    physicsOutput.Close();
}

float ComputeInvMass(RVecF pt, RVecF eta, RVecF phi, RVecF E) {
    Math::PtEtaPhiEVector photon_1 = Math::PtEtaPhiEVector(pt[0], eta[0], phi[0], E[0]);
    Math::PtEtaPhiEVector photon_2 = Math::PtEtaPhiEVector(pt[1], eta[1], phi[1], E[1]);
    
    auto photon = photon_1 + photon_2;
    
    return photon.mass()/1000.;
}