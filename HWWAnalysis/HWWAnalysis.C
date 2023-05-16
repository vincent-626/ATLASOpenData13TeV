// H to Ww Analysis using RDataFrame
// Run the code with root -b -l -q HWWAnalysis.C

#include <iostream>
#include <TROOT.h>

using namespace ROOT;

void Analyze(TString, TString);
RVecI GoodLepton(RVecB, RVecF, RVecF, RVecF, RVecF, RVecF, RVecF, RVecI, RVecF, RVecF, RVecF);
float ComputeInvMass(RVecF, RVecF, RVecF, RVecF);
float ComputePhiLL(RVecF, RVecF, RVecF, RVecF);
float ComputeDPhiLL(RVecF, RVecF, RVecF, RVecF);
float ComputePtLL(RVecF, RVecF, RVecF, RVecF);
float ComputeMt(RVecF, RVecF, RVecF, RVecF, float, float);

void HWWAnalysis() {
    TString inDir = "/Users/vincent626/Code/CERN/opendata/ATLASOpenData13TeV/Samples/2lep/";

    // Data Samples
    Analyze(inDir + "Data/data_A.2lep.root", "data_A.2lep");
    Analyze(inDir + "Data/data_B.2lep.root", "data_B.2lep");
    Analyze(inDir + "Data/data_C.2lep.root", "data_C.2lep");
    Analyze(inDir + "Data/data_D.2lep.root", "data_D.2lep");

    // MC Samples
    // Diboson
    Analyze(inDir + "MC/mc_363356.ZqqZll.2lep.root", "ZqqZll.2lep");
    Analyze(inDir + "MC/mc_363358.WqqZll.2lep.root", "WqqZll.2lep");
    Analyze(inDir + "MC/mc_363359.WpqqWmlv.2lep.root", "WpqqWmlv.2lep");
    Analyze(inDir + "MC/mc_363360.WplvWmqq.2lep.root", "WplvWmqq.2lep");
    Analyze(inDir + "MC/mc_363489.WlvZqq.2lep.root", "WlvZqq.2lep");

    Analyze(inDir + "MC/mc_363490.llll.2lep.root", "llll.2lep");
    Analyze(inDir + "MC/mc_363491.lllv.2lep.root", "lllv.2lep");
    Analyze(inDir + "MC/mc_363492.llvv.2lep.root", "llvv.2lep");
    Analyze(inDir + "MC/mc_363493.lvvv.2lep.root", "lvvv.2lep");

    //Single Top
    Analyze(inDir + "MC/mc_410011.single_top_tchan.2lep.root", "single_top_tchan.2lep");
    Analyze(inDir + "MC/mc_410012.single_antitop_tchan.2lep.root", "single_antitop_tchan.2lep");
    Analyze(inDir + "MC/mc_410025.single_top_schan.2lep.root", "single_top_schan.2lep");
    Analyze(inDir + "MC/mc_410026.single_antitop_schan.2lep.root", "single_antitop_schan.2lep");
    Analyze(inDir + "MC/mc_410013.single_top_wtchan.2lep.root", "single_top_wtchan.2lep");
    Analyze(inDir + "MC/mc_410014.single_antitop_wtchan.2lep.root", "single_antitop_wtchan.2lep");

    // Z+Jets Inclusive
    Analyze(inDir + "MC/mc_361106.Zee.2lep.root", "Zee.2lep");
    Analyze(inDir + "MC/mc_361107.Zmumu.2lep.root", "Zmumu.2lep");

    Analyze(inDir + "MC/mc_361108.Ztautau.2lep.root", "Ztautau.2lep");

    // ttbar
    Analyze(inDir + "MC/mc_410000.ttbar_lep.2lep.root", "ttbar_lep.2lep");

    // W+Jets Inclusive
    Analyze(inDir + "MC/mc_361100.Wplusenu.2lep.root", "Wplusenu.2lep");
    Analyze(inDir + "MC/mc_361101.Wplusmunu.2lep.root", "Wplusmunu.2lep");
    Analyze(inDir + "MC/mc_361102.Wplustaunu.2lep.root", "Wplustaunu.2lep");
    Analyze(inDir + "MC/mc_361103.Wminusenu.2lep.root", "Wminusenu.2lep");
    Analyze(inDir + "MC/mc_361104.Wminusmunu.2lep.root", "Wminusmunu.2lep");
    Analyze(inDir + "MC/mc_361105.Wminustaunu.2lep.root", "Wminustaunu.2lep");

    //ggH
    Analyze(inDir + "MC/mc_345324.ggH125_WW2lep.2lep.root", "ggH125_WW2lep.2lep");

    //VBF
    Analyze(inDir + "MC/mc_345323.VBFH125_WW2lep.2lep.root", "VBFH125_WW2lep.2lep");
}

void Analyze(TString inFile, TString sampleName) {
    // Parallel Processing
    EnableImplicitMT();

    // Printing Status
    std::cout << "Analyzing " << sampleName << "..." << endl;

    // Construct RDataFrame
    TFile *f = TFile::Open(inFile);
    RDataFrame df("mini", f);

    // Define Weight
    auto weighted = df.Define("weight", "1.0");

    if (!sampleName.Contains("data")) {
        weighted = df.Define("weight", "scaleFactor_ELE * scaleFactor_MUON * scaleFactor_LepTRIGGER * scaleFactor_PILEUP * mcWeight");
    }

    ////////////////////////////////////////////////
    // Event Selection
    ////////////////////////////////////////////////

    // This example shows an one-liner event selection
    // The comments below explain the RDataFrame operations line by line

    // Electron/Muon Trigger
    // Good Lepton (Tight ID, Tight Isolation, pT > 15 GeV and some criteria for electron and muon)
    // 2 Good Leptons (Different Flavour, Opposite Sign, pT1 > 22 GeV)
    // Change MET to GeV
    // MET > 20 GeV     
    // Good jet
    // 0/1 jet
    // Good B-jet
    // 0 b-jet
    // del phi(ll, MET) > pi/2
    // Define del phi(ll)
    // del phi(ll) < 1.8
    // Define pT(ll)
    // pT(ll) > 30 GeV
    // Define mll
    // 10 GeV < mll < 55 GeV
    // Define mt

    auto selected = weighted.Filter("trigE || trigM")
                    .Define("good_lep", "GoodLepton(lep_isTightID, lep_ptcone30, lep_etcone20, lep_pt, lep_eta, lep_trackd0pvunbiased, lep_tracksigd0pvunbiased, lep_type, lep_z0, lep_phi, lep_E)")
                    .Filter("Sum(good_lep) == 2 && lep_type[good_lep][0] != lep_type[good_lep][1] && lep_charge[good_lep][0] * lep_charge[good_lep][1] < 0. && lep_pt[good_lep][0] > 22000.")
                    .Redefine("met_et", "met_et/1000.")
                    .Filter("met_et > 20.")
                    .Define("good_jet", "jet_pt > 30000. && abs(jet_eta) < 2.5")
                    .Filter("(Sum(good_jet) == 0 || Sum(good_jet) == 1)")
                    .Define("good_bjet", "jet_MV2c10[good_jet] > 0.1758475 && jet_pt[good_jet] > 20000.")
                    .Filter("Sum(good_bjet) == 0")
                    .Filter("TMath::Abs(ComputePhiLL(lep_pt[good_lep], lep_eta[good_lep], lep_phi[good_lep], lep_E[good_lep]) - met_phi) > TMath::Pi()/2.")
                    .Define("ll_dphi", "TMath::Abs(ComputeDPhiLL(lep_pt[good_lep], lep_eta[good_lep], lep_phi[good_lep], lep_E[good_lep]))")
                    .Filter("ll_dphi < 1.8")
                    .Define("ll_pt", "ComputePtLL(lep_pt[good_lep], lep_eta[good_lep], lep_phi[good_lep], lep_E[good_lep])")
                    .Filter("ll_pt > 30.")
                    .Define("ll_m", "ComputeInvMass(lep_pt[good_lep], lep_eta[good_lep], lep_phi[good_lep], lep_E[good_lep])")
                    .Filter("ll_m > 10. && ll_m < 55.")
                    .Define("mt", "ComputeMt(lep_pt[good_lep], lep_eta[good_lep], lep_phi[good_lep], lep_E[good_lep], met_et, met_phi)");

    auto c = selected.Count();
    std::cout << "Number of Events: " << *c << std::endl;


    ////////////////////////////////////////////////
    // Histogramming
    ////////////////////////////////////////////////

    // Fill Histograms
    auto histDPhi = selected.Histo1D({"dphi", "ll_dphi", 20u, 0., 3.}, "ll_dphi", "weight");
    auto histPt = selected.Histo1D({"pt", "ll_pt", 30u, 0., 200.}, "ll_pt", "weight");
    auto histMET = selected.Histo1D({"met", "met_et", 20u, 0., 200.}, "met_et", "weight");
    auto histMass = selected.Histo1D({"mt", "mt", 15u, 50., 200.}, "mt", "weight");

    // Save Histograms
    TString output_name = "Output/" + sampleName + ".root";
    TFile physicsOutput(output_name, "recreate");
    physicsOutput.cd();
    histDPhi->Write();
    histPt->Write();
    histMET->Write();
    histMass->Write();
    physicsOutput.Close();

    std::cout << "Finish Analyzing " << sampleName << "." << std::endl;
}

RVecI GoodLepton(RVecB lep_isTightID, RVecF lep_ptcone30, RVecF lep_etcone20, RVecF lep_pt, RVecF lep_eta, RVecF lep_trackd0pvunbiased, RVecF lep_tracksigd0pvunbiased, RVecI lep_type, RVecF lep_z0, RVecF lep_phi, RVecF lep_E) {
    RVecI good_lep;
    int num = lep_isTightID.size();

    for (int i=0; i < num; i++) {
        if (!lep_isTightID.at(i)) {
            good_lep.push_back(0);
            continue;
        }

        if (lep_pt.at(i) <= 15000.) {
            good_lep.push_back(0);
            continue;
        }

        if (lep_ptcone30.at(i)/lep_pt.at(i) >= 0.15 || lep_etcone20.at(i)/lep_pt.at(i) >= 0.15) {
            good_lep.push_back(0);
            continue;
        }

        Math::PtEtaPhiEVector leptemp = Math::PtEtaPhiEVector(lep_pt.at(i)/1000., lep_eta.at(i), lep_phi.at(i), lep_E.at(i)/1000.);

        if (lep_type.at(i) == 11) {
            if (!(TMath::Abs(lep_eta.at(i)) < 2.47) && (TMath::Abs(lep_eta.at(i)) < 1.37 || TMath::Abs(lep_eta.at(i)) > 1.52)) {
                good_lep.push_back(0);
                continue;
            }

            if (!(TMath::Abs(lep_trackd0pvunbiased.at(i))/lep_tracksigd0pvunbiased.at(i) < 5. && TMath::Abs(lep_z0.at(i)*TMath::Sin(leptemp.Theta())) < 0.5)) {
                good_lep.push_back(0);
                continue;
            }
        }

        else if (lep_type.at(i) == 13) {
            if (!(TMath::Abs(lep_eta.at(i)) < 2.5)) {
                good_lep.push_back(0);
                continue;
            }
            
            if (!(TMath::Abs(lep_trackd0pvunbiased.at(i))/lep_tracksigd0pvunbiased.at(i) < 3. && TMath::Abs(lep_z0.at(i)*TMath::Sin(leptemp.Theta())) < 0.5)) {
                good_lep.push_back(0);
                continue;
            }
        }

        good_lep.push_back(1);
    }

    return good_lep;
}

float ComputeInvMass(RVecF lep_pt, RVecF lep_eta, RVecF lep_phi, RVecF lep_E) {
    Math::PtEtaPhiEVector lepton_1 = Math::PtEtaPhiEVector(lep_pt[0], lep_eta[0], lep_phi[0], lep_E[0]);
    Math::PtEtaPhiEVector lepton_2 = Math::PtEtaPhiEVector(lep_pt[1], lep_eta[1], lep_phi[1], lep_E[1]);

    auto lepton = lepton_1 + lepton_2;

    return lepton.mass()/1000.;
}

float ComputePhiLL(RVecF lep_pt, RVecF lep_eta, RVecF lep_phi, RVecF lep_E) {
    Math::PtEtaPhiEVector lepton_1 = Math::PtEtaPhiEVector(lep_pt[0], lep_eta[0], lep_phi[0], lep_E[0]);
    Math::PtEtaPhiEVector lepton_2 = Math::PtEtaPhiEVector(lep_pt[1], lep_eta[1], lep_phi[1], lep_E[1]);

    auto lepton = lepton_1 + lepton_2;

    return lepton.Phi();
}

float ComputeDPhiLL(RVecF lep_pt, RVecF lep_eta, RVecF lep_phi, RVecF lep_E) {
    Math::PtEtaPhiEVector lepton_1 = Math::PtEtaPhiEVector(lep_pt[0], lep_eta[0], lep_phi[0], lep_E[0]);
    Math::PtEtaPhiEVector lepton_2 = Math::PtEtaPhiEVector(lep_pt[1], lep_eta[1], lep_phi[1], lep_E[1]);

    auto lepton = lepton_1 + lepton_2;

    return TMath::Abs(lepton_1.Phi() - lepton_2.Phi());
}

float ComputePtLL(RVecF lep_pt, RVecF lep_eta, RVecF lep_phi, RVecF lep_E) {
    Math::PtEtaPhiEVector lepton_1 = Math::PtEtaPhiEVector(lep_pt[0], lep_eta[0], lep_phi[0], lep_E[0]);
    Math::PtEtaPhiEVector lepton_2 = Math::PtEtaPhiEVector(lep_pt[1], lep_eta[1], lep_phi[1], lep_E[1]);

    auto lepton = lepton_1 + lepton_2;

    return lepton.Pt()/1000.;
}

float ComputeMt(RVecF lep_pt, RVecF lep_eta, RVecF lep_phi, RVecF lep_E, float met_et, float met_phi) {
    Math::PtEtaPhiEVector lepton_1 = Math::PtEtaPhiEVector(lep_pt[0], lep_eta[0], lep_phi[0], lep_E[0]);
    Math::PtEtaPhiEVector lepton_2 = Math::PtEtaPhiEVector(lep_pt[1], lep_eta[1], lep_phi[1], lep_E[1]);

    auto lepton = lepton_1 + lepton_2;

    Math::PtEtaPhiEVector met = Math::PtEtaPhiEVector(met_et, 0, met_phi, met_et);

    return (lepton + met).Mt()/1000.;
}