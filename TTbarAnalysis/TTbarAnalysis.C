// TTbar Analysis using RDataFrame
// Run the code with root -b -l -q TTbarAnalysis.C

#include <iostream>
#include <TROOT.h>
using namespace ROOT;

void Analyze(TString, TString);
RVecU GoodLep(unsigned int, RVecB, RVecF, RVecF, RVecF, unsigned int, RVecF, RVecF, RVecF);
float ComputeMTW(RVecF, float, RVecF, float);
RVecU GoodJet(unsigned int, RVecF, RVecF, RVecF);
float ComputeSystPt(unsigned int, RVecF, RVecF);
float ComputeSystMet(float, float, float);
float ComputeSystMTW();


void TTbarAnalysis() {
    TString inDir = "/home/vincent626/opendata/atlas2020/1lep/";

    // Data Samples
    Analyze(inDir + "Data/data_A.1lep.root", "data_A.1lep");
    Analyze(inDir + "Data/data_B.1lep.root", "data_B.1lep");
    Analyze(inDir + "Data/data_C.1lep.root", "data_C.1lep");
    Analyze(inDir + "Data/data_D.1lep.root", "data_D.1lep");

    // MC Samples
    Analyze(inDir + "Data/data_A.1lep.root", "data_A.1lep");
    Analyze(inDir + "Data/data_B.1lep.root", "data_B.1lep");
    Analyze(inDir + "Data/data_C.1lep.root", "data_C.1lep");
    Analyze(inDir + "Data/data_D.1lep.root", "data_D.1lep");
}

void Analyze(TString inFile, TString sampleName) {
    // Parallel Processing
    EnableImplicitMT();

    // Printing Status
    std::cout << "Analyzing " << sampleName << "..." << std::endl;

    // Construct RDataFrame
    TFile *f = TFile::Open(inFile);
    RDataFrame df("mini", f);

    // Define Weight
    auto weighted = df.Define("weight", "1.0");

    if (!sampleName.Contains("data")) {
        weighted = df.Define("weight", "scaleFactor_ELE * scaleFactor_MUON * scaleFactor_LepTRIGGER * scaleFactor_PILEUP * scaleFactor_BTAG * mcWeight");
    }

    ////////////////////////////////////////////////
    // Event Selection
    ////////////////////////////////////////////////

    // e/mu Trigger
    auto trigger = weighted.Filter("trigE || trigM");

    auto cTrigger = trigger.Count();
    std::cout << "Trigger: " << *cTrigger << std::endl;

    // Define Good Lepton (Pt > 30GeV & Tight ID & Isolation)
    // Exactly 1 Good Lepton
    auto goodLep = trigger.Define("good_lep", "GoodLep(lep_n, lep_isTightID, lep_pt, lep_ptcone30, lep_etcone20, lep_type, lep_eta, lep_trackd0pvunbiased, lep_tracksigd0pvunbiased)")
                            .Filter("Sum(good_lep) == 1");

    auto cLep = goodLep.Count();
    std::cout << "Good Lepton: " << *cLep << std::endl;

    // MET > 30GeV
    auto metCut = goodLep.Filter("met_et > 30000.");

    auto cMet = metCut.Count();
    std::cout << "MET: " << *cMet << std::endl;

    // Define Good Jet (Pt > 30GeV & JVT Cleaning)
    // >= 4 Jets
    auto goodJet = metCut.Define("good_jet", "GoodJet(jet_n, jet_pt, jet_eta, jet_jvt)")
                            .Filter("Sum(good_jet) >= 4");

    auto cJet = goodJet.Count();
    std::cout << "Jet: " << *cJet << std::endl;

    // Define MTW
    // MTW > 30GeV
    auto mtwCut = goodJet.Define("mtw", "ComputeMTW(lep_pt, met_et, lep_phi, met_phi)")
                            .Filter("mtw > 30000.");

    auto cMtw = mtwCut.Count();
    std::cout << "MTW: " << *cMtw << std::endl;

    // Define Good B-Jet (70% WP)
    // >= 2 Good B-jets
    auto bJet = mtwCut.Define("good_bjet", "good_jet && jet_MV2c10 > 0.8244273")
                        .Filter("Sum(good_bjet) >= 2");
    
    auto cBJet = bJet.Count();
    std::cout << "B-jet: " << *cBJet << std::endl;

    // Define Leading Jet Pt
    auto leadPt = bJet.Define("leadJet_pt", "jet_pt[good_jet][0]");

    // Systematic Variation and Cut
    auto syst = leadPt.Define("systPt", "ComputeSystPt(jet_n, jet_pt[good_jet], jet_pt_syst[good_jet])")
                        .Define("systMet", "ComputeSystMet(met_et, met_phi, met_et_syst)")
                        .Define("systMTW", "ComputeSystMTW(systPt, systMet, lep_phi[good_lep], met_phi)")
                        .Filter("systPt > 30000.")
                        .Filter("systMet > 30000.")
                        .Filter("systMTW > 30000.");

    auto cVar = syst.Count();
    std::cout << "Systematic Variation: " << *cVar << std::endl;


    ////////////////////////////////////////////////
    // Histogramming
    ////////////////////////////////////////////////

    // Fill Histograms
    auto histJetPt = leadPt.Histo1D("leadJet_pt", "weight");
    auto histSystPt = syst.Histo1D("systPt", "weight");
    auto histMet = leadPt.Histo1D("met_et", "weight");
    auto histSystMet = syst.Histo1D("systMet", "weight");
    auto histMtw = leadPt.Histo1D("mtw", "weight");
    auto histSystMTW = syst.Histo1D("systMTW", "weight");
    auto histNjet = leadPt.Histo1D("jet_n", "weight");

    // Save Histograms
    TString output_name = "Output/" + sampleName + ".root";
    TFile physicsOutput(output_name, "recreate");
    physicsOutput.cd();
    histJetPt->Write();
    histSystPt->Write();
    histMet->Write();
    histSystMet->Write();
    histMtw->Write();
    histSystMTW->Write();
    histNjet->Write();
    physicsOutput.Close();

    std::cout << "Finish Analyzing " << sampleName << "." << std::endl << std::endl;
}

RVecU GoodLep(unsigned int lep_n, RVecB lep_isTightID, RVecF lep_pt, RVecF lep_ptcone30, RVecF lep_etcone20, RVecU lep_type, RVecF lep_eta, RVecF lep_trackd0pvunbiased, RVecF lep_tracksigd0pvunbiased) {
    RVecU good_lep;

    if (!lep_n) return good_lep;

    for (int i=0; i<lep_n; i++) {
        if (!(lep_isTightID.at(i) && lep_pt.at(i) > 30000. && lep_ptcone30.at(i)/lep_pt.at(i) < 0.15 && lep_etcone20.at(i)/lep_pt.at(i) < 0.15)) {
            good_lep.push_back(0);
            continue;
        }

        if (lep_type.at(i) == 11) {
            if (!(TMath::Abs(lep_eta.at(i) < 2.47 && (TMath::Abs(lep_eta.at(i)) < 1.37 || TMath::Abs(lep_eta.at(i)) > 1.52)))) {
                good_lep.push_back(0);
                continue;
            }

            if (!(TMath::Abs(lep_trackd0pvunbiased.at(i)/lep_tracksigd0pvunbiased.at(i)) < 5. && TMath::Abs(lep_trackd0pvunbiased.at(i)/lep_tracksigd0pvunbiased.at(i)) < 0.5)) {
                good_lep.push_back(0);
                continue;
            }
        }

        else if (lep_type.at(i) == 13) {
            if (!(TMath::Abs(lep_eta.at(i)) < 2.5)) {
                good_lep.push_back(0);
                continue;
            }

            if (!(TMath::Abs(lep_trackd0pvunbiased.at(i)/lep_tracksigd0pvunbiased.at(i)) < 3. && TMath::Abs(lep_trackd0pvunbiased.at(i)/lep_tracksigd0pvunbiased.at(i)) < 0.5)) {
                good_lep.push_back(0);
                continue;
            }
        }

        good_lep.push_back(1);
    }

    return good_lep;
}

float ComputeMTW(RVecF lep_pt, float met_et, RVecF lep_phi, float met_phi) {
    float mtw;

    mtw = TMath::Sqrt(2 * lep_pt[0] * met_et * (1 - TMath::Cos(TMath::Abs(lep_phi[0] - met_phi))));

    return mtw;
}

RVecU GoodJet(unsigned int jet_n, RVecF jet_pt, RVecF jet_eta, RVecF jet_jvt) {
    RVecU good_jet;

    if (!jet_n) return good_jet;

    for (int i=0; i<jet_n; i++) {
        if (jet_pt.at(i) <= 30000.) {
            good_jet.push_back(0);
            continue;
        }
        
        if (jet_pt.at(i) < 60000. && TMath::Abs(jet_eta.at(i)) < 2.4 && jet_jvt.at(i) < 0.59) {
            good_jet.push_back(0);
            continue;
        }

        good_jet.push_back(1);
    }
    
    return good_jet;
}

float ComputeSystPt(unsigned int n, RVecF pt, RVecF pt_syst) {
    float systPt;

    TRandom3 *gRand = new TRandom3(0);

    systPt = gRand->Gaus(pt[0], pt_syst[0]);

    return systPt;
}

float ComputeSystMet(float met_et, float met_phi, float met_et_syst) {
    float systMet;

    TRandom3* gRand = new TRandom3(0);

	systMet = gRand->Gaus(met_et, met_et_syst);

    return systMet;
}

float ComputeSystMTW(float systPt, float systMet, RVecF lep_phi, float met_phi) {
    float SystMTW;

    SystMTW = TMath::Sqrt(2 * systPt * systMet * (1 - TMath::Cos(TMath::Abs(lep_phi[0] - met_phi))));

    return SystMTW;
}