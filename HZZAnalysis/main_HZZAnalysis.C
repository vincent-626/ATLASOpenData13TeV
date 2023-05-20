#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TProof.h"

void main_HZZAnalysis() {
    TString path_4lep = "/Users/vincent626/Code/CERN/opendata/ATLASOpenData13TeV/Samples/4lep/";
    TString path_2lep = "/Users/vincent626/Code/CERN/opendata/ATLASOpenData13TeV/Samples/2lep/";

    // Parallel processing
    TProof::Open("");

    // ------
    // Data
    // ------

    TChain* chain_data = new TChain("mini");

    chain_data->AddFile(path_4lep + "Data/data_A.4lep.root");
    chain_data->AddFile(path_4lep + "Data/data_B.4lep.root");
    chain_data->AddFile(path_4lep + "Data/data_C.4lep.root");
    chain_data->AddFile(path_4lep + "Data/data_D.4lep.root");

    chain_data->SetProof();
    chain_data->Process("HZZAnalysis.C+", "data");

    // // ----------
    // // Signal MC
    // // ----------

    // // Higgs
    // TChain* chain_higgs = new TChain("mini");

    // chain_higgs->AddFile(path_4lep + "MC/mc_345060.ggH125_ZZ4lep.4lep.root");
    // chain_higgs->AddFile(path_4lep + "MC/mc_341947.ZH125_ZZ4lep.4lep.root");
    // chain_higgs->AddFile(path_4lep + "MC/mc_341964.WH125_ZZ4lep.4lep.root");
    // chain_higgs->AddFile(path_4lep + "MC/mc_344235.VBFH125_ZZ4lep.4lep.root");

    // chain_higgs->SetProof();
    // chain_higgs->Process("HZZAnalysis.C+", "higgs");

    // // --------------
    // // Background MC
    // // --------------

    // // TT bar
    // TChain* chain_ttbar = new TChain("mini");

    // chain_ttbar->AddFile(path_4lep + "MC/mc_410000.ttbar_lep.4lep.root");

    // chain_ttbar->SetProof();
    // chain_ttbar->Process("HZZAnalysis.C+", "ttbar");

    // // Single top
    // TChain* chain_singlet = new TChain("mini");

    // chain_singlet->AddFile(path_4lep + "MC/mc_410011.single_top_tchan.4lep.root");
    // chain_singlet->AddFile(path_4lep + "MC/mc_410012.single_antitop_tchan.4lep.root");
    // chain_singlet->AddFile(path_4lep + "MC/mc_410025.single_top_schan.4lep.root");
    // chain_singlet->AddFile(path_4lep + "MC/mc_410026.single_antitop_schan.4lep.root");
    // chain_singlet->AddFile(path_4lep + "MC/mc_410013.single_top_wtchan.4lep.root");
    // chain_singlet->AddFile(path_4lep + "MC/mc_410014.single_antitop_wtchan.4lep.root");

    // chain_singlet->SetProof();
    // chain_singlet->Process("HZZAnalysis.C+", "singleTop");

    // // W+jets
    // TChain* chain_wjets = new TChain("mini");

    // chain_wjets->AddFile(path_2lep + "MC/mc_361100.Wplusenu.2lep.root");
    // chain_wjets->AddFile(path_2lep + "MC/mc_361101.Wplusmunu.2lep.root");
    // chain_wjets->AddFile(path_2lep + "MC/mc_361102.Wplustaunu.2lep.root");
    // chain_wjets->AddFile(path_2lep + "MC/mc_361103.Wminusenu.2lep.root");
    // chain_wjets->AddFile(path_2lep + "MC/mc_361104.Wminusmunu.2lep.root");
    // chain_wjets->AddFile(path_2lep + "MC/mc_361105.Wminustaunu.2lep.root");

    // chain_wjets->SetProof();
    // chain_wjets->Process("HZZAnalysis.C+", "w+jets");

    // // Z+jets
    // TChain* chain_zjets = new TChain("mini");

    // chain_zjets->AddFile(path_4lep + "MC/mc_361106.Zee.4lep.root");
    // chain_zjets->AddFile(path_4lep + "MC/mc_361107.Zmumu.4lep.root");
    // chain_zjets->AddFile(path_4lep + "MC/mc_361108.Ztautau.4lep.root");

    // chain_zjets->SetProof();
    // chain_zjets->Process("HZZAnalysis.C+", "z+jets");
    
    // // Diboson
    // TChain* chain_diboson = new TChain("mini");

    // chain_diboson->AddFile(path_4lep + "MC/mc_363356.ZqqZll.4lep.root");
    // chain_diboson->AddFile(path_2lep + "MC/mc_363358.WqqZll.2lep.root");
    // chain_diboson->AddFile(path_2lep + "MC/mc_363359.WpqqWmlv.2lep.root");
    // chain_diboson->AddFile(path_2lep + "MC/mc_363360.WplvWmqq.2lep.root");
    // chain_diboson->AddFile(path_2lep + "MC/mc_363489.WlvZqq.2lep.root");
    // chain_diboson->AddFile(path_4lep + "MC/mc_363490.llll.4lep.root");
    // chain_diboson->AddFile(path_4lep + "MC/mc_363491.lllv.4lep.root");
    // chain_diboson->AddFile(path_4lep + "MC/mc_363492.llvv.4lep.root");
    // chain_diboson->AddFile(path_2lep + "MC/mc_363493.lvvv.2lep.root");

    // chain_diboson->SetProof();
    // chain_diboson->Process("HZZAnalysis.C+", "diboson");
}