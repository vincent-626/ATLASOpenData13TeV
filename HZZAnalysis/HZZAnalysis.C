#define HZZAnalysis_cxx
// The class definition in HZZAnalysis.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("HZZAnalysis.C")
// root> T->Process("HZZAnalysis.C","some options")
// root> T->Process("HZZAnalysis.C+")
//


#include "HZZAnalysis.h"
#include <TH1.h>
#include <TStyle.h>
#include <TLorentzVector.h>

void HZZAnalysis::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void HZZAnalysis::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   printf("Starting analysis with process option: %s \n", option.Data());

   nEvent = fChain->GetTree()->GetEntries();
   iEvent = 0;

   // 4l mass histogram
   h1_fourlep = new TH1F("h1_fourlep", "Mass of four-lepton system; m_{4l} [GeV];Events / bin", 24, 80, 170);
}

Bool_t HZZAnalysis::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fChain->GetTree()->GetEntry(entry);

   iEvent++;
   if (iEvent % 50000 == 0) std::cout << "Analyzing " << iEvent << "out of " << nEvent << std::endl;

   // Scale factors
   Float_t scaleFactor = scaleFactor_ELE * scaleFactor_MUON * scaleFactor_LepTRIGGER * scaleFactor_PILEUP;
   
   // MC weight
   Float_t m_mcWeight = mcWeight;

   // Read input option
   TString option = GetOption();
   if (option.Contains("single")) m_mcWeight = (mcWeight/TMath::Abs(mcWeight)); // set to 1 or -1 for single top MCs

   // Total weight
   Float_t weight = scaleFactor * m_mcWeight;

   // For data
   if (option.Contains("data")) weight = 1.;

   // e/mu trigger
   if (!(trigE || trigM)) return kTRUE;

   int goodlep_index[lep_n];
   int goodlep_n = 0;
   int lep_index = 0;

   // Lepton selection
   for (unsigned int i; i < lep_n; i++) {
      TLorentzVector leptemp;
      leptemp.SetPtEtaPhiE(lep_pt->at(i)/1000., lep_eta->at(i), lep_phi->at(i), lep_E->at(i)/1000.);

      // loosely isolated and very soft
      if (!(lep_pt->at(i) >5000. && TMath::Abs(lep_eta->at(i) < 2.5) && (lep_ptcone30->at(i)/lep_pt->at(i)) < 0.3 && (lep_etcone20->at(i) / lep_pt->at(i)) < 0.3)) continue;

      // electron
      if (lep_type->at(i) == 11) {
         if (!(lep_pt->at(i) > 7000. && TMath::Abs(lep_eta->at(i)) < 2.47)) continue;
         if (!(TMath::Abs(lep_trackd0pvunbiased->at(i))/lep_tracksigd0pvunbiased->at(i) < 5 && TMath::Abs(lep_z0->at(i)*TMath::Sin(leptemp.Theta())) < 0.5)) continue;
      }
      
      // muon
      else if (lep_type->at(i) == 13) {
         if (!(TMath::Abs(lep_trackd0pvunbiased->at(i))/lep_tracksigd0pvunbiased->at(i) < 3 && TMath::Abs(lep_z0->at(i)*TMath::Sin(leptemp.Theta())) < 0.5)) continue;
      }

      goodlep_n++;
      goodlep_index[lep_index] = i;
      lep_index++;
   }

   // Exactly 4 leptons
   if (goodlep_n != 4) return kTRUE;

   int goodlep1_index = goodlep_index[0];
   int goodlep2_index = goodlep_index[1];
   int goodlep3_index = goodlep_index[2];
   int goodlep4_index = goodlep_index[3];

   // First lepton pT > 25 GeV, second > 15 GeV and third > 10 GeV
   if (!(lep_pt->at(goodlep1_index) > 25000. && lep_pt->at(goodlep2_index) > 15000. && lep_pt->at(goodlep3_index) > 10000.)) return kTRUE;

   TLorentzVector Lepton_1  = TLorentzVector();
   TLorentzVector Lepton_2  = TLorentzVector();
   TLorentzVector Lepton_3  = TLorentzVector();
   TLorentzVector Lepton_4  = TLorentzVector();

   Lepton_1.SetPtEtaPhiE(lep_pt->at(goodlep1_index), lep_eta->at(goodlep1_index), lep_phi->at(goodlep1_index), lep_E->at(goodlep1_index));
   Lepton_2.SetPtEtaPhiE(lep_pt->at(goodlep2_index), lep_eta->at(goodlep2_index), lep_phi->at(goodlep2_index), lep_E->at(goodlep2_index));
   Lepton_3.SetPtEtaPhiE(lep_pt->at(goodlep3_index), lep_eta->at(goodlep3_index), lep_phi->at(goodlep3_index), lep_E->at(goodlep3_index));
   Lepton_4.SetPtEtaPhiE(lep_pt->at(goodlep4_index), lep_eta->at(goodlep4_index), lep_phi->at(goodlep4_index), lep_E->at(goodlep4_index));
		  
   TLorentzVector FourLepSystem = TLorentzVector();
   FourLepSystem = Lepton_1 + Lepton_2 + Lepton_3 + Lepton_4;
   float FourLepSystem_M = FourLepSystem.Mag()/1000.;

   h1_fourlep->Fill(FourLepSystem_M, weight);

   return kTRUE;
}

void HZZAnalysis::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void HZZAnalysis::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   TString filename_option = GetOption();
   printf("Writting with name option: %s \n", filename_option.Data());

   TFile output("Output/"+filename_option+".root", "recreate");
   h1_fourlep->Write();
   output.Close();
}