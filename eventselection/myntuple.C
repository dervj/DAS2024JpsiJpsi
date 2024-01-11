#define myntuple_cxx
#include "myntuple.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TLorentzVector.h"
#include <iostream> 
#include <fstream>
#include <numeric>
using namespace std;

struct PairedMuonIdx {
	int p11;
	int p12;
	int p21;
	int p22;
};

void myntuple::Loop()
{
	//   In a ROOT session, you can do:
	//      root> .L myntuple.C
	//      root> myntuple t
	//      root> t.GetEntry(12); // Fill t data members with entry number 12
	//      root> t.Show();       // Show values of entry 12
	//      root> t.Show(16);     // Read and show values of entry 16
	//      root> t.Loop();       // Loop on all entries
	//

	//     This is the loop skeleton where:
	//    jentry is the global entry number in the chain
	//    ientry is the entry number in the current Tree
	//  Note that the argument to GetEntry must be:
	//    jentry for TChain::GetEntry
	//    ientry for TTree::GetEntry and TBranch::GetEntry
	//
	//       To read only selected branches, Insert statements like:
	// METHOD1:
	//    fChain->SetBranchStatus("*",0);  // disable all branches
	//    fChain->SetBranchStatus("branchname",1);  // activate branchname
	// METHOD2: replace line
	//    fChain->GetEntry(jentry);       //read all branches
	//by  b_branchname->GetEntry(ientry); //read only this branch

	//define txt file
	ofstream myoutputfile("mJJDataFull6000_15000.txt");

	//define root file
	TFile* myhbk = new TFile ("myhbk.root","recreate");

	//define histogram
	TH1F* myDiMuon1mass = new TH1F("myDiMuon1mass","myDiMuon1mass",80,2.7,3.5);
	TH1F* myDiMuon2mass = new TH1F("myDiMuon2mass","myDiMuon2mass",80,2.7,3.5);
	TH1F* myFourMuonmass = new TH1F("myFourMuonmass","myFourMuonmass",360,6,15);

	TH1F* myDiMuon1mass_PT = new TH1F("myDiMuon1mass_PT","myDiMuon1mass",80,2.7,3.5);
	TH1F* myDiMuon2mass_PT = new TH1F("myDiMuon2mass_PT","myDiMuon2mass",80,2.7,3.5);
	TH1F* myFourMuonmass_PT = new TH1F("myFourMuonmass_PT","myFourMuonmass",360,6,15);

	TH1F* myDiMuon1mass_eta = new TH1F("myDiMuon1mass_eta","myDiMuon1mass",80,2.7,3.5);
	TH1F* myDiMuon2mass_eta = new TH1F("myDiMuon2mass_eta","myDiMuon2mass",80,2.7,3.5);
	TH1F* myFourMuonmass_eta = new TH1F("myFourMuonmass_eta","myFourMuonmass",360,6,15);

	TH1F* myDiMuon1mass_chg = new TH1F("myDiMuon1mass_chg","myDiMuon1mass",80,2.7,3.5);
	TH1F* myDiMuon2mass_chg = new TH1F("myDiMuon2mass_chg","myDiMuon2mass",80,2.7,3.5);
	TH1F* myFourMuonmass_chg = new TH1F("myFourMuonmass_chg","myFourMuonmass",360,6,15);

	TH1F* myDiMuon1mass_mass = new TH1F("myDiMuon1mass_mass","myDiMuon1mass",80,2.7,3.5);
	TH1F* myDiMuon2mass_mass = new TH1F("myDiMuon2mass_mass","myDiMuon2mass",80,2.7,3.5);
	TH1F* myFourMuonmass_mass = new TH1F("myFourMuonmass_mass","myFourMuonmass",360,6,15);

	const double MUON_MASS = 0.1056583745; //GeV
	const double JPSI_MASS = 3.096900; //GeV

	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntries();

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		if (jentry%10000 == 0) cout << "I am running " << jentry << "th entries out of " << nentries << " total entries" << endl;
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;

		bool TrigThreeMuonJpsi = false;
		bool TrigThreeMuonJpsi3p5mu2 = false;

		for (unsigned int i = 0; i != TrigRes->size(); ++i) { 
			if ((TrigNames->at(i).find("HLT_Dimuon0_Jpsi3p5_Muon2_") != string::npos && TrigRes->at(i) == 1)) {
				TrigThreeMuonJpsi3p5mu2 = true;
			}
			if (TrigNames->at(i).find("HLT_Dimuon0_Jpsi_Muon_") != string::npos && TrigRes->at(i) == 1) {
				TrigThreeMuonJpsi = true;
			}
		} //end of trigger loop
		if (!TrigThreeMuonJpsi3p5mu2 && !TrigThreeMuonJpsi){
			continue;
		}
		for (unsigned int myFourMuIdx = 0; myFourMuIdx < nMyFourMuon; myFourMuIdx++) {
			TLorentzVector myFourMuonP4;
			myFourMuonP4.SetXYZM((*MyFourMuonPx)[myFourMuIdx],(*MyFourMuonPy)[myFourMuIdx], (*MyFourMuonPz)[myFourMuIdx], (*MyFourMuonMass)[myFourMuIdx]);

			vector<int> fitMuCharge;
			fitMuCharge.push_back((*muCharge)[(*MyFourMuonMu1Idx)[myFourMuIdx]]);
			fitMuCharge.push_back((*muCharge)[(*MyFourMuonMu2Idx)[myFourMuIdx]]);
			fitMuCharge.push_back((*muCharge)[(*MyFourMuonMu3Idx)[myFourMuIdx]]);
			fitMuCharge.push_back((*muCharge)[(*MyFourMuonMu4Idx)[myFourMuIdx]]);

			vector<TLorentzVector> rawMup4vect, rawMuinFourMuFMp4vect;
			TLorentzVector Rawmu, RawmuinFourMuFM;
			float raw_muPx = (*muPx)[(*MyFourMuonMu1Idx)[myFourMuIdx]];
			float raw_muPy = (*muPy)[(*MyFourMuonMu1Idx)[myFourMuIdx]];
			float raw_muPz = (*muPz)[(*MyFourMuonMu1Idx)[myFourMuIdx]];
			Rawmu.SetXYZM(raw_muPx, raw_muPy,raw_muPz, MUON_MASS);   
			RawmuinFourMuFM = Rawmu; RawmuinFourMuFM.Boost(-myFourMuonP4.BoostVector());
			rawMup4vect.push_back(Rawmu); rawMuinFourMuFMp4vect.push_back(RawmuinFourMuFM);
			raw_muPx = (*muPx)[(*MyFourMuonMu2Idx)[myFourMuIdx]]; raw_muPy = (*muPy)[(*MyFourMuonMu2Idx)[myFourMuIdx]];raw_muPz = (*muPz)[(*MyFourMuonMu2Idx)[myFourMuIdx]];
			Rawmu.SetXYZM(raw_muPx, raw_muPy,raw_muPz, MUON_MASS);  rawMup4vect.push_back(Rawmu);
			RawmuinFourMuFM = Rawmu; RawmuinFourMuFM.Boost(-myFourMuonP4.BoostVector());
			rawMuinFourMuFMp4vect.push_back(RawmuinFourMuFM);
			raw_muPx = (*muPx)[(*MyFourMuonMu3Idx)[myFourMuIdx]]; raw_muPy = (*muPy)[(*MyFourMuonMu3Idx)[myFourMuIdx]];raw_muPz = (*muPz)[(*MyFourMuonMu3Idx)[myFourMuIdx]];
			Rawmu.SetXYZM(raw_muPx, raw_muPy,raw_muPz, MUON_MASS);  rawMup4vect.push_back(Rawmu);
			RawmuinFourMuFM = Rawmu; RawmuinFourMuFM.Boost(-myFourMuonP4.BoostVector());
			rawMuinFourMuFMp4vect.push_back(RawmuinFourMuFM);
			raw_muPx = (*muPx)[(*MyFourMuonMu4Idx)[myFourMuIdx]]; raw_muPy = (*muPy)[(*MyFourMuonMu4Idx)[myFourMuIdx]];raw_muPz = (*muPz)[(*MyFourMuonMu4Idx)[myFourMuIdx]];
			Rawmu.SetXYZM(raw_muPx, raw_muPy,raw_muPz, MUON_MASS);  rawMup4vect.push_back(Rawmu);
			RawmuinFourMuFM = Rawmu; RawmuinFourMuFM.Boost(-myFourMuonP4.BoostVector());
			rawMuinFourMuFMp4vect.push_back(RawmuinFourMuFM);

			// int muonPtPass = 0;
			// for(unsigned int )

			vector<TLorentzVector> fitMup4vect;
			TLorentzVector Fitmu;
			//Muon From X Fit:        
			//Fit Muon 1
			float fit_muPx = (*MyFourMuonMu1Px)[myFourMuIdx];
			float fit_muPy = (*MyFourMuonMu1Py)[myFourMuIdx];
			float fit_muPz = (*MyFourMuonMu1Pz)[myFourMuIdx];
			float fit_muM   = MUON_MASS;
			Fitmu.SetXYZM(fit_muPx, fit_muPy,fit_muPz, fit_muM);
			fitMup4vect.push_back(Fitmu);

			fit_muPx = (*MyFourMuonMu2Px)[myFourMuIdx]; fit_muPy = (*MyFourMuonMu2Py)[myFourMuIdx]; fit_muPz = (*MyFourMuonMu2Pz)[myFourMuIdx];
			Fitmu.SetXYZM(fit_muPx, fit_muPy,fit_muPz, fit_muM);  fitMup4vect.push_back(Fitmu);
			fit_muPx = (*MyFourMuonMu3Px)[myFourMuIdx]; fit_muPy = (*MyFourMuonMu3Py)[myFourMuIdx]; fit_muPz = (*MyFourMuonMu3Pz)[myFourMuIdx];
			Fitmu.SetXYZM(fit_muPx, fit_muPy,fit_muPz, fit_muM);  fitMup4vect.push_back(Fitmu);
			fit_muPx = (*MyFourMuonMu4Px)[myFourMuIdx]; fit_muPy = (*MyFourMuonMu4Py)[myFourMuIdx]; fit_muPz = (*MyFourMuonMu4Pz)[myFourMuIdx];
			Fitmu.SetXYZM(fit_muPx, fit_muPy,fit_muPz, fit_muM);  fitMup4vect.push_back(Fitmu);

			PairedMuonIdx myCombIdx[3];   
			myCombIdx[0].p11 = 0; myCombIdx[0].p12 = 1; myCombIdx[0].p21 = 2; myCombIdx[0].p22 = 3; 
			myCombIdx[1].p11 = 0; myCombIdx[1].p12 = 2; myCombIdx[1].p21 = 1; myCombIdx[1].p22 = 3; 
			myCombIdx[2].p11 = 0; myCombIdx[2].p12 = 3; myCombIdx[2].p21 = 1; myCombIdx[2].p22 = 2; 

			int myNumPatSoftMuon = (*muIsPatSoftMuon)[(*MyFourMuonMu1Idx)[myFourMuIdx]] + (*muIsPatSoftMuon)[(*MyFourMuonMu2Idx)[myFourMuIdx]] + (*muIsPatSoftMuon)[(*MyFourMuonMu3Idx)[myFourMuIdx]] + (*muIsPatSoftMuon)[(*MyFourMuonMu4Idx)[myFourMuIdx]];
			int myNumPatLooseMuon = (*muIsPatLooseMuon)[(*MyFourMuonMu1Idx)[myFourMuIdx]] + (*muIsPatLooseMuon)[(*MyFourMuonMu2Idx)[myFourMuIdx]] + (*muIsPatLooseMuon)[(*MyFourMuonMu3Idx)[myFourMuIdx]] + (*muIsPatLooseMuon)[(*MyFourMuonMu4Idx)[myFourMuIdx]];
			int myNumPatMediumMuon = (*muIsPatMediumMuon)[(*MyFourMuonMu1Idx)[myFourMuIdx]] + (*muIsPatMediumMuon)[(*MyFourMuonMu2Idx)[myFourMuIdx]] + (*muIsPatMediumMuon)[(*MyFourMuonMu3Idx)[myFourMuIdx]] + (*muIsPatMediumMuon)[(*MyFourMuonMu4Idx)[myFourMuIdx]];
			int myNumPatTightMuon = (*muIsPatTightMuon)[(*MyFourMuonMu1Idx)[myFourMuIdx]] + (*muIsPatTightMuon)[(*MyFourMuonMu2Idx)[myFourMuIdx]] + (*muIsPatTightMuon)[(*MyFourMuonMu3Idx)[myFourMuIdx]] + (*muIsPatTightMuon)[(*MyFourMuonMu4Idx)[myFourMuIdx]];
			float DiMuonMass1 = 0.; 
			float DiMuonMass2 = 0.;
                        double m4Muon = 0.;



			if (1
					// soft muon: tracker muon + 1 hit in the muon system 
					&& myNumPatSoftMuon >= 4 
					// requiring all muons to come from the same vertex
					&& (*MyFourMuonVtxCL)[myFourMuIdx] >= 0.005     
					// Here add the selections!
					&& (TrigThreeMuonJpsi3p5mu2 || TrigThreeMuonJpsi) == 1
					// &&  == 1

					// && accumulate(fitMuCharge.begin(), fitMuCharge.end(), 0) == 0
					&& (fitMuCharge[0] + fitMuCharge[1] + fitMuCharge[2] + fitMuCharge[3] == 0)
			) {
				for (int mypidx = 0; mypidx < 3; mypidx++)  {
					int muIdxp11, muIdxp12, muIdxp21, muIdxp22;
					muIdxp11 = myCombIdx[mypidx].p11; muIdxp12 = myCombIdx[mypidx].p12; muIdxp21 = myCombIdx[mypidx].p21; muIdxp22 = myCombIdx[mypidx].p22;

					if(1
                        // Here, require the muon pairs to have muons with opposite charges
						 && fitMuCharge[muIdxp11] + fitMuCharge[muIdxp12] == 0 
						 && fitMuCharge[muIdxp21] + fitMuCharge[muIdxp22] == 0
					  ) {
						DiMuonMass1 = (fitMup4vect[muIdxp11] + fitMup4vect[muIdxp12]).M();
						DiMuonMass2 = (fitMup4vect[muIdxp21] + fitMup4vect[muIdxp22]).M();
						// calculate the 4 muon mass:  M(µ1µ2µ3µ4)-M(µ1µ2)-M(µ3µ4)+2*M(J/psi)  
						m4Muon 	= (
							(fitMup4vect[muIdxp11] + fitMup4vect[muIdxp12] + fitMup4vect[muIdxp21] + fitMup4vect[muIdxp22]).M()
							- (fitMup4vect[muIdxp11] + fitMup4vect[muIdxp12]).M() - (fitMup4vect[muIdxp21] + fitMup4vect[muIdxp22]).M()
							+ 2*JPSI_MASS 
							);
						myDiMuon1mass_chg->Fill(DiMuonMass1);
						myDiMuon2mass_chg->Fill(DiMuonMass2);
						myFourMuonmass_chg->Fill(m4Muon);

						if (1
							&& rawMup4vect[0].Pt() >= 2.
							&& rawMup4vect[1].Pt() >= 2.
							&& rawMup4vect[2].Pt() >= 2.
							&& rawMup4vect[3].Pt() >= 2.
						){
							myDiMuon1mass_PT->Fill(DiMuonMass1);
							myDiMuon2mass_PT->Fill(DiMuonMass2);
							myFourMuonmass_PT->Fill(m4Muon);

							if (1
								&& abs(rawMup4vect[0].Eta()) <= 2.4
								&& abs(rawMup4vect[1].Eta()) <= 2.4
								&& abs(rawMup4vect[2].Eta()) <= 2.4
								&& abs(rawMup4vect[3].Eta()) <= 2.4
							){
								myDiMuon1mass_eta->Fill(DiMuonMass1);
								myDiMuon2mass_eta->Fill(DiMuonMass2);
								myFourMuonmass_eta->Fill(m4Muon);

								// Here require that each DiMuonMass is in the appropriate mass range [2.95,3.25] GeV
								if (1
									&& DiMuonMass1 <= 3.25
									&& DiMuonMass1 >= 2.95
									&& DiMuonMass2 <= 3.25
									&& DiMuonMass2 >= 2.95
								) {
									myDiMuon1mass_mass->Fill(DiMuonMass1);
									myDiMuon2mass_mass->Fill(DiMuonMass2);
									myFourMuonmass_mass->Fill(m4Muon);
								}
							}
						}
					}
				}
			}
		}}
			// if (1
			// 		// soft muon: tracker muon + 1 hit in the muon system 
			// 		&& myNumPatSoftMuon >= 4 
			// 		// requiring all muons to come from the same vertex
			// 		&& (*MyFourMuonVtxCL)[myFourMuIdx] >= 0.005     
			// 		// Here add the selections!
			// 		&& TrigThreeMuonJpsi3p5mu2 == 1
			// 		&& TrigThreeMuonJpsi == 1

			// 		//pt cut
			// 		&& rawMup4vect[0].Pt() >= 2.
			// 		&& rawMup4vect[1].Pt() >= 2.
			// 		&& rawMup4vect[2].Pt() >= 2.
			// 		&& rawMup4vect[3].Pt() >= 2.

			// 		&& abs(rawMup4vect[0].Eta()) <= 2.4
			// 		&& abs(rawMup4vect[1].Eta()) <= 2.4
			// 		&& abs(rawMup4vect[2].Eta()) <= 2.4
			// 		&& abs(rawMup4vect[3].Eta()) <= 2.4

					
			// 	) {
			// 	for (int mypidx = 0; mypidx < 3; mypidx++)  {
			// 		int muIdxp11, muIdxp12, muIdxp21, muIdxp22;
			// 		muIdxp11 = myCombIdx[mypidx].p11; muIdxp12 = myCombIdx[mypidx].p12; muIdxp21 = myCombIdx[mypidx].p21; muIdxp22 = myCombIdx[mypidx].p22;


			// 		if(1
            //                                     // Here, require the muon pairs to have muons with opposite charges
			// 			 && fitMuCharge(muIdxp11) + fitMuCharge(muIdxp12) == 0 
			// 			 && fitMuCharge(muIdxp21) + fitMuCharge(muIdxp22) == 0
			// 		  )
			// 		{
			// 			// Modify the DiMuonMass expression appropriatly. 
			// 			// Use the fitMup4vect and the muIdxpXY indexes defined above.
			// 			DiMuonMass1 = 0; 
			// 			DiMuonMass2 = 0;
			// 			myDiMuon1mass->Fill(DiMuonMass1);
			// 			myDiMuon2mass->Fill(DiMuonMass2);

						// if (1
						// 		// Here require that each DiMuonMass is in the appropriate mass range [2.95,3.25] GeV
						// 								   )
						// {
						// 	// calculate the 4 muon mass:  M(µ1µ2µ3µ4)-M(µ1µ2)-M(µ3µ4)+2*M(J/psi)  
						// 	m4Muon = 0;

						// 	myFourMuonmass->Fill(m4Muon);

					

						// }
					// }

				// }
			

	cout << "I have done " << nentries << " total entries" << endl;
	myhbk->Write();

	TFile *outFile = new TFile("myntuple_hists.root", "RECREATE");

	myDiMuon1mass->Write("myDiMuon1mass");
	myDiMuon2mass->Write("myDiMuon2mass");
	myFourMuonmass->Write("myFourMuonmass");
	myDiMuon1mass_PT->Write("myDiMuon1mass_PT");
	myDiMuon2mass_PT->Write("myDiMuon2mass_PT");
	myFourMuonmass_PT->Write("myFourMuonmass_PT");
	myDiMuon1mass_eta->Write("myDiMuon1mass_eta");
	myDiMuon2mass_eta->Write("myDiMuon2mass_eta");
	myFourMuonmass_eta->Write("myFourMuonmass_eta");
	myDiMuon1mass_chg->Write("myDiMuon1mass_chg");
	myDiMuon2mass_chg->Write("myDiMuon2mass_chg");
	myFourMuonmass_chg->Write("myFourMuonmass_chg");
	myDiMuon1mass_mass->Write("myDiMuon1mass_mass");
	myDiMuon2mass_mass->Write("myDiMuon2mass_mass");
	myFourMuonmass_mass->Write("myFourMuonmass_mass");

	outFile->Close();

}

