#define myntuple_cxx
#include "myntuple.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TLorentzVector.h"
#include <iostream> 
#include <fstream>
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
			// std::cout << "no cuts" <<std::endl;
			if (1
					// soft muon: tracker muon + 1 hit in the muon system 
					&& myNumPatSoftMuon >= 4 
					// requiring all muons to come from the same vertex
					&& (*MyFourMuonVtxCL)[myFourMuIdx] >= 0.005     
					// Here add the selections! 
                    && (TrigThreeMuonJpsi3p5mu2 == 1 || TrigThreeMuonJpsi == 1)
					
					&& (fitMuCharge[0] + fitMuCharge[1] + fitMuCharge[2] + fitMuCharge[3] == 0)

                    && rawMup4vect[0].Pt() >= 2.
                    && rawMup4vect[1].Pt() >= 2.
                    && rawMup4vect[2].Pt() >= 2.
                    && rawMup4vect[3].Pt() >= 2.

                    && abs(rawMup4vect[0].Eta()) <= 2.4
                    && abs(rawMup4vect[1].Eta()) <= 2.4
                    && abs(rawMup4vect[2].Eta()) <= 2.4
                    && abs(rawMup4vect[3].Eta()) <= 2.4

				) {
					// std::cout << "first pass" <<std::endl;
				for (int mypidx = 0; mypidx < 3; mypidx++)  {
					int muIdxp11, muIdxp12, muIdxp21, muIdxp22;
					muIdxp11 = myCombIdx[mypidx].p11; muIdxp12 = myCombIdx[mypidx].p12; muIdxp21 = myCombIdx[mypidx].p21; muIdxp22 = myCombIdx[mypidx].p22;


					if(1
                        // Here, require the muon pairs to have muons with opposite charges
                        && fitMuCharge[muIdxp11] + fitMuCharge[muIdxp12] == 0 
						&& fitMuCharge[muIdxp21] + fitMuCharge[muIdxp22] == 0
					  )
					{
						// std::cout << "second pass" <<std::endl;
						// Modify the DiMuonMass expression appropriatly. 
						// Use the fitMup4vect and the muIdxpXY indexes defined above.
                        DiMuonMass1 = (fitMup4vect[muIdxp11] + fitMup4vect[muIdxp12]).M();
						DiMuonMass2 = (fitMup4vect[muIdxp21] + fitMup4vect[muIdxp22]).M();
						// calculate the 4 muon mass:  M(µ1µ2µ3µ4)-M(µ1µ2)-M(µ3µ4)+2*M(J/psi)  
						

						myDiMuon1mass->Fill(DiMuonMass1);
						myDiMuon2mass->Fill(DiMuonMass2);

						if (1
								// Here require that each DiMuonMass is in the appropriate mass range [2.95,3.25] GeV
                                && DiMuonMass1 <= 3.25
                                && DiMuonMass1 >= 2.95
                                && DiMuonMass2 <= 3.25
                                && DiMuonMass2 >= 2.95
														   )
						{
							// calculate the 4 muon mass:  M(µ1µ2µ3µ4)-M(µ1µ2)-M(µ3µ4)+2*M(J/psi)  
							m4Muon 	= (
							(fitMup4vect[muIdxp11] + fitMup4vect[muIdxp12] + fitMup4vect[muIdxp21] + fitMup4vect[muIdxp22]).M()
							- DiMuonMass1 - DiMuonMass2
							+ 2*JPSI_MASS );

							myFourMuonmass->Fill(m4Muon);

							myoutputfile
								<< std::fixed
								<< m4Muon
								<< endl;

						}
					}

				}
			}

		}

	}

	cout << "I have done " << nentries << " total entries" << endl;
	myhbk->Write();
}
