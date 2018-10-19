#include <TMatrixD.h>
#include <TTree.h>

class p0dCCEvent {
	public:
		TMatrixD* WeightsMatrix;//15,400
		Int_t topology;
		Int_t nu_pdg;
		Float_t nu_trueE;
		Float_t truelepton_mom;
		Float_t truelepton_costheta;
		Float_t* vtx_truepos;//4
		Int_t IsOnWater;
		Float_t weightHL;

		p0dCCEvent() {
			WeightsMatrix = new TMatrixD(15,400);
			topology=-1;
			nu_pdg=-1;
			nu_trueE=-999.;
			truelepton_mom=-999.;
			truelepton_costheta=-999.;
			vtx_truepos = new Float_t[4];
			IsOnWater=-1;
			weightHL=-1;
		}
		~p0dCCEvent() {
			delete []vtx_truepos;
			//TODO how to delete: delete WeightsMatrixInit;
		}

		void SetBranchAddresses(TTree* inT, bool isMC, bool isGENIE) {
			if(!isGENIE){
				inT->SetBranchAddress("WeightsMatrix",&WeightsMatrix);
				inT->SetBranchAddress("weight",&weightHL);
			}
			inT->SetBranchAddress("topology",&topology);
			inT->SetBranchAddress("nu_pdg",&nu_pdg);
			inT->SetBranchAddress("nu_trueE",&nu_trueE);
			inT->SetBranchAddress("truelepton_mom",&truelepton_mom);
			inT->SetBranchAddress("truelepton_costheta",&truelepton_costheta);
			inT->SetBranchAddress("IsOnWater",&IsOnWater);
			if(isMC){
				inT->SetBranchAddress("truevtx_pos",vtx_truepos);
			} else {
				inT->SetBranchAddress("selvtx_truepos",vtx_truepos);
			}
		}
};
