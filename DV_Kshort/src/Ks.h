#ifndef DDLSTUDIES_KSHORT_DDL_H
#define DDLSTUDIES_KSHORT_DDL_H 

#include "AthAnalysisBaseComps/AthAnalysisAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/ITHistSvc.h"
 


#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/Vertex.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"

using namespace std;

// forward declaration
class TTree;
class TEfficiency;



namespace DDL {

    class Ks: public AthAnalysisAlgorithm 
    {
 
		private: 
			// histogram service
        		//ServiceHandle<ITHistSvc> m_histSvc;
        		//std::string m_hist_name;

			int event_counter=0;

			TTree* m_kshort_tree = nullptr;
			ITHistSvc * m_thistSvc = nullptr;
			StoreGateSvc* m_storeGate = nullptr;

			//variables for the TTree
			Double_t m_piplus_pt = 0; 	// transverse momentum of pi+ track
			Double_t m_piplus_p  = 0; 	// momentum magnitude of pi+ track
			Double_t m_piplus_px = 0;	// x component of momentum of pi+ track
			Double_t m_piplus_py = 0;	// y component of momentum of pi+ track
			Double_t m_piplus_pz = 0;	// z component of momentum of pi+ track
			Double_t m_piplus_e  = 0;	// energy of pi+ track
			Double_t m_piplus_z0 = 0;	// z0 of pi+ track 
			Double_t m_piplus_d0 = 0;	// d0 of pi+ track
			Double_t m_piplus_eta= 0;	// pseudo rapidity of pi+ track 

			Double_t m_piminus_pt= 0;	// transverse momentum of pi- track
			Double_t m_piminus_p = 0;	// momentum magnitude of pi- track
			Double_t m_piminus_px= 0;	// x component of momentum of pi- track
			Double_t m_piminus_py= 0;	// y component of momentum of pi- track
			Double_t m_piminus_pz= 0;	// z component of momentum of pi- track
			Double_t m_piminus_e = 0;	// energy of pi- track
			Double_t m_piminus_z0= 0;	// z0 of pi- track 
			Double_t m_piminus_d0= 0;	// d0 of pi- track
			Double_t m_piminus_eta=0;	// pseudo rapidity of pi- track 
	
			Double_t m_kshort_mass=0;	// the mass of the kshort particle

        	public: 
            		Ks( const std::string& name, ISvcLocator* pSvcLocator );
            		virtual ~Ks() = default; 

            		StatusCode  initialize() override;
            		StatusCode  execute() override;
            		StatusCode  finalize() override;

            		// for good Ks ?
            		//virtual bool goodKs () const;
			bool isPi(float piMass);
			StatusCode finding_right_ks();
	};

} 

#endif // DDLSTUDIES_KSHORT_DDL_H  
