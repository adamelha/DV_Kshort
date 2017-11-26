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

#include <TLorentzVector.h>

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
			TTree* m_primary_vertices_tree = nullptr;

			ITHistSvc * m_thistSvc = nullptr;
			StoreGateSvc* m_storeGate = nullptr;

			// Variables of pi+ for the TTree
			Double_t m_piplus_pt 	= 0; 	// Transverse momentum of pi+ track
			Double_t m_piplus_p  	= 0; 	// Momentum magnitude of pi+ track
			Double_t m_piplus_px 	= 0;	// x component of momentum of pi+ track
			Double_t m_piplus_py 	= 0;	// y component of momentum of pi+ track
			Double_t m_piplus_pz 	= 0;	// z component of momentum of pi+ track
			Double_t m_piplus_e  	= 0;	// Energy of pi+ track
			Double_t m_piplus_z0 	= 0;	// z0 of pi+ track 
			Double_t m_piplus_d0 	= 0;	// d0 of pi+ track
			Double_t m_piplus_eta	= 0;	// Pseudo rapidity of pi+ track 
			// Variables of pi- for the TTree
			Double_t m_piminus_pt	= 0;	// Transverse momentum of pi- track
			Double_t m_piminus_p 	= 0;	// Momentum magnitude of pi- track
			Double_t m_piminus_px	= 0;	// x component of momentum of pi- track
			Double_t m_piminus_py	= 0;	// y component of momentum of pi- track
			Double_t m_piminus_pz	= 0;	// z component of momentum of pi- track
			Double_t m_piminus_e 	= 0;	// Energy of pi- track
			Double_t m_piminus_z0	= 0;	// z0 of pi- track 
			Double_t m_piminus_d0	= 0;	// d0 of pi- track
			Double_t m_piminus_eta	= 0;	// Pseudo rapidity of pi- track 
			// Variables of Ks for the TTree
			Double_t m_kshort_mass	  = 0;	// Invariant mass of Ks vertex (To be used, the correct one)
			Double_t m_kshort_invMass = 0;  // Indirect invariant mass calculation of Ks using 4-momenta of tracks constituting its vertex (i.e. pion tracks). NOT to be used
			Double_t m_kshort_rDV	  = 0;  // r_DV of Ks (distance between primary(0,0,0) and secondary vertices of Ks in xy-plane) 
			Double_t m_kshort_theta	  = 0;  // Theta calculation with Ks vertex (angle between Ks momentum and beam axis [which in the direction of +ve z axis]) 
			Double_t m_kshort_eta	  = 0;  // Pseudo rapidity of Ks (calculated using 'Theta')
			Double_t m_kshort_e	  = 0;  // Energy of Ks (calculated using the info. of Ks mass and momentum)
			Double_t m_kshort_pt	  = 0;  // Transverse momemntum of Ks (NOT to be used, general addition of the pT of individual tracks)
			Double_t m_kshort_p	  = 0;  // Momentum magnitude of Ks
			Double_t m_kshort_px	  = 0;  // x component of momentum of Ks
			Double_t m_kshort_py	  = 0;  // y component of momentum of Ks
			Double_t m_kshort_pz	  = 0;  // z component of momentum of Ks
			Double_t m_kshort_alpha   = 0;	// Alpha calculation with Ks vertex (the angle between r_DV [which is on xy plane] and momentum of Ks) 			
			Double_t m_kshort_pTCalc  = 0;  // Transverse momentum calculation using px and py of Ks (To be used, the correct one)
 			// Variables related to primary vertices, stored in a different TTree
			Double_t m_primary_vertex_x = 0;// x position of Ks primary vertex (from Primary vertex container, which has different size from the secondary one)
			Double_t m_primary_vertex_y = 0;// y position of Ks primary vertex (from Primary vertex container, which has different size from the secondary one)
			Double_t m_primary_vertex_z = 0;// z position of Ks primary vertex (from Primary vertex container, which has different size from the secondary one)


        	public: 
            		Ks( const std::string& name, ISvcLocator* pSvcLocator );
            		virtual ~Ks() = default; 

            		StatusCode  initialize() override;
            		StatusCode  execute() override;
            		StatusCode  finalize() override;

            		// For mass of the charged pions that constitute the Ks vertex
			bool isPi(float piMass);
			StatusCode finding_right_ks();
	};

} 

#endif // DDLSTUDIES_KSHORT_DDL_H  
