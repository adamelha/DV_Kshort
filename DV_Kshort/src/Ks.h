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
#include <vector>
#include <TLorentzVector.h>

using namespace std;

// forward declaration
class TTree;
class TEfficiency;



namespace DDL 
{

    class Ks: public AthAnalysisAlgorithm 
    {
 
		private: 
			// histogram service
        		//ServiceHandle<ITHistSvc> m_histSvc;
        		//std::string m_hist_name;

			// Event number
			int m_event_counter		 = 0;
			// TTree
			TTree* m_kshort_tree 		 = nullptr;
			// Description??			
			ITHistSvc * m_thistSvc 		 = nullptr;
			StoreGateSvc* m_storeGate 	 = nullptr;
			// Description?? 
			size_t m_total_rec_ks 		 = 0;
			size_t m_total_rec_ks_hits 	 = 0;
			// The PDG number for Ks
			const Int_t m_KS_PDG		 = 310;			
			// Sizes of MC (simulated) and real data samples
			size_t m_sim_counter 		 = 0;
			size_t m_data_counter 		 = 0;


			// Variables of pi+ for the TTree (from Secondary vertices container)
			Double_t m_piplus_pt 		 = 0; // Transverse momentum of pi+ track
			Double_t m_piplus_p  		 = 0; // Momentum magnitude of pi+ track
			Double_t m_piplus_px 		 = 0; // x component of momentum of pi+ track
			Double_t m_piplus_py 		 = 0; // y component of momentum of pi+ track
			Double_t m_piplus_pz 		 = 0; // z component of momentum of pi+ track
			Double_t m_piplus_e  		 = 0; // Energy of pi+ track
			Double_t m_piplus_z0 		 = 0; // z0 of pi+ track 
			Double_t m_piplus_d0 		 = 0; // d0 of pi+ track
			Double_t m_piplus_eta		 = 0; // Pseudo rapidity of pi+ track 
			// Variables of pi- for the TTree (from Secondary vertices container)
			Double_t m_piminus_pt		 = 0; // Transverse momentum of pi- track
			Double_t m_piminus_p 		 = 0; // Momentum magnitude of pi- track
			Double_t m_piminus_px		 = 0; // x component of momentum of pi- track
			Double_t m_piminus_py		 = 0; // y component of momentum of pi- track
			Double_t m_piminus_pz		 = 0; // z component of momentum of pi- track
			Double_t m_piminus_e 		 = 0; // Energy of pi- track
			Double_t m_piminus_z0		 = 0; // z0 of pi- track 
			Double_t m_piminus_d0		 = 0; // d0 of pi- track
			Double_t m_piminus_eta		 = 0; // Pseudo rapidity of pi- track 
			// Variables of Ks for the TTree (from Secondary vertices container)
			Double_t m_kshort_mass	  	 = 0; // Invariant mass of Ks vertex (To be used, the correct one)
			Double_t m_kshort_invMass 	 = 0; // Indirect invariant mass calculation of Ks using 4-momenta of tracks constituting its vertex (i.e. pion tracks). NOT to be used
			Double_t m_kshort_rDV	  	 = 0; // r_DV of Ks (distance between primary(0,0,0) and secondary vertices of Ks in xy-plane) 
			Double_t m_kshort_theta	  	 = 0; // Theta calculation with Ks vertex (angle between Ks momentum and beam axis [which in the direction of +ve z axis]) 
			Double_t m_kshort_eta	  	 = 0; // Pseudo rapidity of Ks (calculated using 'Theta')
			Double_t m_kshort_e	  	 = 0; // Energy of Ks (calculated using the info. of Ks mass and momentum)
			Double_t m_kshort_pt	  	 = 0; // Transverse momemntum of Ks (NOT to be used, general addition of the pT of individual tracks)
			Double_t m_kshort_p	  	 = 0; // Momentum magnitude of Ks
			Double_t m_kshort_px	  	 = 0; // x component of momentum of Ks
			Double_t m_kshort_py	  	 = 0; // y component of momentum of Ks
			Double_t m_kshort_pz	  	 = 0; // z component of momentum of Ks
			Double_t m_kshort_x	  	 = 0; // x coordinate of Ks vertex
			Double_t m_kshort_y	  	 = 0; // y coordinate of Ks vertex
			Double_t m_kshort_z	  	 = 0; // z coordinate of Ks vertex
			Double_t m_kshort_alpha   	 = 0; // Alpha calculation with Ks vertex (the angle between r_DV [which is on xy plane] and momentum of Ks) 
			Double_t m_kshort_modified_alpha = 0; // Alpha calculation with Ks vertex (the angle between r_DV-rpv_t (of most energetic PV) and momentum of Ks)			
			Double_t m_kshort_pTCalc  	 = 0; // Transverse momentum calculation using px and py of Ks (To be used, the correct one)
 			Double_t m_z_sv_pv		 = 0; // The distance between the first primary vertex and the extrapolation of the momentum of the secondary vertex along z-axis
			Double_t m_kshort_phi		 = 0; // The angle between rdv and the x axis.


			/*
			// Variables related to error matrix of Ks vertices (from Secondary vertices container)
			// It's a 3x3 matrix, and all its nine elements are saved individually along with the square root of the diagonal elements 
			Double_t  m_covariance00; // Element in 0th row and 0th column (diagonal element)
			Double_t  m_covariance01; // Element in 0th row and 1st column
			Double_t  m_covariance02; // Element in 0th row and 2nd column
			Double_t  m_covariance10; // Element in 1st row and 0th column
			Double_t  m_covariance11; // Element in 1st row and 1st column (diagonal element)
			Double_t  m_covariance12; // Element in 1st row and 2nd column
			Double_t  m_covariance20; // Element in 2nd row and 0th column
			Double_t  m_covariance21; // Element in 2nd row and 1st column
			Double_t  m_covariance22; // Element in 2nd row and 2nd column (diagonal element)
			Double_t  m_sr_covariance00; // sqrt(Element in 0th row and 0th column (diagonal element))
			Double_t  m_sr_covariance11; // sqrt(Element in 1st row and 1st column (diagonal element))
			Double_t  m_sr_covariance22; // sqrt(Element in 2nd row and 2nd column (diagonal element))
			*/

			// Variables related to the most energetic primary vertices 
			Double_t m_primary_vertex_x 	 = 0; // x position of primary vertex (from Primary vertex container, which has different size from the secondary one)
			Double_t m_primary_vertex_y 	 = 0; // y position of primary vertex (from Primary vertex container, which has different size from the secondary one)
			Double_t m_primary_vertex_z 	 = 0; // z position of primary vertex (from Primary vertex container, which has different size from the secondary one)
			Double_t m_primary_vertex_pt	 = 0; // Sum of transverse momenta of all the tracks of the most energetic primary vertex 
			Int_t m_most_energetic_primary_vertex_index  = 0; // index of the most energetic primary vertex

			// Variables related to the first primary vertices	
			Double_t m_primary_vertex_0x     = 0; // x position of the first primary vertex (from Primary vertex container)
			Double_t m_primary_vertex_0y     = 0; // y position of the first primary vertex (from Primary vertex container)
			Double_t m_primary_vertex_0z     = 0; // z position of the first primary vertex (from Primary vertex container)
			Double_t m_zeroth_primary_vertex_pt    = 0; // Sum of transverse momenta of all tracks of the first primary vertex		

			// Variables related to Truth checking
			Int_t m_truth_piplus_pdgid	 = 0; // PDG id of pi+ (If TRUE, the numerical value should be equal to 211)
			Int_t m_truth_piminus_pdgid	 = 0; // PDG id of pi- (If TRUE, the numerical value should be equal to -211)
			Int_t m_truth_kshort_pdgid	 = 0; // PDG id of Ks  (If TRUE, the numerical value should be equal to 310)


        	public: 
            		Ks( const std::string& name, ISvcLocator* pSvcLocator );
            		virtual ~Ks() = default; 

            		StatusCode  initialize() override;
            		StatusCode  execute() override;
            		StatusCode  finalize() override;

            		// For mass of the charged pions that constitute the Ks vertex
			bool isPi(float piMass);
			// Finding right Ks w/ the Truth information
			StatusCode finding_right_ks();
			Int_t get_most_energetic_vertex_index(const xAOD::VertexContainer* vertices);
			Double_t get_sum_pt(const xAOD::Vertex* vertex_ptr);



	};

} 

#endif // DDLSTUDIES_KSHORT_DDL_H  
