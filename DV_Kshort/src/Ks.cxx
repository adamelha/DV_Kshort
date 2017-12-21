// DDLStudies includes
#include "Ks.h"

// framework includes
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ITHistSvc.h"
				
// EDM includes: - if move to header file will not compile?
#include "xAODEventInfo/EventInfo.h"
//#include "xAODTracking/TrackMeasurementValidationContainer.h"
//#include "xAODTracking/TrackStateValidationContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
//#include "xAODTracking/TrackParticleAuxContainer.h"
#include "xAODTracking/VertexContainer.h"

#include "xAODTracking/Vertex.h"
//#include "xAODTracking/VertexAuxContainer.h"

#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODTruth/TruthVertexContainer.h"
#include "xAODTruth/xAODTruthHelpers.h"


//ROOT
#include "TEfficiency.h"
#include "TFile.h"


#include "GaudiKernel/ITHistSvc.h"

#include<iostream>
#include<bitset>
#include<vector>
#include<string>
#include <stdint.h>
#include <algorithm>
#include <math.h>
#include <functional>
#include <cmath>



DDL::Ks::Ks( const std::string& name, ISvcLocator* pSvcLocator ) : 
   AthAnalysisAlgorithm( name, pSvcLocator ) 
{

  //declareProperty("HistSvcName", m_hist_name = "Ks" ); 

}

StatusCode DDL::Ks::initialize() 
{ 
	std::cout << "Initializing" << std::endl;

	StatusCode sc = service("StoreGateSvc", m_storeGate);
	if (sc.isFailure())
	{
		//msg(MSG::ERROR) << "Unable to retrieve pointer to StoreGateSvc" << endreq;
		std::cout << "Unable to retrieve pointer to StoreGateSvc" << std::endl;
		return sc;
	}

	sc = service("THistSvc", m_thistSvc);
	if (sc.isFailure())
	{
		//msg(MSG::ERROR) << "Unable to retrieve pointer to THistSvc" << endreq;
		std::cout << "Unable to retrieve pointer to THistSvc" << std::endl;
		return sc;
	}
	// Histograming...
	ServiceHandle<ITHistSvc> histSvc("THistSvc", name());
	CHECK(histSvc.retrieve());

	if (sc.isFailure())
	{
		//msg(MSG::ERROR) << "Unable to register histogramming service" << endreq;
		std::cout << "Unable to register histogramming service" << std::endl;
		return sc;
	}


	//// ATH_MSG_INFO ("Initializing " << name() << "...");


	m_kshort_tree = new TTree("KsTTree", "Ks Tree");
	CHECK(histSvc->regTree("/Ks/m_kshort_tree", m_kshort_tree));
	// For detailed description of the branches, please refer to Ks.h (the header file)	
        // Branches for pi+
	m_kshort_tree->Branch("piplus_pt",   &m_piplus_pt,   "piplus_pt/D"  );
	m_kshort_tree->Branch("piplus_p",    &m_piplus_p,    "piplus_p/D"   );
	m_kshort_tree->Branch("piplus_px",   &m_piplus_px,   "piplus_px/D"  );
	m_kshort_tree->Branch("piplus_py",   &m_piplus_py,   "piplus_py/D"  );
	m_kshort_tree->Branch("piplus_pz",   &m_piplus_pz,   "piplus_pz/D"  );
	m_kshort_tree->Branch("piplus_e",    &m_piplus_e,    "piplus_e/D"   );
	m_kshort_tree->Branch("piplus_z0",   &m_piplus_z0,   "piplus_z0/D"  );
	m_kshort_tree->Branch("piplus_d0",   &m_piplus_d0,   "piplus_d0/D"  );
	m_kshort_tree->Branch("piplus_eta",  &m_piplus_eta,  "piplus_eta/D" );
	// Branches for pi-
	m_kshort_tree->Branch("piminus_pt",  &m_piminus_pt,  "piminus_pt/D" );
	m_kshort_tree->Branch("piminus_p",   &m_piminus_p,   "piminus_p/D"  );
	m_kshort_tree->Branch("piminus_px",  &m_piminus_px,  "piminus_px/D" );
	m_kshort_tree->Branch("piminus_py",  &m_piminus_py,  "piminus_py/D" );
	m_kshort_tree->Branch("piminus_pz",  &m_piminus_pz,  "piminus_pz/D" );
	m_kshort_tree->Branch("piminus_e",   &m_piminus_e,   "piminus_e/D"  );
	m_kshort_tree->Branch("piminus_z0",  &m_piminus_z0,  "piminus_z0/D" );
	m_kshort_tree->Branch("piminus_d0",  &m_piminus_d0,  "piminus_d0/D" );
	m_kshort_tree->Branch("piminus_eta", &m_piminus_eta, "piminus_eta/D");
	// Branches for Ks
	m_kshort_tree->Branch("kshort_mass",   &m_kshort_mass,   "kshort_mass/D"   );
	m_kshort_tree->Branch("kshort_rDV",    &m_kshort_rDV,    "kshort_rDV/D"    );
	m_kshort_tree->Branch("kshort_theta",  &m_kshort_theta,  "kshort_theta/D"  );
	m_kshort_tree->Branch("kshort_eta",    &m_kshort_eta,    "kshort_eta/D"    );
	m_kshort_tree->Branch("kshort_e",      &m_kshort_e,      "kshort_e/D"      );
	m_kshort_tree->Branch("kshort_pt",     &m_kshort_pt,     "kshort_pt/D"     );
	m_kshort_tree->Branch("kshort_p",      &m_kshort_p,      "kshort_p/D"      );
	m_kshort_tree->Branch("kshort_px",     &m_kshort_px,     "kshort_px/D"     );
	m_kshort_tree->Branch("kshort_py",     &m_kshort_py,     "kshort_py/D"     );
	m_kshort_tree->Branch("kshort_pz",     &m_kshort_pz,     "kshort_pz/D"     );
        m_kshort_tree->Branch("kshort_x",      &m_kshort_x,      "kshort_x/D"      );
	m_kshort_tree->Branch("kshort_y",      &m_kshort_y,      "kshort_y/D"      );
	m_kshort_tree->Branch("kshort_z",      &m_kshort_z,      "kshort_z/D"      );
	m_kshort_tree->Branch("kshort_pTCalc", &m_kshort_pTCalc, "kshort_pTCalc/D" );
	m_kshort_tree->Branch("kshort_invMass",&m_kshort_invMass,"kshort_invMass/D");
	m_kshort_tree->Branch("kshort_alpha",  &m_kshort_alpha,  "kshort_alpha/D"  );
	// Brnaches for primary vertex
	m_kshort_tree->Branch("primary_vertex_x", &m_primary_vertex_x, "primary_vertex_x/D" );
	m_kshort_tree->Branch("primary_vertex_y", &m_primary_vertex_y, "primary_vertex_y/D" );
	m_kshort_tree->Branch("primary_vertex_z", &m_primary_vertex_z, "primary_vertex_z/D" );
	m_kshort_tree->Branch("primary_vertex_pt",&m_primary_vertex_pt,"primary_vertex_pt/D");
	m_kshort_tree->Branch("most_energetic_primary_vertex_index",&m_most_energetic_primary_vertex_index , "most_energetic_primary_vertex_index/I"); 
	// Branches for Truth checking
	m_kshort_tree->Branch("truth_piplus_pdgid", &m_truth_piplus_pdgid, "truth_piplus_pdgid/I" );
	m_kshort_tree->Branch("truth_piminus_pdgid",&m_truth_piminus_pdgid,"truth_piminus_pdgid/I");
	m_kshort_tree->Branch("truth_kshort_pdgid", &m_truth_kshort_pdgid, "truth_kshort_pdgid/I" );
	//Branches for error matrix related to Ks vertices
	m_kshort_tree->Branch("covariance00", &m_covariance00, "covariance00/D");
	m_kshort_tree->Branch("covariance01", &m_covariance01, "covariance01/D");
	m_kshort_tree->Branch("covariance02", &m_covariance02, "covariance02/D");
	m_kshort_tree->Branch("covariance10", &m_covariance10, "covariance10/D");
	m_kshort_tree->Branch("covariance11", &m_covariance11, "covariance11/D");
	m_kshort_tree->Branch("covariance12", &m_covariance12, "covariance12/D");
	m_kshort_tree->Branch("covariance20", &m_covariance20, "covariance20/D");
	m_kshort_tree->Branch("covariance21", &m_covariance21, "covariance21/D");
	m_kshort_tree->Branch("covariance22", &m_covariance22, "covariance22/D");
	// Related to both Ks and primary vertices 
	m_kshort_tree->Branch("z_sv_pv",&m_z_sv_pv,"z_sv_pv/D");


	return StatusCode::SUCCESS;
}

StatusCode DDL::Ks::finalize() 
{

	// ATH_MSG_INFO ("Finalizing " << name() << "...");

	std::cout << "Total number of K Short particles reconstructed: " << m_total_rec_ks << "\n";
	if(m_sim_counter>0)
	{
		std::cout << "Out of the reconstructed K Short particles " << m_total_rec_ks_hits << " were actually true\n";
		if(m_total_rec_ks>0)
		{
			std::cout << "Success rate: " << (double)m_total_rec_ks_hits / (double)m_total_rec_ks << "\n";
		}
  	}
	std::cout << "Simulation counter = " << m_sim_counter << ", Data counter = " << m_data_counter << std::endl;
	return StatusCode::SUCCESS;
}

StatusCode DDL::Ks::execute() 
{
	this->m_event_counter++;
	std::cout << "Event Number : " << this->m_event_counter << std::endl;

	const xAOD::EventInfo* ei_ptr = 0;
  	//evtStore() returns the EventStore of the current event in the event loop. it represents the event.
    	CHECK(evtStore()->retrieve(ei_ptr,"EventInfo"));
  	//we want to check if the current event is a Data event or a Monte Carlo Simulation event 
  	if(ei_ptr->eventType(xAOD::EventInfo::IS_SIMULATION))
  	{
		std::cout << "Simulation" << std::endl;
		m_sim_counter++;
		const xAOD::TruthVertexContainer* truthVertexContainer = nullptr;
		CHECK(m_storeGate->retrieve(truthVertexContainer, "TruthVertices"));
		std::cout << "Truth Vertex Container size = " << truthVertexContainer->size() << std::endl;

  	}
  	else
  	{
		std::cout << "Data" << std::endl;
		m_data_counter++;
  	} 


	StatusCode sc = StatusCode::SUCCESS;  
	sc = finding_right_ks();
	if (sc.isFailure())
	{
		msg(MSG::ERROR) << "Finding right ks failed" << endreq;
		return StatusCode::SUCCESS;
	}
	return sc;
  	// ATH_MSG_DEBUG ("Executing " << name() << "...");
}


// Function for charged pion mass (equating float to some number inside a condition is not recommended, hence this function was made) 
bool DDL::Ks::isPi(float piMass)
{
        float actual_pi_mass = 139.57061; // MeV (taken from recent PDG, http://pdglive.lbl.gov/Particle.action?init=0&node=S008)
        //float margin_of_error = 0.00024; // MeV (taken from recent PDG, http://pdglive.lbl.gov/Particle.action?init=0&node=S008)
		float margin_of_error = 0.02; // MeV (taken from recent PDG, http://pdglive.lbl.gov/Particle.action?init=0&node=S008)
        return ( (piMass >= (actual_pi_mass - margin_of_error)) &&  (piMass <= (actual_pi_mass + margin_of_error)) );
}

Int_t DDL::Ks::get_most_energetic_vertex_index(const xAOD::VertexContainer* vertices)
{
	
	size_t max_vertex_index=0, i, num_of_vertices;
	Double_t max_vertex_pt;
	Double_t current_vertex_pt;

	if(vertices != nullptr && vertices->size()>0)
	{
		num_of_vertices = vertices->size();
		max_vertex_pt = get_sum_pt((*vertices)[0]);
		for(i=1; i<num_of_vertices; i++ )
		{
			current_vertex_pt = get_sum_pt((*vertices)[i]);
			if(current_vertex_pt > max_vertex_pt)
			{
				max_vertex_pt = current_vertex_pt;
				max_vertex_index = i;
			}
		}
	}
	return max_vertex_index;
	
}

Double_t DDL::Ks::get_sum_pt(const xAOD::Vertex* vertex_ptr)
{
	Double_t sum_pt = 0;
	size_t i, num_of_tracks;
	if(vertex_ptr != nullptr)
	{
		num_of_tracks = vertex_ptr->nTrackParticles();
		for (i = 0; i < num_of_tracks; i++) 
		{
			sum_pt += vertex_ptr->trackParticle(i)->pt();
		}
	}
	return sum_pt;	
}


StatusCode DDL::Ks::finding_right_ks()
{

	StatusCode sc = StatusCode::SUCCESS;
	//----------------------------------------------------------------------------

	const xAOD::VertexContainer* primary_vertices = nullptr;
	const xAOD::VertexContainer* secondary_vertices = nullptr;

	CHECK(evtStore()->retrieve(primary_vertices, "PrimaryVertices")); // For primary vertices
	CHECK(evtStore()->retrieve(secondary_vertices , "VrtSecInclusive_SecondaryVertices")); // For secondary vertices

	// "TrackParticle": pointer to a given track that was used in vertex reconstruction
	const xAOD::TrackParticle* piplus_track = nullptr;
	const xAOD::TrackParticle* piminus_track = nullptr;

	TLorentzVector p4_sum;
	Double_t rdvDotPt;
	Double_t avg_num_of_tracks = 0;

	// Pointer for the first primary vertex from the primary vertices container
	const xAOD::Vertex* most_energetic_vertex = nullptr;
	
	if(primary_vertices->size()>0)
	{	
		
		m_most_energetic_primary_vertex_index = get_most_energetic_vertex_index(primary_vertices);
		most_energetic_vertex = (*primary_vertices)[m_most_energetic_primary_vertex_index];

		
		// Position coordinates of the first primary vertex
		m_primary_vertex_x = most_energetic_vertex ->x();
		m_primary_vertex_y = most_energetic_vertex ->y();
		m_primary_vertex_z = most_energetic_vertex ->z();
		
		m_primary_vertex_pt = get_sum_pt(most_energetic_vertex);
		
		for(const xAOD::Vertex* vertex_ptr : *primary_vertices)
		{
			avg_num_of_tracks += vertex_ptr->nTrackParticles();
		}
		avg_num_of_tracks = avg_num_of_tracks/(primary_vertices->size());
		std::cout << "Average tracks number for primary vertices : " << avg_num_of_tracks << std::endl;
	}


	// Looping over the vertices to find Kshort vertices which decay to pi+ and pi-
	// Defining the iterator in a short and efficient way
	for (xAOD::Vertex* vertex_ptr : *secondary_vertices )
	{
		// Vertices with only two tracks 
		if (vertex_ptr->nTrackParticles() == 2)
		{
			// K short pidType is 310, (or PDG::pidType::K_S0)
			piplus_track = vertex_ptr->trackParticle(0);
			piminus_track = vertex_ptr->trackParticle(1);

			if (piplus_track->charge() < piminus_track->charge()) //in case track#0 is pi- and track#1 is pi+
			{
				std::swap(piminus_track, piplus_track);
			}
			// Checking whether the masses are equal to the charged pion mass or not
			if (isPi(piplus_track->m()) && isPi(piminus_track->m()))
			{
				// Checking if the tracks have charges of -1 and +1
				if (piplus_track->charge() + piminus_track->charge() == 0 && piplus_track->charge() == 1)
				{
					
					// We decide that these pions have originated from Ks
					m_total_rec_ks++;
					
					
					//pi+ track info
					m_piplus_pt   = piplus_track->pt();
					m_piplus_p    = piplus_track->p4().P();
					m_piplus_px   = piplus_track->p4().Px();
					m_piplus_py   = piplus_track->p4().Py();
					m_piplus_pz   = piplus_track->p4().Pz();
					m_piplus_e    = piplus_track->e();
					m_piplus_z0   = piplus_track->z0();
					m_piplus_d0   = piplus_track->d0();
					m_piplus_eta  = piplus_track->eta();
					// pi- track info
					m_piminus_pt  = piminus_track->pt();
					m_piminus_p   = piminus_track->p4().P();
					m_piminus_px  = piminus_track->p4().Px();
					m_piminus_py  = piminus_track->p4().Py();
					m_piminus_pz  = piminus_track->p4().Pz();
					m_piminus_e   = piminus_track->e();
					m_piminus_z0  = piminus_track->z0();
					m_piminus_d0  = piminus_track->d0();
					m_piminus_eta = piminus_track->eta();
					// Info. related to both pi+ and pi tracks
					p4_sum = piplus_track->p4() + piminus_track->p4();
					m_kshort_invMass = p4_sum.Mag();
					// Ks info
					m_kshort_mass = vertex_ptr->auxdataConst<float>("vtx_mass");
					m_kshort_rDV  = sqrt(pow(vertex_ptr->x(),2) + pow(vertex_ptr->y(),2)); 
					m_kshort_e  = sqrt(pow(m_kshort_mass ,2)+pow(m_kshort_p ,2)); 
					m_kshort_px = vertex_ptr->auxdataConst<float>("vtx_px");
					m_kshort_py = vertex_ptr->auxdataConst<float>("vtx_py");
					m_kshort_pz = vertex_ptr->auxdataConst<float>("vtx_pz");
					m_kshort_pt = vertex_ptr->auxdataConst<float>("pT");
					m_kshort_pTCalc = sqrt(pow(m_kshort_px ,2) + pow(m_kshort_py ,2));
					m_kshort_p  = sqrt(pow(m_kshort_px,2) + pow(m_kshort_py,2) + pow(m_kshort_pz,2));
					m_kshort_theta=acos(m_kshort_pz/m_kshort_p);
					m_kshort_eta=log(1.0/(tan(m_kshort_theta/2)));
					rdvDotPt=vertex_ptr->x()*m_kshort_px+vertex_ptr->y()*m_kshort_py;
					m_kshort_alpha = acos(rdvDotPt/(m_kshort_rDV*m_kshort_pTCalc));
					m_kshort_x = vertex_ptr->x();
					m_kshort_y = vertex_ptr->y();
					m_kshort_z = vertex_ptr->z();

					// Related to both Ks and primary vertices
					m_z_sv_pv = abs( m_kshort_z - ( (m_kshort_rDV * m_kshort_pz) / (m_kshort_pTCalc)) - m_primary_vertex_z );



					//init the truth vars
					m_truth_piplus_pdgid = -999;
					m_truth_piminus_pdgid = -999;
					m_truth_kshort_pdgid = -999;

					const xAOD::EventInfo* ei_ptr = 0;
    					CHECK(evtStore()->retrieve(ei_ptr,"EventInfo"));
					if(ei_ptr->eventType(xAOD::EventInfo::IS_SIMULATION))
  					{

						// For MC Truth				
						// Get truth particle from pi+ and pi-
						const xAOD::TruthParticle *piplus_truth = xAOD::TruthHelpers::getTruthParticle(*piplus_track);
						const xAOD::TruthParticle *piminus_truth = xAOD::TruthHelpers::getTruthParticle(*piminus_track);
	
						if(piplus_truth && piminus_truth)
						{
							m_truth_piplus_pdgid = piplus_truth->pdgId();	
							m_truth_piminus_pdgid = piminus_truth->pdgId();
	
							if (piplus_truth->nParents() == 1 && piminus_truth->nParents() == 1) 
							{
								if (piplus_truth->parent(0) == piminus_truth->parent(0)) 
								{
									const xAOD::TruthParticle *potential_truth_ks = piplus_truth->parent(0);
									m_truth_kshort_pdgid = potential_truth_ks->pdgId();
									if(m_truth_kshort_pdgid == m_KS_PDG)
									{
										m_total_rec_ks_hits++;	
									}
	
								}	
							}
						}
					}

					// Elements of the 3x3 error matrix related to Ks vertices
                    			AmgSymMatrix(3) covariance_matrix = vertex_ptr->covariancePosition();
                    			m_covariance00 = covariance_matrix(0,0);
                    			m_covariance01 = covariance_matrix(0,1);
                    			m_covariance02 = covariance_matrix(0,2);
                    			m_covariance10 = covariance_matrix(1,0);
                    			m_covariance11 = covariance_matrix(1,1);
                    			m_covariance12 = covariance_matrix(1,2);
                    			m_covariance20 = covariance_matrix(2,0);
                    			m_covariance21 = covariance_matrix(2,1);
                    			m_covariance22 = covariance_matrix(2,2);


					this->m_kshort_tree->Fill();

				}
			}
		}
	}	
	
	return sc;	
}



