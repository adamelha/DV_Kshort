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

//#include "xAODTruth/TruthParticleContainer.h"
//#include "xAODTruth/TruthParticle.h"
//#include "xAODTruth/TruthVertex.h"
//#include "xAODTruth/TruthVertexContainer.h"
//#include "xAODTruth/xAODTruthHelpers.h"


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
	std::cout << "initialize" << std::endl;

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


	m_kshort_tree = new TTree("KsTTree", "Ks Tree!");


	CHECK(histSvc->regTree("/Ks/m_kshort_tree", m_kshort_tree));

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
	m_kshort_tree->Branch("kshort_mass", &m_kshort_mass, "kshort_mass/D");
	m_kshort_tree->Branch("kshort_rDV",  &m_kshort_rDV,  "kshort_rDV/D" );
	m_kshort_tree->Branch("kshort_theta",&m_kshort_theta,  "kshort_theta/D" );
	m_kshort_tree->Branch("kshort_eta",  &m_kshort_eta,  "kshort_eta/D" );
	m_kshort_tree->Branch("kshort_e",    &m_kshort_e,    "kshort_e/D"   );
	m_kshort_tree->Branch("kshort_pt",   &m_kshort_pt,   "kshort_pt/D"  );
	m_kshort_tree->Branch("kshort_p",    &m_kshort_p,    "kshort_p/D"   );
	m_kshort_tree->Branch("kshort_px",   &m_kshort_px,   "kshort_px/D"  );
	m_kshort_tree->Branch("kshort_py",   &m_kshort_py,   "kshort_py/D"  );
	m_kshort_tree->Branch("kshort_pz",   &m_kshort_pz,   "kshort_pz/D"  );
	m_kshort_tree->Branch("kshort_pTCalc",&m_kshort_pTCalc,"kshort_pTCalc/D");
	m_kshort_tree->Branch("kshort_invMass",&m_kshort_invMass,"kshort_invMass/D");
	m_kshort_tree->Branch("kshort_alpha",&m_kshort_alpha, "kshort_alpha/D");

	m_kshort_tree->Branch("primary_vertex_x",&m_primary_vertex_x,"primary_vertex_x/D");
	m_kshort_tree->Branch("primary_vertex_y",&m_primary_vertex_y,"primary_vertex_y/D");
	m_kshort_tree->Branch("primary_vertex_z",&m_primary_vertex_z,"primary_vertex_z/D");


	return StatusCode::SUCCESS;
}

StatusCode DDL::Ks::finalize() 
{

  //// ATH_MSG_INFO ("Finalizing " << name() << "...");

  return StatusCode::SUCCESS;
}

StatusCode DDL::Ks::execute() 
{
	this->event_counter++;
	//std::cout << "Event Number : " << this->event_counter << std::endl;
	StatusCode sc = StatusCode::SUCCESS;  
	sc = finding_right_ks();
	if (sc.isFailure())
	{
		msg(MSG::ERROR) << "Finding right ks failed" << endreq;
		return StatusCode::SUCCESS;
	}
	return sc;
  //// ATH_MSG_DEBUG ("Executing " << name() << "...");
}


// Function for charged pion mass (equating float to some number inside a condition is not recommended, hence this function was made) 
bool DDL::Ks::isPi(float piMass)
{
        float actual_pi_mass = 139.57061; // MeV (taken from recent PDG, http://pdglive.lbl.gov/Particle.action?init=0&node=S008)
        //float margin_of_error = 0.00024; // MeV (taken from recent PDG, http://pdglive.lbl.gov/Particle.action?init=0&node=S008)
		float margin_of_error = 0.02; // MeV (taken from recent PDG, http://pdglive.lbl.gov/Particle.action?init=0&node=S008)
        return ( (piMass >= (actual_pi_mass - margin_of_error)) &&  (piMass <= (actual_pi_mass + margin_of_error)) );
}


StatusCode DDL::Ks::finding_right_ks()
{

	StatusCode sc = StatusCode::SUCCESS;
	//----------------------------------------------------------------------------


	const xAOD::VertexContainer* secondary_vertices = nullptr;
	const xAOD::VertexContainer* primary_vertices = nullptr;

	CHECK(evtStore()->retrieve(primary_vertices, "PrimaryVertices"));
	CHECK(evtStore()->retrieve(secondary_vertices , "VrtSecInclusive_SecondaryVertices")); // Is VrtSecInclusive_SecondaryVertices OK or any name can replace it ?

	std::cout << primary_vertices->size() << " : " << secondary_vertices->size() << std::endl;


	// "TrackParticle": pointer to a given track that was used in vertex reconstruction
	const xAOD::TrackParticle* piplus_track = nullptr;
	const xAOD::TrackParticle* piminus_track = nullptr;

	TLorentzVector p4_sum;
	Double_t r;	

	// Looping over the vertices to find Kshort vertices which decay to pi+ and pi-
	// Defining the iterator in a short and efficient way
	for (xAOD::Vertex* vertex_ptr : *secondary_vertices )
	{
		// Vertices with only two tracks 
		if (vertex_ptr->nTrackParticles() == 2)
		{
			
			piplus_track = vertex_ptr->trackParticle(0);
			piminus_track = vertex_ptr->trackParticle(1);

			if (piplus_track->charge() < piminus_track->charge()) //in case track 0 is piMinus and track 1 is piPlus
			{
				std::swap(piminus_track, piplus_track);
			}
			// Checking whether the masses are equal to the charged pion mass
			if (isPi(piplus_track->m()) && isPi(piminus_track->m()))
			{
				// Checking if they have charges of -1 and +1
				if (piplus_track->charge() + piminus_track->charge() == 0 && piplus_track->charge() == 1)
				{
					// Saving individual track info in ttree
					//std::cout << "Pi+ and Pi found" << std::endl;					
					//std::cout << "Pi+:pt = " << piplus_track->pt() << ", Pi-:pt = "  << piminus_track->pt() << std::endl;
					//std::cout << "Pi+:z0 = " << piplus_track->z0() << ", Pi-:z0 = "  << piminus_track->z0() << std::endl;
					//std::cout << "Pi+:d0 = " << piplus_track->d0() << ", Pi-:d0 = "  << piminus_track->d0() << std::endl;
					m_piplus_pt   = piplus_track->pt();
					m_piplus_p    = piplus_track->p4().P();
					m_piplus_px   = piplus_track->p4().Px();
					m_piplus_py   = piplus_track->p4().Py();
					m_piplus_pz   = piplus_track->p4().Pz();
					m_piplus_e    = piplus_track->e();
					m_piplus_z0   = piplus_track->z0();
					m_piplus_d0   = piplus_track->d0();
					m_piplus_eta  = piplus_track->eta();
				

					m_piminus_pt  = piminus_track->pt();
					m_piminus_p   = piminus_track->p4().P();
					m_piminus_px  = piminus_track->p4().Px();
					m_piminus_py  = piminus_track->p4().Py();
					m_piminus_pz  = piminus_track->p4().Pz();
					m_piminus_e   = piminus_track->e();
					m_piminus_z0  = piminus_track->z0();
					m_piminus_d0  = piminus_track->d0();
					m_piminus_eta = piminus_track->eta();
					
					p4_sum = piplus_track->p4() + piminus_track->p4();
					m_kshort_invMass = p4_sum.Mag();

					m_kshort_mass = vertex_ptr->auxdataConst<float>("vtx_mass");
					m_kshort_rDV  = sqrt (pow(vertex_ptr->x(),2) + pow(vertex_ptr->y(),2)); 
					m_kshort_e  = sqrt(pow(m_kshort_mass ,2)+pow(m_kshort_p ,2)); 
					m_kshort_px = vertex_ptr->auxdataConst<float>("vtx_px");
					m_kshort_py = vertex_ptr->auxdataConst<float>("vtx_py");
					m_kshort_pz = vertex_ptr->auxdataConst<float>("vtx_pz");

					m_kshort_pt = vertex_ptr->auxdataConst<float>("pT");
					m_kshort_pTCalc = sqrt (pow(m_kshort_px ,2) + pow(m_kshort_py ,2));

					m_kshort_p  = sqrt(pow(m_kshort_px,2) + pow(m_kshort_py,2) + pow(m_kshort_pz,2));
					
					m_kshort_theta=acos(m_kshort_pz/m_kshort_p);
					m_kshort_eta=log(1.0/(tan(m_kshort_theta/2)));

					//r = sqrt (pow(vertex_ptr->x(),2) + pow(vertex_ptr->y(),2) + pow(vertex_ptr->z(),2));
					m_kshort_alpha = acos(abs(vertex_ptr->x()*m_kshort_px+vertex_ptr->y()*m_kshort_py)/(m_kshort_rDV*m_kshort_pTCalc));

					this->m_kshort_tree->Fill();
				}
			}
		}
	}
	for (xAOD::Vertex* vertex_ptr : *primary_vertices )
	{
		m_primary_vertex_x = vertex_ptr->x();
		m_primary_vertex_y = vertex_ptr->y();
		m_primary_vertex_z = vertex_ptr->z();
		this->m_kshort_tree->Fill();
	}

	return sc;	
}
