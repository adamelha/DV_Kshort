// Copy this file to DDLStudies/src


// DDLStudies includes
#include "Kshort_DDL.h"

// framework includes
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ITHistSvc.h"
#include "StoreGate/StoreGate.h"
#include "StoreGate/StoreGateSvc.h"
#include "StoreGate/DataHandle.h"

				
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

#include "StoreGate/StoreGate.h"
#include "StoreGate/StoreGateSvc.h"
#include "StoreGate/DataHandle.h"

//Event includes
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "TrigDecisionTool/ChainGroup.h"
#include "TrigDecisionTool/FeatureContainer.h"
#include "TrigDecisionTool/Feature.h"


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



DDL::Kshort_DDL::Kshort_DDL( const std::string& name, ISvcLocator* pSvcLocator ) : 
   AthAnalysisAlgorithm( name, pSvcLocator )//,m_nTracks(0)   
{

  //declareProperty("HistSvcName", m_hist_name = "Kshort_DDL" ); 

}

StatusCode DDL::Kshort_DDL::initialize() 
{ 
	std::cout << "initialize!!!" << std::endl;

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


	m_kshort_tree = new TTree("KshortTTree", "Kshort Tree!");


	CHECK(histSvc->regTree("/Kshort_DDL/m_kshort_tree", m_kshort_tree));
	//CHECK( histSvc->regHist("/HNL/m_nTracks",m_nTracks) );


	m_kshort_tree->Branch("piplus_pt", &piplus_pt, "piplus_pt/f");
	m_kshort_tree->Branch("piplus_p", &piplus_p, "piplus_p/f");
	m_kshort_tree->Branch("piplus_px", &piplus_px, "piplus_px/f");
	m_kshort_tree->Branch("piplus_py", &piplus_py, "piplus_py/f");
	m_kshort_tree->Branch("piplus_pz", &piplus_pz, "piplus_pz/f");
	m_kshort_tree->Branch("piplus_e", &piplus_e, "piplus_e/f");
	m_kshort_tree->Branch("piplus_z0", &piplus_z0, "piplus_z0/f");
	m_kshort_tree->Branch("piplus_d0", &piplus_d0, "piplus_d0/f");
	m_kshort_tree->Branch("piplus_eta", &piplus_eta, "piplus_eta/f");

	m_kshort_tree->Branch("piminus_pt", &piminus_pt, "piminus_pt/f");
	m_kshort_tree->Branch("piminus_p", &piminus_p, "piminus_p/f");
	m_kshort_tree->Branch("piminus_px", &piminus_px, "piminus_px/f");
	m_kshort_tree->Branch("piminus_py", &piminus_py, "piminus_py/f");
	m_kshort_tree->Branch("piminus_pz", &piminus_pz, "piminus_pz/f");
	m_kshort_tree->Branch("piminus_e", &piminus_e, "piminus_e/f");
	m_kshort_tree->Branch("piminus_z0", &piminus_z0, "piminus_z0/f");
	m_kshort_tree->Branch("piminus_d0", &piminus_d0, "piminus_d0/f");
	m_kshort_tree->Branch("piminus_eta", &piminus_eta, "piminus_eta/f");
	return StatusCode::SUCCESS;
}

StatusCode DDL::Kshort_DDL::finalize() 
{

  //// ATH_MSG_INFO ("Finalizing " << name() << "...");

  return StatusCode::SUCCESS;
}

StatusCode DDL::Kshort_DDL::execute() 
{
	this->event_counter++;
	std::cout << "Event Number : " << this->event_counter << std::endl;
	StatusCode sc = StatusCode::SUCCESS;
	/*
	sc = fillKs();
    	if ( sc.isFailure() ) 
	{
       		msg(MSG::ERROR) << "Ks filling failed" << endreq;
       		return StatusCode::SUCCESS;
    	}
	*/  
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
bool DDL::Kshort_DDL::isPi(float piMass)
{
        float actual_pi_mass = 139.57061; // MeV (taken from recent PDG, http://pdglive.lbl.gov/Particle.action?init=0&node=S008)
        //float margin_of_error = 0.00024; // MeV (taken from recent PDG, http://pdglive.lbl.gov/Particle.action?init=0&node=S008)
		float margin_of_error = 0.02; // MeV (taken from recent PDG, http://pdglive.lbl.gov/Particle.action?init=0&node=S008)
        return ( (piMass >= (actual_pi_mass - margin_of_error)) &&  (piMass <= (actual_pi_mass + margin_of_error)) );
}


StatusCode DDL::Kshort_DDL::finding_right_ks()
{

	StatusCode sc = StatusCode::SUCCESS;

	//----------------------------------------------------------------------------


	const xAOD::VertexContainer* vertices = nullptr;
	CHECK(evtStore()->retrieve(vertices, "PrimaryVertices"));
	//CHECK(evtStore()->retrieve(vertices, "VrtSecIncludive_SecondaryVertices")); // Is VrtSecIncludive_SecondaryVertices OK or any name can replace it ?

	// "TrackParticle": pointer to a given track that was used in vertex reconstruction
	const xAOD::TrackParticle* piplus_track = nullptr;
	const xAOD::TrackParticle* piminus_track = nullptr;

	// Vertex container for Kshorts
	//xAOD::VertexContainer* kshortVertices = new VertexContainer();


	// Looping over the vertices to find Kshort vertices which decay to pi+ and pi-
	// Defining the iterator in a short and efficient way
	for (xAOD::Vertex* vertex_ptr : *vertices)
	{
		if (vertex_ptr->nTrackParticles() >= 1) {
			cout << "Partical 0 mass: " << vertex_ptr->trackParticle(0)->m() << std::endl;
		}
		cout << "There are " << vertex_ptr->nTrackParticles() << " in this vertex" << std::endl;
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
				cout << "Both are Pis! piplus_track->charge() = " << piplus_track->charge() << " AND piminus_track->charge() = " << piminus_track->charge() << std::endl;
				// Checking if they have charges of -1 and +1
				if (piplus_track->charge() + piminus_track->charge() == 0 && piplus_track->charge() == 1)
				{
					// How to save the mass of such vertices and other parameters ?
					// The idea is to save all these variables in form of a Tree.	
					// Mass of the Kshort vertices after confirming that they are indeed Kshort vertices
					// Also a storage for the Kshort vertices
					//KshortVertices->push_back(vertex_ptr); //not sure we need it



					// Saving individual track info in ttree

					piplus_pt = piplus_track->pt();
					piplus_p = piplus_track->p4().P();
					piplus_px = piplus_track->p4().Px();
					piplus_py = piplus_track->p4().Py();
					piplus_pz = piplus_track->p4().Pz();
					piplus_e = piplus_track->e();
					piplus_z0 = piplus_track->z0();
					piplus_d0 = piplus_track->d0();
					piplus_eta = piplus_track->eta();

					piminus_pt = piminus_track->pt();
					piminus_p = piminus_track->p4().P();
					piminus_px = piminus_track->p4().Px();
					piminus_py = piminus_track->p4().Py();
					piminus_pz = piminus_track->p4().Pz();
					piminus_e = piminus_track->e();
					piminus_z0 = piminus_track->z0();
					piminus_d0 = piminus_track->d0();
					piminus_eta = piminus_track->eta();

					std::cout << "Updating the tree" << std::endl;
					this->m_kshort_tree->Fill();
					std::cout << "Updated" << std::endl;
				}
			}
		}

	}
	return sc;
}
