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
  	/*
	StatusCode sc = service("StoreGateSvc", m_storeGate);
	if (sc.isFailure()) 
	{
		msg(MSG::ERROR) << "Unable to retrieve pointer to StoreGateSvc" << endreq;
		return sc;
	}

	sc = m_trigDec.retrieve();
	if (sc.isFailure())
	{
		msg(MSG::ERROR) << "Can't get handle on TrigDecisionTool" << endreq;
	}

	sc = service("THistSvc", m_thistSvc);
	if (sc.isFailure()) {
		msg(MSG::ERROR) << "Unable to retrieve pointer to THistSvc" << endreq;
		return sc;
	}
	
	// Histogramming...
	ServiceHandle<ITHistSvc> histSvc("THistSvc",name());
	CHECK( histSvc.retrieve() );

	if (sc.isFailure()) 
	{
		msg(MSG::ERROR) << "Unable to register histogramming service" << endreq;
		return sc;
	}
	*/

	//// ATH_MSG_INFO ("Initializing " << name() << "...");

	// Create a tree
	//m_Ks_tree = new TTree("Kshort_DDL", "");

	//m_vxtree->SetTree(m_Ks_tree);
	// Register the tree
	//std::string hist_path = "/" + m_hist_name + "/";
	//CHECK(m_histSvc->regTree(hist_path + m_Ks_tree->GetName(), m_Ks_tree));


	TTree *kshort_tree = new TTree("KshortTTree","Kshort Tree!");


	kshort_tree->Branch("piplus_pt",&piplus_pt,"piplus_pt/f");
	kshort_tree->Branch("piplus_pt",&piplus_p,"piplus_pt/f");
	kshort_tree->Branch("piplus_pt",&piplus_px,"piplus_pt/f");
	kshort_tree->Branch("piplus_pt",&piplus_py,"piplus_pt/f");
	kshort_tree->Branch("piplus_pt",&piplus_pz,"piplus_pt/f");
	kshort_tree->Branch("piplus_pt",&piplus_e,"piplus_pt/f");
	kshort_tree->Branch("piplus_pt",&piplus_z0,"piplus_pt/f");
	kshort_tree->Branch("piplus_pt",&piplus_d0,"piplus_pt/f");
	kshort_tree->Branch("piplus_pt",&piplus_eta,"piplus_pt/f");
	
	kshort_tree->Branch("piplus_pt",&piminus_pt,"piplus_pt/f");
	kshort_tree->Branch("piplus_pt",&piminus_p,"piplus_pt/f");
	kshort_tree->Branch("piplus_pt",&piminus_px,"piplus_pt/f");
	kshort_tree->Branch("piplus_pt",&piminus_py,"piplus_pt/f");
	kshort_tree->Branch("piplus_pt",&piminus_pz,"piplus_pt/f");
	kshort_tree->Branch("piplus_pt",&piminus_e,"piplus_pt/f");
	kshort_tree->Branch("piplus_pt",&piminus_z0,"piplus_pt/f");
	kshort_tree->Branch("piplus_pt",&piminus_d0,"piplus_pt/f");
	kshort_tree->Branch("piplus_pt",&piminus_eta,"piplus_pt/f");


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
      	return sc;     
  //// ATH_MSG_DEBUG ("Executing " << name() << "...");
}


// Function for charged pion mass (equating float to some number inside a condition is not recommended, hence this function was made) 
bool DDL::Kshort_DDL::isPi(float piMass)
{
        float actual_pi_mass = 139.57061; // MeV (taken from recent PDG, http://pdglive.lbl.gov/Particle.action?init=0&node=S008)
        float margin_of_error = 0.00024; // MeV (taken from recent PDG, http://pdglive.lbl.gov/Particle.action?init=0&node=S008)

        return ( (piMass >= (actual_pi_mass - margin_of_error)) &&  (piMass <= (actual_pi_mass + margin_of_error)) );
}


StatusCode DDL::Kshort_DDL::finding_right_ks() 
{

	StatusCode sc = StatusCode::SUCCESS; 

//----------------------------------------------------------------------------


	const xAOD::VertexContainer* vertices = nullptr;
	CHECK(evtStore()->retrieve(vertices, "VrtSecIncludive_SecondaryVertices")); // Is VrtSecIncludive_SecondaryVertices OK or any name can replace it ?
  
	// "TrackParticle": pointer to a given track that was used in vertex reconstruction
	const xAOD::TrackParticle* piplus_track = nullptr;
	const xAOD::TrackParticle* piminus_track = nullptr;

	// Vertex container for Kshorts
	//xAOD::VertexContainer* kshortVertices = new VertexContainer();
  
  
	// Looping over the vertices to find Kshort vertices which decay to pi+ and pi-
	// Defining the iterator in a short and efficient way
	for(xAOD::Vertex* vertex_ptr : *vertices)
	{
        // Vertices with only two tracks 
		if(vertex_ptr->nTrackParticles()==2)
		{
			piplus_track = vertex_ptr->trackParticle(0);
			piminus_track = vertex_ptr->trackParticle(1);
			if(piplus_track->charge()<piminus_track->charge()) //in case track 0 is piMinus and track 1 is piPlus
			{
				std::swap(piminus_track,piplus_track);
			}
			// Checking whether the masses are equal to the charged pion mass
			if( isPi(piplus_track->m()) && isPi(piminus_track->m()) )
			{
				// Checking if they have charges of -1 and +1
				if(piplus_track->charge()+piminus_track->charge()==0 && piplus_track->charge()==1)
                		{      
					// How to save the mass of such vertices and other parameters ?
					// The idea is to save all these variables in form of a Tree.	
					// Mass of the Kshort vertices after confirming that they are indeed Kshort vertices
					// Also a storage for the Kshort vertices
					//KshortVertices->push_back(vertex_ptr); //not sure we need it
   				

					
					// Saving individual track info in ttree
					
					piplus_pt 	= 	piplus_track->pt();
					piplus_p  	= 	piplus_track->p4().P();
					piplus_px	= 	piplus_track->p4().Px();
					piplus_py	= 	piplus_track->p4().Py();
					piplus_pz	= 	piplus_track->p4().Pz();
					piplus_e	= 	piplus_track->e();
					piplus_z0	=	piplus_track->z0();
					piplus_d0	=	piplus_track->d0();
					piplus_eta	=	piplus_track->eta();

					piminus_pt 	= 	piminus_track->pt();
					piminus_p  	= 	piminus_track->p4().P();
					piminus_px	= 	piminus_track->p4().Px();
					piminus_py	= 	piminus_track->p4().Py();
					piminus_pz	= 	piminus_track->p4().Pz();
					piminus_e	= 	piminus_track->e();
					piminus_z0	=	piminus_track->z0();
					piminus_d0	=	piminus_track->d0();
					piminus_eta	=	piminus_track->eta();

					kshort_tree->Fill();
				
				}
			}			
		}	
	}

  // Track container for all tracks 
  /*
  const xAOD::TrackParticleContainer* recoTracks = 0;
  CHECK( evtStore()->retrieve( recoTracks, "InDetTrackParticles" ) ); 
  
  // pi+ tracks
  xAOD::TrackParticleContainer* piplus_tracks = new TrackParticleContainer();
  xAOD::TrackParticleContainer* piminus_tracks = new TrackParticleContainer();
  //CHECK( evtStore()->retrieve( recoTracks_piplus, "InDetTrackParticles" ) );
  for(xAOD:TrackParticle* track : recoTracks)
  {
	if(track->m()==0.1396))
	{
		if(track->charge()==1)
		{
			piplus_tracks.push_back(track);
		}
		else if(track->charge()==-1)
		{
			pminus_tracks.push_back(track);
		}
		
	}
  } 
  */

  // How to combine 'recoTracks_piplus' and 'recoTracks_piminus'

 //----------------------------------------------------------------------------

  // Vertex container 

 //----------------------------------------------------------------------------

  // Truth particle container
	/*
	const xAOD::TruthParticleContainer* truthParticles;

	isMC = false;
    sc=m_storeGate->retrieve(truthParticles, "TruthParticles");
    if(!sc.isFailure())
	{  
      isMC = true;
    }
	*/
	/*
	// Number of all reconstructed tracks
	m_nTracks->Fill(recoTracks->size());
	// Number of reconstructed pi+ tracks
  
	// Number of reconstructed pi- tracks

	// How to get to two vertices of pi+ and pi- using VertexContainer
	double counterDV = 0; 
	float vtxDist = 0;
	float vtxMass = -999;
	double vtxX = -999; 
	double vtxY = -999; 
	double vtxZ = -999;
	double vtxType = -999; 
	double vtxNbTracks = -999; // number of tracks that make the vertex 
	double vtxNbTrackParticles = -999; 
  
	for(xAOD::VertexContainer::const_iterator rVtx_itr = vertexContainer->begin();rVtx_itr != vertexContainer->end(); ++rVtx_itr)
	{//Loop on RECO vertex

      const xAOD::Vertex *dv = (*rVtx_itr);
      vtxDist = sqrt(((*rVtx_itr)->x()*(*rVtx_itr)->x())+((*rVtx_itr)->y()*(*rVtx_itr)->y()));
      vtxMass = (*rVtx_itr)->auxdataConst<float>("vtx_mass");//Works only on RRV ?
      // how to use 'vtxMass' for our purpose?  
      vtxX = (*rVtx_itr)->x();
      vtxY = (*rVtx_itr)->y();
      vtxZ = (*rVtx_itr)->z();
      vtxType = (*rVtx_itr)->vertexType();
      vtxNbTrackParticles = (*rVtx_itr)->nTrackParticles();
         
       if(((*rVtx_itr)->auxdataConst<float>("vtx_mass")/1000.0)>1){//If DVmass>1GeV
                h_vtxMass_cut1G->Fill((*rVtx_itr)->auxdataConst<float>("vtx_mass")/1000.0);
                h_vtxdR_cut1G->Fill(delta_R((*mu_itr)->eta(),(*mu_itr)->phi(),dvTrk->eta(),dvTrk->phi()));
             }//End of DVmass1G     
  
      if((*rVtx_itr)->nTrackParticles()>2){continue;}//If vtx has more than 2 tracks, stop.

    }
	*/
	/*
	double counterKs = 0;
	for(size_t tk = 0; tk < (*rVtx_itr)->nTrackParticles(); ++tk)
	{//loop over all the tracks from the recoVertex
	}	
	*/
	return sc;
}
