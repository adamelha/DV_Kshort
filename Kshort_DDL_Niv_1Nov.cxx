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
#include "xAODTracking/TrackMeasurementValidationContainer.h"
#include "xAODTracking/TrackStateValidationContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/TrackParticleAuxContainer.h"
#include "xAODTracking/VertexContainer.h"

#include "xAODTracking/Vertex.h"
#include "xAODTracking/VertexAuxContainer.h"

#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODTruth/TruthVertexContainer.h"
#include "xAODTruth/xAODTruthHelpers.h"

#include "StoreGate/StoreGate.h"
#include "StoreGate/StoreGateSvc.h"
#include "StoreGate/DataHandle.h"

//Event includes
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "TrigDecisionTool/ChainGroup.h"
#include "TrigDecisionTool/FeatureContainer.h"
#include "TrigDecisionTool/Feature.h"


//ROOT
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
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
   AthAnalysisAlgorithm( name, pSvcLocator ),
   m_nTracks(0)   
{

  declareProperty( "Property", m_nProperty ); //example property declaration

}


// Kshort_DDL::~Kshort_DDL() {}


StatusCode DDL::Kshort_DDL::initialize() { 

  using namespace std;
  
  StatusCode sc = service("StoreGateSvc", m_storeGate);
  if (sc.isFailure()) {
     msg(MSG::ERROR) << "Unable to retrieve pointer to StoreGateSvc" << endreq;
     return sc;
  }

  sc = m_trigDec.retrieve();
  if (sc.isFailure()){
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

  if (sc.isFailure()) {
     msg(MSG::ERROR) << "Unable to register histogramming service" << endreq;
     return sc;
  }

  // Creation of new histograms : pi+, pi-, Kshort histograms
  h_Name = new TH1I("h_Name", "; Title",n_bins,n_min,n_max);
  h_Name2 = new TH2F("h_Name2", "; Title1; Title2", n_bins1,n_min1,n_max1,n_bins2,n_min2,n_max2);
 
  // Efficiency histograms
  h_Name_Eff = new TH1F("h_Name_Eff", " Quantity1 vs. Quantity2 ", bins, q_min, q_max);

  // Truth particle histograms 
 
  // Vertex Histograms
 
  // Cutflow histograms      

  //// ATH_MSG_INFO ("Initializing " << name() << "...");

  return StatusCode::SUCCESS;
}

StatusCode DDL::Kshort_DDL::finalize() {

  //// ATH_MSG_INFO ("Finalizing " << name() << "...");

  return StatusCode::SUCCESS;
}

StatusCode DDL::Kshort_DDL::execute() {
      StatusCode sc = StatusCode::SUCCESS;
      sc = fillKs();
      if ( sc.isFailure() ) {
         msg(MSG::ERROR) << "Ks filling failed" << endreq;
         return StatusCode::SUCCESS;
      }  
      return sc;     
  //// ATH_MSG_DEBUG ("Executing " << name() << "...");
}


bool DDL::Kshort_DDL::isPi(float particalMass)
{
	float theoretical_pi_mass = 139.5706; // MeV
	float margin_of_error = 0.00024
	
	return (particalMass >= theoretical_pi_mass - margin_of_error) && (particalMass <= theoretical_pi_mass + margin_of_error);
}

StatusCode DDL::Kshort_DDL::fillKs() {

  using namespace std;
  int counterEvent = 0;
  StatusCode sc = StatusCode::SUCCESS; 

//----------------------------------------------------------------------------



  const xAOD::VertexContainer* vertices = nullptr;
  CHECK(evtStore()->retrieve(vertices, "VrtSecIncludive_SecondaryVertices");
  
  xAOD::TrackParticle* first_track = nullptr;
  xAOD::TrackParticle* second_track = nullptr;

  //vertex container of kshorts
  xAOD::VertexContainer* kshortVertices = new VertexContainer();
  
  
  //we loop over the vertices to find kshort vertices that decay to pi+ and pi-
  for(xAOD::Vertex* vertex_ptr : vertices)
  {
	if(vertex_ptr->nTrackParticle()==2)
	{
		first_track = vertex_ptr->trackParticle(0);
		second_track = vertex_ptr->trackParticle(1);
		//lets check that the mass are equal and equal to pi mass
		if( isPi(first_track->m()) && isPi(second_track->m()))
		{
			//lets check if they have charges of -1 and +1.
			if(first_track->charge()+second_track->charge()==0 && std::abs(first_track->charge())==1)
			{
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
  const xAOD::VertexContainer* vertexContainer = 0;
  StatusCode sc=m_storeGate->retrieve(vertexContainer, "VrtSecInclusive_SecondaryVertices");//If RRV sample
  if(sc==StatusCode::FAILURE)return;

 //----------------------------------------------------------------------------

  // Truth particle container
  const xAOD::TruthParticleContainer* truthParticles;

  isMC = false;
    sc=m_storeGate->retrieve(truthParticles, "TruthParticles");
    if(!sc.isFailure()){  
      isMC = true;
    }

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
  
  for(xAOD::VertexContainer::const_iterator rVtx_itr = vertexContainer->begin();rVtx_itr != vertexContainer->end(); ++rVtx_itr){//Loop on RECO vertex
     
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

   double counterKs = 0;
   for(size_t tk = 0; tk < (*rVtx_itr)->nTrackParticles(); ++tk){//loop over all the tracks from the recoVertex
     

   }

}

// Do we really need this? 
// Looping on tracks
void DDL::Kshort_DDL::fillTrack(){

  // Get to the track container 
  const xAOD::TrackParticleContainer* recoTracks = 0;
  
  StatusCode sc=m_storeGate->retrieve(recoTracks, "InDetTrackParticles");

  // obtain pi+ tracks, make one loop ; then obtain pi- tracks, make anothe loop in the pi+ loop
  for(xAOD::TrackParticleContainer::const_iterator recoTrk_itr = recoTracks->begin();recoTrk_itr != recoTracks->end(); recoTrk_itr++){ 

    for(xAOD::TrackParticleContainer::const_iterator recoTrk_itr2 = recoTrk_itr+1;recoTrk_itr2 != recoTrk_itr->end(); recoTrk_itr2++){

    } 
    (*recoTrk_itr)->TrackParticle()->charge();
    (*recoTrk_itr)->TrackParticle()->m();
    (*recoTrk_itr)->TrackParticle()->m2();
    (*recoTrk_itr)->TrackParticle()->px();
    (*recoTrk_itr)->TrackParticle()->py();
    (*recoTrk_itr)->TrackParticle()->pz();
    (*recoTrk_itr)->TrackParticle()->p();
    (*recoTrk_itr)->TrackParticle()->p2();
    (*recoTrk_itr)->TrackParticle()->eta();
    (*recoTrk_itr)->TrackParticle()->e();
    (*recoTrk_itr)->TrackParticle()->phi();
    (*recoTrk_itr)->TrackParticle()->pt();
    (*recoTrk_itr)->TrackParticle()->reconstructedVertex();          

  }



}
