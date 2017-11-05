#ifndef DDLSTUDIES_KSHORT_DDL_H
#define DDLSTUDIES_KSHORT_DDL_H 

#include "AthenaBaseComps/AthAnalysisAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/ITHistSvc.h"
#include "StoreGate/StoreGateSvc.h" 

// DDL
#include "DDLBase/IDiLepDVCuts.h"
#include "DDLBase/IDVCuts.h"

#include "TrigDecisionTool/TrigDecisionTool.h"
#include "JetCalibTools/IJetCalibrationTool.h"

#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/Vertex.h"
#include "TTree.h"

// forward declaration
class TH1I;
class TH1F;
class TH2F;
class TH1D;
class TEfficiency;



namespace DDL {

    class Kshort_DDL: public ::AthAnalysisAlgorithm 
    {
 
        public: 
            Kshort_DDL( const std::string& name, ISvcLocator* pSvcLocator );
            virtual ~Kshort_DDL() = default; 

            StatusCode  initialize() override;
            StatusCode  execute() override;
            StatusCode  finalize() override;

            // for good Kshort ?
            virtual bool goodKs () const;

        private: 

   };

} 

#endif // DDLSTUDIES_KSHORT_DDL_H  
