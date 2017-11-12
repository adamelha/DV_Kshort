##############################################################################################
#                                                                                            #
#     Job option file for a full athena environment for Kshort reconstruction in DDL package #
#                                                                                            #
##############################################################################################

# Enable file reading
import AthenaPoolCnvSvc.ReadAthenaPool
from AthenaCommon.AthenaCommonFlags import athenaCommonFlags as acf


# Get input file(s) for running on the local environment. If want to run on GRID -- 
from glob import glob
acf.FilesInput =["/afs/cern.ch/user/n/nkonstan/event_files/DAOD_RPVLL.10626802._032446.pool.root.1"]	

#import our algorithm
#from DDL.DDLStudiesConf import DDL__Kshort_DDL

# Config rec
from RecExConfig.RecFlags import rec
rec.OutputLevel = WARNING
rec.doCBNT = False
rec.doFileMetaData.set_Value_and_Lock( False )
rec.doWriteESD.set_Value_and_Lock(False)
rec.doWriteAOD.set_Value_and_Lock(False)
rec.doWriteTAG.set_Value_and_Lock(False)
rec.doPerfMon=False

rec.doHist      = False
rec.doTrigger   = True
rec.doWritexAOD = False

from RecExConfig.RecAlgsFlags import recAlgs
recAlgs.doTrigger = True
from TriggerJobOpts.TriggerFlags import TriggerFlags
TriggerFlags.doTriggerConfigOnly = True

include("RecExCommon/RecExCommon_topOptions.py")

# Set events
theApp.EvtMax = vars().get("maxEvents", -1)
svcMgr.EventSelector.SkipEvents = vars().get("skipEvents", 0)

# Setup trigger tool
ToolSvc += CfgMgr.Trig__TrigDecisionTool("TrigDecisionTool",
                                         OutputLevel = WARNING,
                                         TrigDecisionKey = "xTrigDecision")

# Setup THistSvc
svc_name = "Kshort_DDL"
svcMgr += CfgMgr.THistSvc()
svcMgr.THistSvc.Output += ["%s DATAFILE='Kshort_DDL_output.root' OPT='RECREATE'" % (svc_name)]

# Required for Kshort reconstruction?
trackPtCut = 10000.

# Setup PionCuts


# Setup KshortCuts


# Setup DVCuts
# DVCuts necessary for Kshort reconstrauction?
ToolSvc += CfgMgr.DDL__DVCuts("DiLepBaseCuts",
                              rDVMax                = 300.,
                              zDVMax                = 300.,
                              chisqPerDofMax        = 5.,
                              distMin               = 4.,
                              LowMass               = 6000.,
                              DVMassMin             = 10000.,
                              MaterialMapName       = "map",
                              MaterialMapFile       = "materialMap3D_Run2_v2.1.1.root",
                              DisabledModuleMapName = "PIXVetoMap",
                              DisabledModuleMapFile = "DisabledModuleMap_Run2_v2.root")

ToolSvc += CfgMgr.JetCalibrationTool("myJESTool",
                                     IsData=False,
                                     ConfigFile="JES_MC15Prerecommendation_April2015.config",
                                     CalibSequence="JetArea_Residual_Origin_EtaJES_GSC",
                                     JetCollection="AntiKt4EMTopo")

# Add algorithm
algSeq = CfgMgr.AthSequencer("AthAlgSeq")
#algSeq += CfgMgr.DDL__LepReco(HistSvcName = svc_name)

# Include your algorithms 
#from DDL.DDLStudiesConf import DDL__Kshort_DDL
#algSeq += CfgMgr.DDL__KshortDDL("KshortDDL")
algSeq += CfgMgr.DDL__HNL("HNL")