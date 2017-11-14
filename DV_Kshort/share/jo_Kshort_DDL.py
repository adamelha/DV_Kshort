##############################################################################################
#                                                                                            #
#     Job option file for a full athena environment for Kshort reconstruction in DDL package #
#                                                                                            #
##############################################################################################

# Place this file in DDLStudies/share/

# Enable file reading
import AthenaPoolCnvSvc.ReadAthenaPool
from AthenaCommon.AthenaCommonFlags import athenaCommonFlags as acf


# Get input file(s) for running on the local environment. If want to run on GRID -- 
from glob import glob
acf.FilesInput =["/afs/cern.ch/work/j/jbiswal/public/DAOD_RPVLL.11691312._000011.pool.root.1"]
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
#theApp.EvtMax = vars().get("maxEvents", -1)
theApp.EvtMax = 50
svcMgr.EventSelector.SkipEvents = vars().get("skipEvents", 0)

# Setup trigger tool
#ToolSvc += CfgMgr.Trig__TrigDecisionTool("TrigDecisionTool",
#                                         OutputLevel = WARNING,
#                                         TrigDecisionKey = "xTrigDecision")

# Setup THistSvc
svc_name = "Kshort_DDL"
svcMgr += CfgMgr.THistSvc()
svcMgr.THistSvc.Output += ["%s DATAFILE='hist_Kshort_DDL_2.root' OPT='RECREATE'" % (svc_name)]

# Required for Kshort reconstruction?
#trackPtCut = 10000.

# Setup PionCuts


# Setup KshortCuts


# Setup DVCuts - This was deleted for Kshort
#                              DisabledModuleMapFile = "DisabledModuleMap_Run2_v2.root")
#ToolSvc += CfgMgr.JetCalibrationTool("myJESTool",
#                                     IsData=False,
#                                     ConfigFile="JES_MC15Prerecommendation_April2015.config",
#                                     CalibSequence="JetArea_Residual_Origin_EtaJES_GSC",
#                                     JetCollection="AntiKt4EMTopo")

# Add algorithm
algSeq = CfgMgr.AthSequencer("AthAlgSeq")

# Include your algorithms 
algSeq += CfgMgr.DDL__Kshort_DDL("Kshort_DDL")
