##############################################################################################
#                                                                                            #
#     Job option file for a full athena environment for Ks reconstruction in DDL package #
#                                                                                            #
##############################################################################################

#GRL
from GoodRunsLists.GoodRunsListsConf import *
ToolSvc += GoodRunsListSelectorTool()
GoodRunsListSelectorTool.GoodRunsListVec = [ '/afs/cern.ch/user/n/nkonstan/Our_DDL/WorkArea/run/data16_13TeV_GRL.xml' ]

from GoodRunsListsUser.GoodRunsListsUserConf import *



# Enable file reading
import AthenaPoolCnvSvc.ReadAthenaPool
from AthenaCommon.AthenaCommonFlags import athenaCommonFlags as acf


# Get input file(s) for running on the local environment. If want to run on GRID -- 
from glob import glob
acf.FilesInput = []
acf.FilesInput += ["/afs/cern.ch/work/j/jbiswal/public/DAOD_RPVLL.11691312._000011.pool.root.1"]
#acf.FilesInput += ["/afs/cern.ch/work/j/jbiswal/public/data16_13TeV/DAOD_RPVLL.10626848._038273.pool.root.1"]
#acf.FilesInput += ["/afs/cern.ch/work/j/jbiswal/public/data16_13TeV/DAOD_RPVLL.10626848._038287.pool.root.1"]


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

# Setup THistSvc
svc_name = "Ks"
svcMgr += CfgMgr.THistSvc()
svcMgr.THistSvc.Output += ["%s DATAFILE='hist_Ks.root' OPT='RECREATE'" % (svc_name)]


# Add algorithm
#algSeq = CfgMgr.AthSequencer("AthAlgSeq")

# Include your algorithms 
#algSeq += CfgMgr.DDL__Ks("Ks")

job = AlgSequence()
job += CfgMgr.AthSequencer("ModSequence1")
job.ModSequence1 += GRLTriggerSelectorAlg('GRLTriggerAlg1')
job.ModSequence1.GRLTriggerAlg1.GoodRunsListArray = ['data16_13TeV_GRL']

job.ModSequence1+=CfgMgr.DDL__Ks("Ks")

