from Gaudi.Configuration import *

from Configurables import LcioEvent, EventDataSvc, MarlinProcessorWrapper
from k4MarlinWrapper.parseConstants import *
algList = []
evtsvc = EventDataSvc()


CONSTANTS = {
}

parseConstants(CONSTANTS)

read = LcioEvent()
read.OutputLevel = INFO
read.Files = ["input_file.slcio"]
algList.append(read)

MyAIDAProcessor = MarlinProcessorWrapper("MyAIDAProcessor")
MyAIDAProcessor.OutputLevel = INFO 
MyAIDAProcessor.ProcessorType = "AIDAProcessor" 
MyAIDAProcessor.Parameters = {
                              "Compress": ["1"],
                              "FileName": ["output_def_tau_finder"],
                              "FileType": ["root"]
                              }

TauFinder = MarlinProcessorWrapper("TauFinder")
TauFinder.OutputLevel = INFO 
TauFinder.ProcessorType = "TauFinder" 
TauFinder.Parameters = {
                          "PFOCollection": ["SelectedPandoraPFOs"],
                          "TauRecCollection": ["TauPFOs"],
                          "TauRecRestCollection": ["RestPFOs"],
                          "TauRecLinkCollectionName": ["TauRecLink"],
                          "FileName_Signal": ["tau_signal.root"]
                          }

algList.append(MyAIDAProcessor)
algList.append(TauFinder)

from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = algList,
                EvtSel = 'NONE',
                EvtMax   = 10,
                ExtSvc = [evtsvc],
                OutputLevel=INFO
              )
