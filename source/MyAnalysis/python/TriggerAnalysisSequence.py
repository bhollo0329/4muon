from AnaAlgorithm.AnaAlgSequence import AnaAlgSequence
from AnaAlgorithm.DualUseConfig import createAlgorithm, createPublicTool, addPrivateTool, createService

"""Create a basic trigger analysis algorithm sequence w/ 
      Keyword arguments:
      datatype - the data type to run on ("data" or "mc")
"""

def makeTriggerAnalysisSequence( dataType ):

    algSeq = AnaAlgSequence( "TriggerAnalysisSequence" )

    # Set up the systematics loader/handler service
    sysService = createService( 'CP::SystematicsSvc', 'SystematicsSvc', sequence = algSeq )
    sysService.sigmaRecommended = 1
    
    #create tdt tool and add to sequence 
    tdt = createPublicTool( "Trig::TrigDecisionTool", "TrigDecisionTool" )
    algSeq.addPublicTool( tdt ) 
    addPrivateTool(tdt, "ConfigTool", "TrigConf::xAODConfigTool")

    tdt.NavigationFormat = "TrigComposite"
    tdt.HLTSummary = "HLTNav_Summary_AODSlimmed"

    alg = createAlgorithm( 'MyxAODAnalysis', 'AnalysisAlg')

    #assign tdt to algorithm as string
    alg.TrigDecisionTool = "%s/%s" % ( tdt.getType(), tdt.getName() )

    #add algorithm to sequence
    algSeq += alg

    return algSeq
