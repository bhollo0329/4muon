#!/usr/bin/env python

# Read the submission directory as a command line argument.
import optparse
parser = optparse.OptionParser()
parser.add_option( '-s', '--submission-dir', dest = 'submission_dir',
                   action = 'store', type = 'string', default = 'submitDir',
                   help = 'Submission directory for EventLoop' )
( options, args ) = parser.parse_args()

# Set up (Py)ROOT.
import ROOT
ROOT.xAOD.Init().ignore()

inputFilePath = '/data/bhollo/SampleT2022/'+'/user.abarton.valid1.800852.P8BEG_23lo_ggX18p4_Upsilon1Smumu_4mu_3pt2.recon.AOD.e8304_s3775_r13670.BPHY4.1_EXT0/'
#inputFilePath = '/data_ceph/onyisi/4mu_22/'+'/valid1.800852.P8BEG_23lo_ggX18p4_Upsilon1Smumu_4mu_3pt2.recon.AOD.e8304_s3775_r13670/'
#inputFilePath = '/data/bhollo/SampleT2022/'+'mc20_13TeV.800852.P8BEG_23lo_ggX18p4_Upsilon1Smumu_4mu_3pt2.merge.AOD.e8304_s3681_d1735_r13308_r13307'

# Set up the sample handler object.
import os
sh = ROOT.SH.SampleHandler()
sh.setMetaString( 'nc_tree', 'CollectionTree' )

#for one specific sample use ROOT.SH.ScanDir().filePattern( 'file.root').scan ( sh, inputFilepath)
ROOT.SH.ScanDir().scan( sh, inputFilePath )
sh.printContent()

dataType = "mc"

# Create an EventLoop job.
job = ROOT.EL.Job()
job.sampleHandler( sh )
job.options().setDouble( ROOT.EL.Job.optMaxEvents, 100000 )
job.options().setString( ROOT.EL.Job.optSubmitDirMode, 'unique-link')

#Add trigger algorithm to jobs
from MyAnalysis.TriggerAnalysisSequence import makeTriggerAnalysisSequence
algSeq = makeTriggerAnalysisSequence( dataType )

for alg in algSeq:
    job.algsAdd( alg )
"""
# Add our algorithm to the job
job.algsAdd( alg)
"""     
# Run the job using the direct driver.
driver = ROOT.EL.DirectDriver()
driver.submit( job, options.submission_dir )

