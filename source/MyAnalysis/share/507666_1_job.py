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

# Set up the sample handler object.
import os
sh = ROOT.SH.SampleHandler()
sh.setMetaString( 'nc_tree', 'CollectionTree' )

inputFilePath = '/data/bhollo/SampleT2022/'+'/user.abarton.valid1.507666.MGPy8EG_23lo_P4b18p4NJ0_Upsi1S2mu_4mu_3pt2.recon.AOD.e8304_s3775_r13670.BPHY4.1_EXT0/'

ROOT.SH.ScanDir().scan( sh, inputFilePath )
"""
# .filePattern( 'file.root' ) is to access a specfic sample; not neccesary if want all samples
ROOT.SH.ScanDir().filePattern( 'AOD.27774087._000001.pool.root.1' ).scan( sh, inputFilePath )
"""
sh.printContent()

dataType = "mc"

# Create an EventLoop job.
job = ROOT.EL.Job()
job.sampleHandler( sh )
job.options().setDouble( ROOT.EL.Job.optMaxEvents, 70000 )
job.options().setString( ROOT.EL.Job.optSubmitDirMode, 'unique-link')

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

