from ROOT import TCanvas, TH1D, TLegend, TPad
import ROOT
import math
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1)
ROOT.gROOT.LoadMacro("AtlasStyle.C")
ROOT.gROOT.LoadMacro("AtlasLabels.C")
ROOT.SetAtlasStyle()

ROOT.xAOD.Init().isSuccess()

path = 'root://eosatlas.cern.ch//eos/atlas/user/d/daits/mc16_13TeV/DAOD_BPHY4/mc16.999031.P8BEG_23lo_ggX18p4_Upsilon1Smumu_4mu_3pt2.deriv.DAOD_BPHY4.e8304_a875_r10724_r10726_p3712_pUM999999/'
#path = 'root://eosatlas.cern.ch//eos/atlas/user/d/daits/mc16_13TeV/DAOD_BPHY4/mc16.988015.MGPy8EG_23lo_S4b18p4NP01_Upsi1S2mu_4mu.deriv.DAOD_BPHY4.e8304_a875_r10724_r10726_p3712_pUM999999/'
#path = 'root://eosatlas.cern.ch//eos/atlas/user/d/daits/mc16_13TeV/DAOD_BPHY4/mc16.988016.MGPy8EG_23lo_P4b18p4NP01_Upsi1S2mu_4mu.deriv.DAOD_BPHY4.e8304_a875_r10724_r10726_p3712_pUM999999/'
#path = 'root://eosatlas.cern.ch//eos/atlas/user/d/daits/mc16_13TeV/DAOD_BPHY4/mc16.988017.MGPy8EG_23lo_S4b18p4NP01ptj09p2_Upsi1S2mu_4mu.deriv.DAOD_BPHY4.e8304_a875_r10724_r10726_p3712_pUM999999/'
#path = 'root://eosatlas.cern.ch//eos/atlas/user/d/daits/mc16_13TeV/DAOD_BPHY4/mc16.988018.MGPy8EG_23lo_P4b18p4NP01ptj09p2_Upsi1S2mu_4mu.deriv.DAOD_BPHY4.e8304_a875_r10724_r10726_p3712_pUM999999/'

filelist = [path+'DAOD_BPHY4.999031._000001.pool.root.1', path+'DAOD_BPHY4.999031._000002.pool.root.1']
#filelist = [path+'DAOD_BPHY4.988015._000001.pool.root.1', path+'DAOD_BPHY4.988015._000002.pool.root.1']
#filelist = [path+'DAOD_BPHY4.988016._000001.pool.root.1', path+'DAOD_BPHY4.988016._000002.pool.root.1']
#filelist = [path+'DAOD_BPHY4.988017._000001.pool.root.1', path+'DAOD_BPHY4.988017._000002.pool.root.1']
#filelist = [path+'DAOD_BPHY4.988018._000001.pool.root.1', path+'DAOD_BPHY4.988018._000002.pool.root.1']

fer = open('runevnt.txt', 'w+') #open run-event number text file to read and write

for f in filelist:
        fi = ROOT.TFile.Open(f)
        t = ROOT.xAOD.MakeTransientTree(fi, 'CollectionTree')

        nEntries = t.GetEntries()

        runlist = []
        evntlist = []

        for event in xrange(nEntries):
                t.GetEntry(event)

                runNum = t.EventInfo.auxdataConst('uint')('runNumber')
                runlist.append(runNum)

                evntNum = t.EventInfo.auxdataConst('ulonglong')('eventNumber')
                evntlist.append(evntNum)
                
                fer.write('%d %d\n' % (runNum, evntNum))  # makes a placeholder for two numbers and writes the value in each event

fer.close()

#print runlist, evntlist

ROOT.xAOD.ClearTransientTrees()

