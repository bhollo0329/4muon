from ROOT import TCanvas, TH1D, TH2D, TLegend, TPad
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1)
ROOT.gROOT.LoadMacro("AtlasStyle.C")
ROOT.gROOT.LoadMacro("AtlasLabels.C")
#ROOT.SetAtlasStyle()

ROOT.xAOD.Init().isSuccess()

c1 = TCanvas('c1',' ' , 1200, 1000)

path = 'root://eosatlas.cern.ch//eos/atlas/user/d/daits/mc16_13TeV/recon/mc16.988018.MGPy8EG_23lo_P4b18p4NP01ptj09p2_Upsi1S2mu_4mu.merge.AOD.e8304_a875_r10724_r10726_pUM999999/'


f = ROOT.TFile.Open(path+'merge.AOD.988018._00001.pool.root.1')
t = ROOT.xAOD.MakeTransientTree(f, 'CollectionTree')

nEntries = t.GetEntries()

for event in xrange(nEntries):
    t.GetEntry(event)

    for j in t.TruthParticles:

        if j.nParents()==0:
            continue

        if j.absPdgId()==13:

            if j.parent(0).pdgId() == 553 or j.parent(0).absPdgId() == 13:
                continue

            print j.parent(0).pdgId()

ROOT.xAOD.ClearTransientTrees()
