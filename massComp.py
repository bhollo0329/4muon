from ROOT import TCanvas, TH1D, TLegend, TPad
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1)
ROOT.gROOT.LoadMacro("AtlasStyle.C")
ROOT.gROOT.LoadMacro("AtlasLabels.C")
ROOT.SetAtlasStyle()

ROOT.xAOD.Init().isSuccess()

c1 = TCanvas('c1',' ' , 1000, 1000)
c1.Draw()

h1 = TH1D('h1', '4 Muon Mass Best Chi^2', 120, 0, 30)
h2 = TH1D('h2', '4 Mu Mass Worse Chi^2', 120, 0, 30)
h3 = TH1D('h3', '4 Mu Mass', 120, -30, 30)

f = ROOT.TFile.Open('root://eosatlas.cern.ch//eos/atlas/user/d/daits/mc16_13TeV/DAOD_BPHY4/mc16.999031.P8BEG_23lo_ggX18p4_Upsilon1Smumu_4mu_3pt2.deriv.DAOD_BP\
HY4.e8304_a875_r10724_r10726_p3712_pUM999999/DAOD_BPHY4.999031._000001.pool.root.1')
t = ROOT.xAOD.MakeTransientTree(f, 'CollectionTree')

nEntries = t.GetEntries()

for event in xrange(nEntries):
    t.GetEntry(event)

    for quad in t.BPHY4Quads:
        value1 = quad.auxdata('float')('QUAD_mass')
        h1.Fill(value1/1000)

        value2 = quad.auxdataConst('string')('CombinationCode')

        p = ROOT.TLorentzVector()

        if len(value2) == 4:
            for c in value2: #seperates string into each character
                muon = t.Muons[int(c)] 
                p += muon.p4()

    BPHY4Quads = sorted(t.BPHY4Quads, key = lambda x: x.chiSquared()) #list that iterates over objects in BPHY4Quads branch and sorts them by their chi squared value


    if len(BPHY4Quads) > 0: # as long there's at least one
        value3 = BPHY4Quads[0].auxdata('float')('QUAD_mass') #use the first element of the list, which will have minimal value of chiSquard

    h2.Fill(p.M()/1000)
    diff = value3 - p.M()
    h3.Fill(diff/1000)

#h1.SetLineColor(ROOT.kBlue)
#h2.SetLineColor(ROOT.kRed)
#h1.GetXaxis().SetTitle('4 Muon Mass [GeV]')
#leg = TLegend(0.63, 0.60, 0.78, 0.75)
#leg.AddEntry(h1, 'Quad Mass', 'l')
#leg.AddEntry(h2, '4 Momentum Mass', 'l')
#h1.Draw()
#h2.Draw('same')
#leg.Draw()
h3.Draw()
h3.GetXaxis().SetTitle('4 Muon Mass [GeV]')

ROOT.xAOD.ClearTransientTrees()
f.Close()

#c1.SaveAs('muMass_check.png')
c1.SaveAs('muMass_diff.png')
