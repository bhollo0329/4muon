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
h3 = TH1D('h3', '4 Mu Mass', 120, 0, 30)

f = ROOT.TFile.Open('root://eosatlas.cern.ch//eos/atlas/user/d/daits/mc16_13TeV/DAOD_BPHY4/mc16.999031.P8BEG_23lo_ggX18p4_Upsilon1Smumu_4mu_3pt2.deriv.DAOD_BPHY4.e8304_a875_r10724_r10726_p3712_pUM999999/DAOD_BPHY4.999031._000001.pool.root.1')
t = ROOT.xAOD.MakeTransientTree(f, 'CollectionTree')

nEntries = t.GetEntries()

for event in xrange(nEntries):
    t.GetEntry(event)
    
    for quad in t.BPHY4Quads:
        value2 = quad.auxdata('float')('QUAD_mass')
        value = quad.auxdataConst('string')('CombinationCode')
        h2.Fill(value2/1000)
        if len(value) > 4:
            print(value)

    BPHY4Quads = sorted(t.BPHY4Quads, key = lambda x: x.chiSquared()) #return a list of quads sorted by the chiSquared value

    
    if len(BPHY4Quads) > 0: # as long there's at least one
        value1 = BPHY4Quads[0].auxdata('float')('QUAD_mass') #use the first element of the list, which will have minimal value of chiSquard
        h1.Fill(value1/1000)

   # for quad in BPHY4Quads[1:]: #use non-0th-elements in the list, which will be the non best chiSquared
    #    value2 = quad.auxdata('float')('QUAD_mass')
     #   h2.Fill(value2/1000)

         
    h1.SetLineColor(ROOT.kBlue)
    h2.SetLineColor(ROOT.kRed)

h2.GetXaxis().SetTitle('4 Muon Mass [GeV]')
leg = TLegend(0.63, 0.60, 0.78, 0.75)
leg.AddEntry(h1, 'Best chi^2', 'l')
leg.AddEntry(h2, 'No chi^2 selection', 'l')
h2.Draw()
h1.Draw('same')
leg.Draw()

ROOT.xAOD.ClearTransientTrees()
f.Close()

#c1.SaveAs('quadmass_Compare.png')
c1.SaveAs('quadmass_chiCom2.png')
