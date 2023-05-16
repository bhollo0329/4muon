from ROOT import TCanvas, TH1D, TH2D, TLegend, TPad
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1)
ROOT.gROOT.LoadMacro("AtlasStyle.C")
ROOT.gROOT.LoadMacro("AtlasLabels.C")
ROOT.SetAtlasStyle()

ROOT.xAOD.Init().isSuccess()

c1 = TCanvas('c1',' ' , 1000, 1000)
c1.Draw()

h1 = TH2D('h1', 'Combined vs Segement Tagged', 5, -0.5, 4.5, 5, -0.5, 4.5)
h0 = TH1D('h0', 'Number of 0 Charge Quads', 60, 0, 30)
h2 = TH1D('h2', 'Number of +-2 Charge Quads', 60, 0, 30)
h4 = TH1D('h4', 'Number of +-4 Charge Quads', 60, 0, 30)

f = ROOT.TFile.Open('root://eosatlas.cern.ch//eos/atlas/user/d/daits/mc16_13TeV/DAOD_BPHY4/mc16.999031.P8BEG_23lo_ggX18p4_Upsilon1Smumu_4mu_3pt2.deriv.DAOD_BPHY4.e8304_a875_r10724_r10726_p3712_pUM999999/DAOD_BPHY4.999031._000001.pool.root.1')
t = ROOT.xAOD.MakeTransientTree(f, 'CollectionTree')

nEntries = t.GetEntries()

n0Quad = 0 #set values to 0 BEFORE event so it can add to the value after each event instead of resetting to 0 every event
n2Quad = 0
n4Quad = 0

for event in xrange(nEntries):
    t.GetEntry(event)

    bphy4_muon_lookup = { muon.auxdataConst('int')('BPHY4MuonIndex'): muon
                          for muon in t.Muons if muon.auxdataConst('int')('BPHY4MuonIndex') != -1 } #Python dicitonary that gives the muon object given the BPHY4MuonIndex value


    for quad in t.BPHY4Quads:
        tCharge = 0
        value = quad.auxdataConst('string')('CombinationCode') #accesses the muon decorator CombinationCode
        value2 = quad.auxdataConst('string')('ChargeCode') #accesses the muon decorator ChargeCode
        value3 = quad.auxdata('float')('QUAD_mass')

        for c in value2:
            if c == '+': #adds 1 to tCharge for every + character
                tCharge += 1
            elif c == '-': #subtracts 1 from tCharge for every - character
                tCharge -= 1

    if tCharge == 0: #must be outside of 'for c' loop so it can define tCharge once for each quad
        n0Quad += 1
        h0.Fill(value3/1000)
    elif tCharge == -2 or tCharge == +2:
        n2Quad += 1
        h2.Fill(value3/1000)
    elif tCharge == -4 or tCharge == +4:
        n4Quad += 1
        h4.Fill(value3/1000)
                
        #nComb = 0
        #nSegTag = 0
        if len(value) == 4:
            for c in value: #seperates string into each character
                muon = bphy4_muon_lookup[int(c)] #converts characters to strings and indexes object Muon
                print(muon.muonType()) #call the type of muon: combined or segment
                
                #if muon.muonType() == 0:
                    #nComb += 1
                #elif muon.muonType() == 2:
                    #nSegTag += 1
               # else:
                   # print(muon.muonType())
            #h1.Fill(nComb, nSegTag)
                
    
            #print(nComb)

t = (n0Quad, n2Quad, n4Quad)
tQuad = float(sum(t))
frac0 = n0Quad/tQuad
frac2 = n2Quad/tQuad
frac4 = n4Quad/tQuad

#print frac0, frac2, frac4
#print n4Quad

h0.SetLineColor(ROOT.kBlue)
h2.SetLineColor(ROOT.kRed)
h4.SetLineColor(ROOT.kGreen)
h0.GetXaxis().SetTitle('Charge Quad Mass [GeV]')
leg = TLegend(0.63, 0.60, 0.78, 0.75)
leg.AddEntry(h0, 'Charge 0', 'l')
leg.AddEntry(h2, 'Charge 2', 'l')
leg.AddEntry(h4, 'Charge 4', 'l')
h0.Draw()
h2.Draw('same')
h4.Draw('same')
leg.Draw()
c1.SaveAs('chargeQuadmass.png')

h1.GetXaxis().SetTitle('# of Combined Muons')
h1.GetYaxis().SetTitle('# of Segment Tagged Muons')
#h1.Draw('COLZ')
#h1.Print('all')


#c1.SaveAs('Com_SegTag.png')

ROOT.xAOD.ClearTransientTrees()
f.Close()


