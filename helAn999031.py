from ROOT import TCanvas, TH1D, TH2F
import ROOT

ROOT.xAOD.Init().isSuccess()


c1 = TCanvas('c1', '', 1000, 1000)

h1 = TH1D('h1', 'Ratio Cos(Helicity) 999031', 30, -1.5, 1.5)
h2 = TH1D('h2', 'Higgs pT (999031)', 100, 0, 100)
h3 = TH1D('h3', ' ', 30, -1.5, 1.5)
sc = TH2F('sc', 'Cos(Helicity) v. Mu+ - Higgs pT Average (999031)', 20, -1.5, 1.5, 100, 0, 100)
scU = TH2F('sc', 'Cos(Helicity) v. (Mu+ - Upsilon pT)/2 (999031)', 20, -1, 1, 40, -10, 4)
scE = TH2F('sc', 'Cos(Helicity) v. Pseudorapidity (999031)', 20, -1.5, 1.5, 30, -4, 6)

path = 'root://eosatlas.cern.ch//eos/atlas/user/d/daits/mc16_13TeV/recon/mc16.999031.P8BEG_23lo_ggX18p4_Upsilon1Smumu_4mu_3pt2.merge.AOD.e8304_a875_r10724_r10726_pUM999999/'

filelist = [path+'merge.AOD.999031._00001.pool.root.1', path+'merge.AOD.999031._00002.pool.root.1', path+'merge.AOD.999031._00003.pool.root.1', path+'merge.AOD.999031._00004.pool.root.1', path+'merge.AOD.999031._00005.pool.root.1', path+'merge.AOD.999031._00006.pool.root.1', path+'merge.AOD.999031._00007.pool.root.1', path+'merge.AOD.999031._00008.pool.root.1', path+'merge.AOD.999031._00009.pool.root.1', path+'merge.AOD.999031._00010.pool.root.1']

fer = open('runevnt.txt', 'r')
run_evnt_set = set()

for line in fer.readlines():
    spl = line.split()
    run = int(spl[0])
    event = int(spl[1])
    run_evnt_set.add((run,event))

def coshel(particle, parent, gparent):
    particle = ROOT.TLorentzVector(particle)
    gparent = ROOT.TLorentzVector(gparent)

    boosttoparent = -(parent.BoostVector())

    particle.Boost(boosttoparent)
    gparent.Boost(boosttoparent)
    #parent.Boost(boosttoparent)

    particle3 = particle.Vect()
    gparent3 = gparent.Vect()
    numerator = particle3.Dot(gparent3)
    denominator = (particle3.Mag())*(gparent3.Mag())
    temp = numerator/denominator

    #print parent.Px()

    return temp


for f in filelist:
    fi = ROOT.TFile.Open(f)
    t = ROOT.xAOD.MakeTransientTree(fi, 'CollectionTree')

    nEntries = t.GetEntries()

    for evnt in xrange(nEntries):
        t.GetEntry(evnt)

        runNum = t.EventInfo.auxdataConst('uint')('runNumber')
        evntNum = t.EventInfo.auxdataConst('ulonglong')('eventNumber')

        for j in t.TruthParticles:

            if j.pdgId()==35:
                ptH1 = j.pt()/1000
                etaH1 = j.eta()
                phiH1 = j.phi()
                mH1 = j.m()/1000

            if j.nParents()==0:
                continue

            if j.pdgId()==553:
                ptU1 = j.pt()/1000
                etaU1 = j.eta()
                phiU1 = j.phi()
                mU1 = j.m()/1000


            if j.absPdgId()==13 and j.parent(0).pdgId()==553:

                if j.charge()==+1:
                    ptM1 = j.pt()/1000
                    etaM1 = j.eta()
                    phiM1 = j.phi()
                    mM1 = j.m()/1000

        Higgs1 = ROOT.TLorentzVector()
        Upsi1 = ROOT.TLorentzVector()
        muon1 = ROOT.TLorentzVector()
        Higgs1.SetPtEtaPhiM(ptH1, etaH1, phiH1, mH1)
        Upsi1.SetPtEtaPhiM(ptU1, etaU1, phiU1, mU1)
        muon1.SetPtEtaPhiM(ptM1, etaM1, phiM1, mM1)
        
        """
        h3.Fill(coshel(muon1, Upsi1, Higgs1))
        h2.Fill(ptH1)
        sc.Fill(coshel(muon1, Upsi1, Higgs1), (ptM1-ptH1)/2)
        scU.Fill(coshel(muon1, Upsi1, Higgs1), (ptU1)/2)
        scE.Fill(coshel(muon1, Upsi1, Higgs1), etaM1)
        """

        if (runNum, evntNum) not in run_evnt_set: #continue if truth and recon run number, event number set matches
            continue

        for j in t.TruthParticles:
                        
            if j.pdgId()==35:
                ptH = j.pt()/1000
                etaH = j.eta()
                phiH = j.phi()
                mH = j.m()/1000
        
            if j.nParents()==0:
                continue

            if j.pdgId()==553:
                ptU = j.pt()/1000
                etaU = j.eta()
                phiU = j.phi()
                mU = j.m()/1000
            

            if j.absPdgId()==13 and j.parent(0).pdgId()==553:
            
                if j.charge()==+1:
                    ptM = j.pt()/1000
                    etaM = j.eta()
                    phiM = j.phi()
                    mM = j.m()/1000

        Higgs = ROOT.TLorentzVector()
        Upsi = ROOT.TLorentzVector()
        muon = ROOT.TLorentzVector()
        Higgs.SetPtEtaPhiM(ptH, etaH, phiH, mH)
        Upsi.SetPtEtaPhiM(ptU, etaU, phiU, mU)
        muon.SetPtEtaPhiM(ptM, etaM, phiM, mM)
        
        h1.Fill(coshel(muon, Upsi, Higgs))
        scU.Fill(coshel(muon, Upsi, Higgs), (ptM - ptU)/2)
        scE.Fill(coshel(muon, Upsi, Higgs), etaM)

"""
h1.GetXaxis().SetTitle('Cos(Helicity Angle)')
h1.Draw()
c1.SaveAs('HelAn999031.pdf')

h2.GetXaxis().SetTitle('Higgs pT [GeV]')
h2.SetLineColor(ROOT.kRed)
h2.Draw()
c1.SaveAs('999ParentpT_truth.png')
"""
hratio = h1.Clone('ratio')
hratio.Divide(h3)
hratio.GetXaxis().SetTitle('Cos(Helicity Angle)')
hratio.Draw()
c1.SaveAs('999031ratioHelAn.pdf')
"""
sc.GetXaxis().SetTitle('Cos(Helicity Angle)')
sc.GetYaxis().SetTitle('Mu+ - Higgs pT Average [GeV]')
sc.Draw()
c1.SaveAs('999031_HelAnMuPt2.png')

scU.GetXaxis().SetTitle('Cos(Helicity Angle)')
scU.GetYaxis().SetTitle('(Mu+ - Upsilon pT)/2 [GeV]')
scU.Draw()
c1.SaveAs('999031_HelAnUpPt.png')

scE.GetXaxis().SetTitle('Cos(Helicity Angle)')
scE.GetYaxis().SetTitle('Mu+ Eta')
scE.Draw('colz')
c1.SaveAs('999031_HelAnMuEta.png')
"""

ROOT.xAOD.ClearTransientTrees()
