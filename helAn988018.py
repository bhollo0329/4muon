from ROOT import TCanvas, TH1D, TH2F, TProfile
import ROOT

ROOT.xAOD.Init().isSuccess()


c1 = TCanvas('c1', '', 1000, 1000)

h1 = TH1D('h1', 'Ratio Cos(Helicity) 988018', 30, -1.5, 1.5)
h2 = TH1D('h2', 'Tetraquark pT (988018)', 100, 0, 100)
h3 = TH1D('h3', ' ', 30, -1.5, 1.5)
h4 = TH1D('h4', 'Cos(Helicity) 988018 Parent pT Cut Ratio', 30, -1.5, 1.5)
sc = TH2F('sc', 'Cos(Helicity) v. Mu+ - Tetraquark pT Average (988018)', 20, -1.5, 1.5, 30, -10, 4)
scU = TH2F('scU', 'Cos(Helicity) v. (Mu+ - Upsilon pT)/2 (988018)', 20, -1, 1, 40, -8, 3)
scE = TH2F('sc', 'Cos(Helicity) v. Mu+ Pseudorapidity (988018)', 20, -1, 1, 30, -4, 4)
pro = TProfile('hpro', 'Cos(Helicity) v. Muon+-Tetraquark pT Average Profile (988018)', 20, -1, 1)
pro1 = TProfile('hpro1', 'Cos(Helicity) v. Muon+-(Upsilon pT/2) Relativistic Profile (988018)', 20, -1, 1)
pro2 = TProfile('hpro2', 'Cos(Helicity) v. Muon+ Pseudorapidity Profile (988018)', 20, -1, 1)

path = 'root://eosatlas.cern.ch//eos/atlas/user/d/daits/mc16_13TeV/recon/mc16.988018.MGPy8EG_23lo_P4b18p4NP01ptj09p2_Upsi1S2mu_4mu.merge.AOD.e8304_a875_r10724_r10726_pUM999999/'

filelist = [path+'merge.AOD.988018._00001.pool.root.1', path+'merge.AOD.988018._00002.pool.root.1', path+'merge.AOD.988018._00003.pool.root.1', path+'merge.AOD.988018._00004.pool.root.1', path+'merge.AOD.988018._00005.pool.root.1', path+'merge.AOD.988018._00006.pool.root.1', path+'merge.AOD.988018._00007.pool.root.1', path+'merge.AOD.988018._00008.pool.root.1', path+'merge.AOD.988018._00009.pool.root.1']

fer = open('runevnt018.txt', 'r')
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

    particle3 = particle.Vect()
    gparent3 = gparent.Vect()
    numerator = particle3.Dot(gparent3)
    denominator = (particle3.Mag())*(gparent3.Mag())
    temp = numerator/denominator

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

            if j.pdgId()==90000025:
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
        
        Tetra1 = ROOT.TLorentzVector()
        Upsi1 = ROOT.TLorentzVector()
        muon1 = ROOT.TLorentzVector()
        Tetra1.SetPtEtaPhiM(ptH1, etaH1, phiH1, mH1)
        Upsi1.SetPtEtaPhiM(ptU1, etaU1, phiU1, mU1)
        muon1.SetPtEtaPhiM(ptM1, etaM1, phiM1, mM1)

        """
        h3.Fill(coshel(muon1, Upsi1, Tetra1))
        h2.Fill(ptH1)
        sc.Fill(coshel(muon1, Upsi1, Tetra1), (ptM1-ptH1)/2)
        scU.Fill(coshel(muon1, Upsi1, Tetra1), ptU1)
        scE.Fill(coshel(muon1, Upsi1, Tetra1), etaM1)
        pro.Fill(coshel(muon1, Upsi1, Tetra1), (ptM1-ptH1)/2, 1)
        """

        if (runNum, evntNum) not in run_evnt_set: #continue if truth and recon run number, event number set matches
            continue

        for j in t.TruthParticles:

            if j.pdgId()==90000025:
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

        Tetra = ROOT.TLorentzVector()
        Upsi = ROOT.TLorentzVector()
        muon = ROOT.TLorentzVector()
        Tetra.SetPtEtaPhiM(ptH, etaH, phiH, mH)
        Upsi.SetPtEtaPhiM(ptU, etaU, phiU, mU)
        muon.SetPtEtaPhiM(ptM, etaM, phiM, mM)

        h1.Fill(coshel(muon, Upsi, Tetra))
        if ptU > 10:
            scU.Fill(coshel(muon, Upsi, Tetra),(ptM- ptU)/2)
            pro1.Fill(coshel(muon, Upsi, Tetra), ptM-(ptU/2), 1)
        scE.Fill(coshel(muon, Upsi, Tetra), etaM)
        #pro1.Fill(coshel(muon, Upsi, Tetra), ptM-(ptU/2), 1)
        pro2.Fill(coshel(muon, Upsi, Tetra), etaM)
        if ptH<=10 and ptH>=5:
            h4.Fill(coshel(muon, Upsi, Tetra))
        
"""
h1.GetXaxis().SetTitle('Cos(Helicity Angle)')
h1.Draw()
c1.SaveAs('HelAn988018.png')

h2.GetXaxis().SetTitle('Tetraquark pT [GeV]')
h2.SetLineColor(ROOT.kRed)
h2.Draw()
c1.SaveAs('018ParentpT_truth.png')

hratio = h1.Clone('ratio')
hratio.Divide(h3)
hratio.GetXaxis().SetTitle('Cos(Helicity Angle)')
hratio.Draw()
c1.SaveAs('988018ratioHelAn.png')

h4.GetXaxis().SetTitle('Cos(Helicity Angle)')
h4.Draw()
c1.SaveAs('HelAn988018_pTCut.png')

sc.GetXaxis().SetTitle('Cos(Helicity Angle)')
sc.GetYaxis().SetTitle('Mu+ - Tetraquark pT Average [GeV]')
sc.Draw()
c1.SaveAs('988018_HelAnMuPt2.png')

hr = h4.Clone('ratio')
hr.Divide(h1)
hr.GetXaxis().SetTitle('Cos(Helicity Angle)')
hr.Draw()
c1.SaveAs('018helAnRatio_pTcut.png')

pro.GetXaxis().SetTitle('Cos(Helicity Angle)')
pro.Draw()
c1.SaveAs('018prohelAnPt.png')
"""
scE.GetXaxis().SetTitle('Cos(Helicity Angle)')
scE.GetYaxis().SetTitle('Mu+ Eta')
scE.Draw()
c1.SaveAs('988018_HelAnMuEta.png')
"""
scU.GetXaxis().SetTitle('Cos(Helicity Angle)')
scU.GetYaxis().SetTitle('(Mu+ - Upsilon pT)/2 [GeV]')
scU.Draw()
c1.SaveAs('988018_HelAnUpT2.png')

pro1.GetXaxis().SetTitle('Cos(Helicity Angle)')
pro1.Draw()
c1.SaveAs('018ProHelAnUpPtRel.png')

pro2.GetXaxis().SetTitle('Cos(Helicity Angle)')
pro2.Draw()
c1.SaveAs('018ProHelAnMuEta.png')
"""


ROOT.xAOD.ClearTransientTrees()
