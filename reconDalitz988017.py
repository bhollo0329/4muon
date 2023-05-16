from ROOT import TCanvas, TH1D, TH2D, TLegend, TPad
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1)
ROOT.gROOT.LoadMacro("AtlasStyle.C")
ROOT.gROOT.LoadMacro("AtlasLabels.C")
#ROOT.SetAtlasStyle()

ROOT.xAOD.Init().isSuccess()

c1 = TCanvas('c1',' ' , 1200, 1000)

hd = TH2D('h2', '988017 Dalitz', 60, 50, 350, 18, 0, 90)
hdr = TH2D('h2', 'Dalitz 988017', 60, 50, 350, 18, 0, 90)

path = 'root://eosatlas.cern.ch//eos/atlas/user/d/daits/mc16_13TeV/recon/mc16.988017.MGPy8EG_23lo_S4b18p4NP01ptj09p2_Upsi1S2mu_4mu.merge.AOD.e8304_a875_r10724_r10726_pUM999999/'

filelist = [path+'merge.AOD.988017._00001.pool.root.1', path+'merge.AOD.988017._00002.pool.root.1', path+'merge.AOD.988017._00003.pool.root.1', path+'merge.AOD.988017._00004.pool.root.1', path+'merge.AOD.988017._00005.pool.root.1', path+'merge.AOD.988017._00006.pool.root.1', path+'merge.AOD.988017._00007.pool.root.1', path+'merge.AOD.988017._00008.pool.root.1', path+'merge.AOD.988017._00009.pool.root.1']

fer = open('runevnt017.txt', 'r')
run_evnt_set = set()

for line in fer.readlines():
    spl = line.split()
    run = int(spl[0])
    event = int(spl[1])
    run_evnt_set.add((run,event))

for f in filelist:
    fi = ROOT.TFile.Open(f)
    t = ROOT.xAOD.MakeTransientTree(fi, 'CollectionTree')

    nEntries = t.GetEntries()

    for event in xrange(nEntries):
        t.GetEntry(event)

        runNum = t.EventInfo.auxdataConst('uint')('runNumber')
        evntNum = t.EventInfo.auxdataConst('ulonglong')('eventNumber')

        for j in t.TruthParticles:

            if j.pdgId()==553: #upsilon 1S
                ptu = j.pt()/1000
                etau = j.eta()
                phiu = j.phi()
                mu = j.m()/1000

            if j.nParents()==0:
                continue

            if j.absPdgId()==13 and j.parent(0).pdgId() == 90000025:

                if j.charge()== +1: #Higgs daughter mu+
                    #print('mu+')
                    ptmp = j.pt()/1000
                    etamp = j.eta()
                    phimp = j.phi()
                    mmp = j.m()/1000

                elif j.charge()== -1: #Higgs daughter mu-
                    #print('mu-')
                    ptmn = j.pt()/1000
                    etamn = j.eta()
                    phimn = j.phi()
                    mmn = j.m()/1000

        pd1 = ROOT.TLorentzVector()
        pd1.SetPtEtaPhiM(ptu, etau, phiu, mu)

        pd2 = ROOT.TLorentzVector()
        pd2.SetPtEtaPhiM(ptmp, etamp, phimp, mmp)

        pd3 = ROOT.TLorentzVector()
        pd3.SetPtEtaPhiM(ptmn, etamn, phimn, mmn)

        sd1 = (pd1 + pd3).M2() #upsilon and mu1 invariant mass squared
        sd2 = (pd2 + pd3).M2() #mu+ and mu- invariant mass squared

        hd.Fill(sd1, sd2)

        if (runNum, evntNum) not in run_evnt_set: #continue if truth and recon run number, event number set matches
            continue

        for j in t.TruthParticles:
            
            if j.pdgId()==553: #upsilon 1S
                ptu2 = j.pt()/1000
                etau2 = j.eta()
                phiu2 = j.phi()
                mu2 = j.m()/1000

            if j.nParents() == 0: #if particle has 0 parents then exclude from loop
                    #print('no muon parents')
                    continue

            if j.absPdgId()==13 and j.parent(0).pdgId()== 90000025: #identifying muons origination from Higgs 0
                        #print('good muon')

                if j.charge()== +1: #Higgs daughter mu+
                    #print('mu+')
                    ptmp2 = j.pt()/1000
                    etamp2 = j.eta()
                    phimp2 = j.phi()
                    mmp2 = j.m()/1000

                elif j.charge()== -1: #Higgs daughter mu-
                    #print('mu-')
                    ptmn2 = j.pt()/1000
                    etamn2 = j.eta()
                    phimn2 = j.phi()
                    mmn2 = j.m()/1000

        pdr1 = ROOT.TLorentzVector()
        pdr1.SetPtEtaPhiM(ptu2, etau2, phiu2, mu2)

        pdr2 = ROOT.TLorentzVector()
        pdr2.SetPtEtaPhiM(ptmp2, etamp2, phimp2, mmp2)

        pdr3 = ROOT.TLorentzVector()
        pdr3.SetPtEtaPhiM(ptmn2, etamn2, phimn2, mmn2)

        sdr1 = (pdr1 + pdr3).M2() #upsilon and mu1 invariant mass squared
        sdr2 = (pdr2 + pdr3).M2() #mu+ and mu- invariant mass squared

        hdr.Fill(sdr1, sdr2)

"""
hratio = hdr.Clone('ratio')
hratio.Divide(hd)

hratio.GetXaxis().SetTitle('s(upsilon, mu-) m^2 [GeV]')
hratio.GetYaxis().SetTitle('s(mu+, mu-) m^2 [GeV]')
hratio.SetAxisRange(0, 0.8, "Z")
hratio.Draw('colz')

c1.SaveAs('ratio_Dalitz988017.png')
"""

hd.GetXaxis().SetTitle('s(upsilon, mu-) m^2 [GeV]')
hd.GetYaxis().SetTitle('s(mu+, mu-) m^2 [GeV]')
hd.Draw('colz')

c1.SaveAs('Dalitz988017.png')

ROOT.xAOD.ClearTransientTrees()

