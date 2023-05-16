#include <AsgMessaging/MessageCheck.h>
#include <MyAnalysis/MyxAODAnalysis.h>
#include <xAODEventInfo/EventInfo.h>
#include <xAODTracking/VertexContainer.h>
#include <xAODMuon/MuonContainer.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <map>
#include <sstream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <TH1.h>
#include <TH2.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TEfficiency.h>
#include <TPaveStats.h>
#include <TStyle.h>

//Trigger include
#include <TrigDecisionTool/TrigDecisionTool.h>

//TDT includes
#include "TrigDecisionTool/ChainGroup.h"
#include "TrigDecisionTool/FeatureContainer.h"
#include "TrigDecisionTool/Feature.h"
#include "TrigDecisionTool/TDTUtilities.h"
#include "TrigDecisionTool/ExpertMethods.h"

//Trigger EDM (Event Data Model) includes
#include "TrigSteeringEvent/TrigRoiDescriptor.h"
#include "TrigSteeringEvent/TrigRoiDescriptorCollection.h"
#include "xAODTrigger/TrigPassBits.h"
#include "xAODTrigger/EmTauRoIContainer.h"
#include <TrigCompositeUtils/TrigCompositeUtils.h>
#include "TrigCompositeUtils/Combinations.h"
#include "TrigCompositeUtils/ChainNameParser.h"
#include "TriggerMatchingTool/R3MatchingTool.h"

using namespace std;

MyxAODAnalysis :: MyxAODAnalysis (const std::string& name,
                                  ISvcLocator *pSvcLocator)
  : EL::AnaAlgorithm (name, pSvcLocator),
    m_trigDec( "Trig::TrigDecisionTool/TrigDecisionTool" ),
    m_evtNr(0), M11_evntNum(0), M14_evntNum(0),
    m_all(0), m_allL1(0), m_allHLT(0), m_allHLT_mu(0), m_2mu4_7M14(0), m_2mu4_7M11(0), m_mu22_mu8noL1(0), m_mu20_mu8noL1(0), m_7M14_L1(0), 
    pEff_M11(0), pEff_M14(0)

{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  This is also where you
  // declare all properties for your algorithm.  Note that things like
  // resetting statistics variables or booking histograms should
  // rather go into the initialize() function.

  declareProperty( "TrigDecisionTool", m_trigDec, "The tool to access TrigDecision");
  declareProperty( "TriggerList", m_chain_names, "List of triggers to analyze");

}



StatusCode MyxAODAnalysis :: initialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  // Retrieve the TDT
  ANA_CHECK( m_trigDec.retrieve() );
  
  //define big chain groups relying on the name convention (ex: all L1 items start their name with L1_)
  m_all = m_trigDec->getChainGroup( ".*" );
  m_allL1 = m_trigDec->getChainGroup( "L1_.*" );
  m_allHLT = m_trigDec->getChainGroup( "HLT_.*" );
  m_allHLT_mu = m_trigDec->getChainGroup( "HLT_.*mu.*" );
  m_2mu4_7M14 = m_trigDec->getChainGroup("HLT_2mu4_L1BPH-7M14-0DR25-MU5VFMU3VF");
  m_2mu4_7M11 = m_trigDec->getChainGroup("HLT_2mu4_L1BPH-7M11-25DR99-2MU3VF");
  m_mu22_mu8noL1 = m_trigDec->getChainGroup("HLT_mu22_mu8noL1_L1MU14FCH");
  m_mu20_mu8noL1 = m_trigDec->getChainGroup("HLT_mu20_ivarmedium_mu8noL1_L1MU14FCH");
  m_7M14_L1 = m_trigDec->getChainGroup("L1_BPH-7M14-0DR25-MU5VFMU3VF");
  
  ANA_CHECK( book( TH1D( "trigFreq_800852", "trigFreq_800852", 406, 0, 1624 )));
  hist ("trigFreq_800852")->GetXaxis()->SetTitle("Fired Triggers");
  hist ("trigFreq_800852")->GetYaxis()->SetTitle("# of Events");

  ANA_CHECK( book( TH1D( "dimu_M11", "HLT_M11 Muon Oppostie-Sign Pair Invariant Mass", 16, 0, 15)));
  hist ("dimu_M11")->GetXaxis()->SetTitle("Invariant Mass [GeV]");
  //hist ("dimu_M11")->GetYaxis()->SetTitle("# of Events");

  ANA_CHECK( book( TH1D( "dimu_M14", "HLT_M14 Muon Oppostie-Sign Pair Invariant Mass", 120, 0, 15)));
  hist ("dimu_M14")->GetXaxis()->SetTitle("Invariant Mass [GeV]");
  //hist ("dimu_M14")->GetYaxis()->SetTitle("# of Events");

  ANA_CHECK( book( TH1D( "dimu_M11_OSU", "HLT_M11 Muon Oppostie-Sign Pair Upsilon Invariant Mass", 60, 0, 15)));
  hist ("dimu_M11_OSU")->GetXaxis()->SetTitle("Invariant Mass [GeV]");
  //hist ("dimu_M11_OSU")->GetYaxis()->SetTitle("# of Events");

  ANA_CHECK( book( TH1D( "dimu_M14_OSU", "HLT_M14 Muon Oppostie-Sign Pair Upsilon Invariant Mass", 150, 0, 15)));
  hist ("dimu_M14_OSU")->GetXaxis()->SetTitle("Invariant Mass [GeV]");
  //hist ("dimu_M14_OSU")->GetYaxis()->SetTitle("# of Events");


  ANA_CHECK( book( TH1D( "dimu_ComNum11", "HLT_M11 Number of Di-Muon Combinations", 11, 0, 10 )));
  hist ("dimu_ComNum11")->GetXaxis()->SetTitle("# of Combinations");
  //hist ("dimu_ComNum")->GetYaxis()->SetTitle("# of Events");

  ANA_CHECK( book( TH1D( "dimu_ComNum14", "HLT_M14 Number of Di-Muon Combinations", 11, 0, 10 )));
  hist ("dimu_ComNum14")->GetXaxis()->SetTitle("# of Combinations");

  //Trigger Efficiency histograms
  ANA_CHECK( book( TH1D( "passM11_pT", "HLT_M11 Passed pT", 30, 0, 50)));
  hist ("passM11_pT")->GetXaxis()->SetTitle("Quad pT [GeV]");

  ANA_CHECK( book( TH1D( "passM14_pT", "HLT_M14 Passed pT", 30, 0, 50)));
  hist ("passM14_pT")->GetXaxis()->SetTitle("Quad pT [GeV]");

  ANA_CHECK( book( TH1D( "passMu20_pT", "HLT_mu20_mu8noL1 Passed pT", 30, 0, 50)));
  hist ("passMu20_pT")->GetXaxis()->SetTitle("Quad pT [GeV]");

  ANA_CHECK( book( TH1D( "passMu22_pT", "HLT_mu22_mu8noL1 Passed pT", 30, 0, 50)));
  hist ("passMu22_pT")->GetXaxis()->SetTitle("Quad pT [GeV]");

  ANA_CHECK( book( TH1D( "passM14L1_pT", "L1_M14 Passed pT", 30, 0, 50)));
  hist ("passM14L1_pT")->GetXaxis()->SetTitle("Quad(X) pT [GeV]");

  ANA_CHECK( book( TH1D( "total_pT", "Total pT", 30, 0, 50)));
  hist ("total_pT")->GetXaxis()->SetTitle("Quad(X) pT [GeV]");

  //Individual quad muons kinematics
  ANA_CHECK( book( TH1D( "Xmu_eta", "Quad Muons Eta", 60, -3, 3)));
  hist ("Xmu_eta")->GetXaxis()->SetTitle("Eta");

  ANA_CHECK( book( TH1D( "quadMupT", "Quad Muons pT", 16, 2, 10)));
  hist ("quadMupT")->GetXaxis()->SetTitle("pT [GeV]");

  ANA_CHECK( book( TH1D( "lowestMupT_X", "Lowest Quad Muon pT, Low X pT ", 16, 2, 10)));
  hist ("lowestMupT_X")->GetXaxis()->SetTitle("pT [GeV]");

  ANA_CHECK( book( TH1D( "lowestMupT", "Lowest Quad Muon pT", 16, 2, 10)));
  hist ("lowestMupT")->GetXaxis()->SetTitle("pT [GeV]"); 
  
  ANA_CHECK( book( TH1D( "SecLowMu_pT", "Second Lowest Quad Muon pT", 16, 2, 10)));
  hist ("SecLowMu_pT")->GetXaxis()->SetTitle("pT [GeV]");

  ANA_CHECK( book( TH1D( "SecHighMu_pT", "Second Highest Quad Muon pT", 16, 2, 10)));
  hist ("SecHighMu_pT")->GetXaxis()->SetTitle("pT [GeV]");

  ANA_CHECK( book( TH1D( "highestMu_pT", "Highest Quad Muon pT", 16, 2, 10)));
  hist ("highestMu_pT")->GetXaxis()->SetTitle("pT [GeV]");

  ANA_CHECK( book(TH1D( "muEta_low3", "Quad Muon Eta < 3 GeV", 12, -3, 3)));
  hist ("muEta_low3")->GetXaxis()->SetTitle("eta");

  ANA_CHECK( book(TH1D( "muEta_up3", "Quad Muon Eta > 3 GeV", 12, -3, 3)));
  hist ("muEta_up3")->GetXaxis()->SetTitle("eta");

  ANA_CHECK( book(TH1D( "muTypeCom", "Combined vs Segment Muons", 5, 0, 6)));
  hist ("muTypeCom")->GetXaxis()->SetTitle("Combined Muons");

  //Muon cuts histograms
  ANA_CHECK( book( TH1D( "X_invMass", "X (Quad) Invariant Mass", 60, 15, 25)));
  hist ("X_invMass")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "X_Mdiff", "X (Quad) Decorator - p4 Invariant Mass", 25, 0, 5)));
  hist ("X_Mdiff")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH2D( "X_mResolP", "4-Momentum Invariant Mass Resolution", 60, 15, 25, 25, 0, 5)));
  hist ("X_mResolP")->GetXaxis()->SetTitle("4-momentum Invariant Mass [GeV]");
  hist ("X_mResolP")->GetYaxis()->SetTitle("Invariant Mass Difference  [GeV]");

  ANA_CHECK( book( TH2D( "X_mResolD", "Decorator Invariant Mass Resolution", 60, 15, 25, 25, 0, 5)));
  hist ("X_mResolD")->GetXaxis()->SetTitle("Decorator Invariant Mass [GeV]");
  hist ("X_mResolD")->GetYaxis()->SetTitle("Invariant Mass Difference  [GeV]");

  ANA_CHECK( book( TH1D( "M_ptCut", "X (Quad) Invariant Mass pT Selection", 60, 15, 25)));
  hist ("M_ptCut")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "M_etaCut", "X (Quad) Invariant Mass Eta Selection", 60, 15, 25)));
  hist ("M_etaCut")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "M_ptCut2", "X (Quad) Invariant Mass Second Stage pT Selection", 60, 15, 25)));
  hist ("M_ptCut2")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "M_ComCut", "X (Quad) Invariant Mass pT + Combination Selection", 60, 15, 25)));
  hist ("M_ComCut")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "M_gQuadCut", "X (Quad) Invariant Mass Good Quad Selection", 60, 15, 25)));
  hist ("M_gQuadCut")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "M_ComCutO", "X (Quad) Invariant Mass Combination Selection", 60, 15, 25)));
  hist ("M_ComCutO")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  //Muon Pair selection histograms
  ANA_CHECK( book( TH1D( "M_2muChi", "X (Quad) Invariant Mass Good Quad + Pair Mass + Chi Squared Selection", 60, 15, 25)));
  hist ("M_2muChi")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "dimuChi", "Muon Pair Chi Squared", 30, 0, 15)));
  hist ("dimuChi")->GetXaxis()->SetTitle("Chi Squared");
  
  ANA_CHECK( book( TH1D( "M_2muUpsi", "X (Quad) Invariant Mass Good Quad + Onia + Upsilon(1S) Selection", 60, 15, 25)));
  hist ("M_2muUpsi")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "M_2muIM", "X (Quad) Invariant Mass Good Quad + Pair Invariant Mass Selection", 60, 15, 25)));
  hist ("M_2muIM")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "dimuIM", "Muon Pair Invariant Mass", 50, 0, 50)));
  hist ("dimuIM")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "M_noOni", "Quad Invariant Mass No Onia Candidates", 60, 15, 25)));
  hist ("M_noOni")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "dimuIM_OS", "Opposite Sign Muon Pair Invariant Mass", 60, 0, 15)));
  hist ("dimuIM_OS")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "M_OniaCut", "X (Quad) Invariant Mass Good Quad + Onia Selection", 60, 15, 25)));
  hist ("M_OniaCut")->GetXaxis()->SetTitle("Invariant Mass [GeV]");
  
  ANA_CHECK( book( TH1D( "M_0Charge", "X (Quad) Invariant Mass Upsilon(1S) Events with Total Charge = 0", 60, 15, 25)));
  hist ("M_0Charge")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "M_2Charge", "X (Quad) Invariant Mass Upsilon(1S) Events with Total Charge = +/- 2", 60, 15, 25)));
  hist ("M_2Charge")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "M_finalIM", "X (Quad) Invariant Mass Upsilon(1S) Events Final Selection", 60, 15, 25)));
  hist ("M_finalIM")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "M_0Charge_2upsi", "X (Quad) Invariant Mass 2 Upsilon(1S) Candidates with Total Charge = 0", 60, 15, 25)));
  hist ("M_0Charge_2upsi")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "M_0Charge_1upsi", "X (Quad) Invariant Mass Single Upsilon(1S) Candidate", 60, 15, 25)));
  hist ("M_0Charge_1upsi")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "X_Onia_noUpsi", "X (Quad) Invariant Mass Onia + No Upsilon(1S) Candidate", 90, 15, 25)));
  hist ("X_Onia_noUpsi")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "X_noUpsiL", "X (Quad) Invariant Mass No Upsilon(1S) Pair < 9.2", 90, 15, 25)));
  hist ("X_noUpsiL")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "X_noUpsiH", "X (Quad) Invariant Mass No Upsilon(1S) Pair > 9.7", 90, 15, 25)));
  hist ("X_noUpsiH")->GetXaxis()->SetTitle("Invariant Mass [GeV]");


  ANA_CHECK( book( TH1D( "dimu_noUpsi", "Onia Pair Invariant Mass No Upsilon(1S)", 140, 1, 15)));
  hist ("dimu_noUpsi")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH2D( "noUpsi_1v2OS", "OS Onia 1 vs Onia 2 Invariant Mass No Upsilon(1S)", 140, 1, 15, 140, 1, 15)));
  hist ("noUpsi_1v2OS")->GetXaxis()->SetTitle("Onia 1 Invariant Mass [GeV]");
  hist ("noUpsi_1v2OS")->GetYaxis()->SetTitle("Onia 2 Invariant Mass [GeV]");

  ANA_CHECK( book( TH2D( "noUpsi_1v2SS", "SS Onia 1 vs Onia 2 Invariant Mass No Upsilon(1S)", 140, 1, 15, 140, 1, 15)));
  hist ("noUpsi_1v2SS")->GetXaxis()->SetTitle("Onia 1 Invariant Mass [GeV]");
  hist ("noUpsi_1v2SS")->GetYaxis()->SetTitle("Onia 2 Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "OS_noUp_M", "OS Onia Invariant Mass No Upsilon", 140, 1, 15)));
  hist ("OS_noUp_M")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "Upsi_pairM", "Upsilon Pairs Invariant Mass", 140, 1, 15)));
  hist ("Upsi_pairM")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "SS_noUp_M", "SS Onia Invariant Mass No Upsilon", 140, 1, 15)));
  hist ("SS_noUp_M")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "X_oneUpsi", "X (Quad) Invariant Mass Onia + One Upsilon(1S) Candidate", 60, 15, 25)));
  hist ("X_oneUpsi")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "X_twoUpsi", "X (Quad) Invariant Mass Onia + Two Upsilon(1S) Candidates", 60, 15, 25)));
  hist ("X_twoUpsi")->GetXaxis()->SetTitle("Invariant Mass [GeV]");

  ANA_CHECK( book( TH1D( "dimu_oneUpsi", "Onia Pair Invariant Mass One Upsilon(1S)", 60, 0, 15)));
  hist ("dimu_oneUpsi")->GetXaxis()->SetTitle("Invariant Mass [GeV]");


  //Opening files to read/write event-run numbers
  /*
  Dren665.open("/data/bhollo/SampleT2022/runevntDAOD507665.txt");
  Dren666.open("/data/bhollo/SampleT2022/runevntDAOD507666.txt");   
  Dren8.open("/data/bhollo/SampleT2022/runevntDAOD800852.txt");
  Dren8.open("/data/bhollo/SampleT2022/800852_list_dup.txt");
  Dren8.open("/data/bhollo/SampleT2022/muHLT_852_10%Pass.txt");
  Dset8.open("/data/bhollo/SampleT2022/runevntDAOD800852.txt");
  */
  M11evnt_list.open("/data/bhollo/SampleT2022/M11evnt_list.txt");
  M14evnt_list.open("/data/bhollo/SampleT2022/M14evnt_list.txt");
  mu22evnt_list.open("/data/bhollo/SampleT2022/mu22evnt_list.txt");
  mu20evnt_list.open("/data/bhollo/SampleT2022/mu20m8evnt_list.txt");
  //M14L1evnt_list.open("/data/bhollo/SampleT2022/M14L1evnt_list.txt");
  
 
  //split each string line and insert it into the run_event_set map
  /* 
  string line;
  while ( getline(Dset8, line) )
    {
      vector<string> v;
      boost::split( v, line, boost::is_any_of(" "));
      
      pair<unsigned int, unsigned long long> numSet;
      numSet.first = stoi(v[0]);
      numSet.second = stoll(v[1]);
                  
      run_event_set8.insert( numSet ); //insert run #-event # set into the run_event_set map
      //ANA_MSG_INFO( "run-event number inserted");
    }
  
  //print the pairs in the run event set  
  for (auto const &var: run_event_set)
    {
      cout << "(" << var.first << "," << var.second << ")" << endl;
    }
  
  string line22;
  while ( getline(mu22evnt_list, line22) )
    {
      vector<string> v;
      boost::split( v, line22, boost::is_any_of(" "));

      pair<unsigned int, unsigned long long> numSet;
      numSet.first = stoi(v[0]);
      numSet.second = stoll(v[1]);

      run_event_set22.insert( numSet ); //insert run #-event # set into the run_event_set map
    }
  
  string line14;
  while ( getline(M14evnt_list, line14) )
    {
      vector<string> v;
      boost::split( v, line14, boost::is_any_of(" "));

      pair<unsigned int, unsigned long long> numSet;
      numSet.first = stoi(v[0]);
      numSet.second = stoll(v[1]);

      run_event_setM14.insert( numSet ); //insert run #-event # set into the run_event_set map
    }

  string line20;
  while ( getline(mu20evnt_list, line20) )
    {
      vector<string> v;
      boost::split( v, line20, boost::is_any_of(" "));

      pair<unsigned int, unsigned long long> numSet;
      numSet.first = stoi(v[0]);
      numSet.second = stoll(v[1]);

      run_event_set20.insert( numSet ); //insert run #-event # set into the run_event_set map
    }

  string line11;
  while ( getline(M11evnt_list, line11) )
    {
      vector<string> v;
      boost::split( v, line11, boost::is_any_of(" "));

      pair<unsigned int, unsigned long long> numSet;
      numSet.first = stoi(v[0]);
      numSet.second = stoll(v[1]);

      run_event_setM11.insert( numSet ); //insert run #-event # set into the run_event_set map
    }
  */
  
  return StatusCode::SUCCESS;
}



StatusCode MyxAODAnalysis :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code ll go.
  
  const xAOD::EventInfo *evntInfo = nullptr;
  ANA_CHECK( evtStore()->retrieve( evntInfo, "EventInfo" ) );

  unsigned int runNum = evntInfo->runNumber();
  runlist.push_back(runNum);

  unsigned long long evntNum = evntInfo->eventNumber();
  eventlist.push_back(evntNum);
  
  const xAOD::VertexContainer *quads = nullptr;
  ANA_CHECK( evtStore()->retrieve( quads, "BPHY4Quads" ));
  
  vector<xAOD::Vertex> Quads;
  for (auto v: *quads)
    {
    Quads.push_back(*v);
    }
  
  sort( Quads.begin(), Quads.end(), [](auto q, auto p) -> bool { return q.chiSquared() < p.chiSquared(); } ); //sort quads by chiSquared

  const xAOD::VertexContainer *pairs = nullptr;
  ANA_CHECK( evtStore()->retrieve( pairs, "BPHY4Pairs" ));

  /*
  //Matching run-event numbers from AOD and xAOD data sets
  pair<unsigned int, unsigned long long> goodSet20;
  goodSet20.first = runNum;
  goodSet20.second = evntNum;

  pair<unsigned int, unsigned long long> goodSet22;
  goodSet22.first = runNum;
  goodSet22.second = evntNum;

  pair<unsigned int, unsigned long long> goodSet11;
  goodSet11.first = runNum;
  goodSet11.second = evntNum;
  
  pair<unsigned int, unsigned long long> goodSet14;
  goodSet14.first = runNum;
  goodSet14.second = evntNum;
  */

  m_evtNr += 1;

  if ( Quads.size() > 0 && pairs->size() > 0 )
    {
      //good 4 muons selection
      vector muonsQ = (Quads[0]).auxdata< vector<ElementLink<xAOD::MuonContainer>>  >("MuonLinks"); //muons directly related to quad as vector link objects
      float quadMass = (Quads[0]).auxdata<float>("QUAD_mass");
          			    
      for (auto mu: muonsQ)
	{
	  double mu_pt = (*mu)->pt();
	  hist ("quadMupT")->Fill( (mu_pt)/1000 );
	  
	  double mu_eta = (*mu)->eta();
          hist ("Xmu_eta")->Fill( mu_eta );

	  if ( mu_pt/1000 <= 3 )
	    {
	      hist ("muEta_low3")->Fill( mu_eta );
	    }

	  if ( mu_pt/1000 > 3 )
            {
              hist ("muEta_up3")->Fill( mu_eta );
            }	  
	}

      const xAOD::Muon* quadMu1 = *muonsQ[0];
      const xAOD::Muon* quadMu2 = *muonsQ[1];
      const xAOD::Muon* quadMu3 = *muonsQ[2];
      const xAOD::Muon* quadMu4 = *muonsQ[3];

      const xAOD::TrackParticle* mu1_IDt = *( quadMu1->inDetTrackParticleLink() );
      const xAOD::TrackParticle* mu2_IDt = *( quadMu2->inDetTrackParticleLink() );
      const xAOD::TrackParticle* mu3_IDt = *( quadMu3->inDetTrackParticleLink() );
      const xAOD::TrackParticle* mu4_IDt = *( quadMu4->inDetTrackParticleLink() );
                  
      double quad_pT = ( quadMu1->p4() + quadMu2->p4() + quadMu3->p4() + quadMu4->p4() ).Pt();
      double X_M = ( quadMu1->p4() + quadMu2->p4() + quadMu3->p4() + quadMu4->p4() ).M();
      double X_charge = quadMu1->charge() + quadMu2->charge() + quadMu3->charge() + quadMu4->charge();

      hist( "X_invMass" )->Fill( X_M/1000 );
      hist( "X_Mdiff" )->Fill( (quadMass-X_M)/1000 );
      hist ( "X_mResolP")->Fill( X_M/1000, (quadMass-X_M)/1000 );
      hist ( "X_mResolD")->Fill( quadMass/1000, (quadMass-X_M)/1000 );

      vector<double> quadMu_pt;
      quadMu_pt.push_back(mu1_IDt->pt()/1000);
      quadMu_pt.push_back(mu2_IDt->pt()/1000);
      quadMu_pt.push_back(mu3_IDt->pt()/1000);
      quadMu_pt.push_back(mu4_IDt->pt()/1000);

      vector<int> quadMu_type;
      quadMu_type.push_back( quadMu1->muonType() );
      quadMu_type.push_back( quadMu2->muonType() );
      quadMu_type.push_back( quadMu3->muonType() );
      quadMu_type.push_back( quadMu4->muonType() );

      sort( quadMu_pt.begin(), quadMu_pt.end(), [](auto q, auto p) -> bool { return q < p; } );
      hist( "lowestMupT" )->Fill( quadMu_pt[0] );
      hist( "SecLowMu_pT")->Fill( quadMu_pt[1] );
      hist( "SecHighMu_pT")->Fill( quadMu_pt[2] );
      hist( "highestMu_pT")->Fill( quadMu_pt[3] );

      int ComMu = 0;  //targeting combined muons, mu->muonType() == 0
      int totCom = count( quadMu_type.begin(), quadMu_type.end(), ComMu ) ;//count how many muons in the quad are combined muons 
      if ( totCom >= 3 ) //at least 3 of the 4 muons are combined muons
	{
	  hist ("M_ComCutO")->Fill( (X_M)/1000 );
	}

      if ( quad_pT/1000 < 10 )
	{
	  hist( "lowestMupT_X" )->Fill( quadMu_pt[0] );
	}
      
      if ( mu1_IDt->pt()/1000 > 3 and mu2_IDt->pt()/1000 > 3 and mu3_IDt->pt()/1000 > 3 and mu4_IDt->pt()/1000 > 3)
        {
          hist ("M_ptCut")->Fill( (X_M)/1000 );

	  if ( quadMu_pt[2] > 4 and quadMu_pt[3] > 4 )
	    {
	      hist ("M_ptCut2")->Fill( (X_M)/1000 ); //all mu_pT > 3 and at least 2 mu_pT > 4
	      
	      int target = 0; //targeting combined muons, mu->muonType() == 0
	      int res = count( quadMu_type.begin(), quadMu_type.end(), target ); //count how many muons in the quad are combined muons

	      if ( res >= 3 ) //at least 3 of the 4 muons are combined muons
		{
		  hist ("M_ComCut")->Fill( (X_M)/1000 );
		  hist ("M_gQuadCut")->Fill( (X_M)/1000 );		  
		  //dimuon pair selections
		  vector<xAOD::Vertex> quadPairs;
		  vector<xAOD::Vertex> Onia_MPairs;
		  vector<xAOD::Vertex> oniPairs;
		  vector<xAOD::Vertex> upsiPairs;
		  vector<xAOD::Vertex> OS_oniPairs;
		  vector<xAOD::Vertex> OS_oniPairs_noUpsi;
		  vector<xAOD::Vertex> SS_oniPairs_noUpsi;

		  for ( auto p: *pairs)
		    {
		      vector muonsP = p->auxdata< vector<ElementLink<xAOD::MuonContainer>> >("MuonLinks"); //muons linked to doublets

		      if ( count( muonsQ.begin(), muonsQ.end(), muonsP[0] ) && count( muonsQ.begin(), muonsQ.end(), muonsP[1] ) ) //match pair to quad
			{
			  quadPairs.push_back(*p);
			}
		    }
		  
		  for ( auto  qP: quadPairs )
		    {
		      vector muonsQP = qP.auxdata< vector<ElementLink<xAOD::MuonContainer>> >("MuonLinks");
                      double qP_M = ( (*muonsQP[0])->p4() + (*muonsQP[1])->p4() ).M();
		      
		      hist ("dimuIM")->Fill( qP_M/1000 ); //pair invaraint mass

		      if ( qP_M/1000 >= 2 && qP_M/1000 <= 50 ) //eliminate pairs not in mass range
			{
			  Onia_MPairs.push_back(qP);
			  hist ("dimuChi")->Fill( qP.chiSquared() ); //chi squared of quad pairs
			  if ( qP.chiSquared() > 10 )
			    {
			      //ANA_MSG_INFO( qP.chiSquared() ); //print chi squared values
			      chi_1K += 1;
			    }
			  if ( qP.chiSquared() < 3 ) //remove pairs
			    {
			      oniPairs.push_back(qP);
			    }
			}
		      else
			{
			  hist ("M_noOni")->Fill( X_M/1000 );
			}
		    }

		  if ( Onia_MPairs.size() > 0 )
		    {
		      hist ("M_2muIM")->Fill( X_M/1000 );
		    }
		  
		  if ( oniPairs.size() > 0 )
		    {
		      hist ("M_2muChi")->Fill( X_M/1000 );
		      hist ("M_OniaCut")->Fill( X_M/1000 );
		    }
		  
		  sort ( oniPairs.begin(), oniPairs.end(), [](auto q, auto p)-> bool{ return q.chiSquared() < p.chiSquared(); } );

		  for ( auto oP: oniPairs )
		    {
		      vector muonsUp = oP.auxdata< vector<ElementLink<xAOD::MuonContainer>> >("MuonLinks");
		      const xAOD::Muon* mu1P = *muonsUp[0];
		      const xAOD::Muon* mu2P = *muonsUp[1];
		      double oP_M = ( mu1P->p4() + mu2P->p4() ).M();

		      if (  mu1P->charge() != mu2P->charge() )
			{
			  hist ("dimuIM_OS")->Fill( oP_M/1000 );
			  OS_oniPairs.push_back(oP);
			  
			  if ( oP_M/1000 >= 9.2 && oP_M/1000 <= 9.7 )
			    {
			      upsiPairs.push_back(oP);
			    }
			  else
			    {
			      OS_oniPairs_noUpsi.push_back(oP);
			    }
			}
		      else
			{
			  SS_oniPairs_noUpsi.push_back(oP);
			}
		    }
		  
		  for ( auto Up: upsiPairs )
		    {
		      vector muons_Up = Up.auxdata< vector<ElementLink<xAOD::MuonContainer>> >("MuonLinks");
		      double Up_M = ( (*muons_Up[0])->p4() + (*muons_Up[1])->p4() ).M();
		      hist ("Upsi_pairM")->Fill( Up_M/1000 );
		    }
	       
		  for ( auto OS_noUp: OS_oniPairs_noUpsi )
		    {
		      vector OS_noUp_mu = OS_noUp.auxdata< vector<ElementLink<xAOD::MuonContainer>> >("MuonLinks");
		      double OS_noUp_M = ( (*OS_noUp_mu[0])->p4() + (*OS_noUp_mu[1])->p4() ).M();
		      hist ("OS_noUp_M")->Fill( OS_noUp_M/1000 );
		    }

		  if ( OS_oniPairs_noUpsi.size() >= 2 )
		    {
		      vector muons_Oni1 = (OS_oniPairs_noUpsi[0]).auxdata< vector<ElementLink<xAOD::MuonContainer>> >("MuonLinks");
		      vector muons_Oni2 = (OS_oniPairs_noUpsi[1]).auxdata< vector<ElementLink<xAOD::MuonContainer>> >("MuonLinks");
		      double oni1_M = ( (*muons_Oni1[0])->p4() + (*muons_Oni1[1])->p4() ).M();
		      double oni2_M = ( (*muons_Oni2[0])->p4() + (*muons_Oni2[1])->p4() ).M();
		      hist ("noUpsi_1v2OS")->Fill( oni1_M/1000, oni2_M/1000 );
		    }

		  for ( auto SS_noUp: SS_oniPairs_noUpsi )
		    {
		      vector SS_noUp_mu = SS_noUp.auxdata< vector<ElementLink<xAOD::MuonContainer>> >("MuonLinks");
		      double SS_noUp_M = ( (*SS_noUp_mu[0])->p4() + (*SS_noUp_mu[1])->p4() ).M();
		      hist ("SS_noUp_M")->Fill( SS_noUp_M/1000 );
		    }

		  if ( SS_oniPairs_noUpsi.size() >= 2 )
		    {
		      vector muons_Oni1 = (SS_oniPairs_noUpsi[0]).auxdata< vector<ElementLink<xAOD::MuonContainer>> >("MuonLinks");
		      vector muons_Oni2 = (SS_oniPairs_noUpsi[1]).auxdata< vector<ElementLink<xAOD::MuonContainer>> >("MuonLinks");
		      double oni1_M = ( (*muons_Oni1[0])->p4() + (*muons_Oni1[1])->p4() ).M();
		      double oni2_M = ( (*muons_Oni2[0])->p4() + (*muons_Oni2[1])->p4() ).M();
		      hist ("noUpsi_1v2SS")->Fill( oni1_M/1000, oni2_M/1000 );
		    }
		  
		  if ( upsiPairs.size() < 1 )
		    {
		      noUpsi += 1;
		      hist ("X_Onia_noUpsi")->Fill( X_M/1000 );
		      //plotting mass of pairs for no upsilon events
		      for ( auto oP_noUp: oniPairs )
			{
			  vector muons_noUp = oP_noUp.auxdata< vector<ElementLink<xAOD::MuonContainer>> >("MuonLinks");
			  double oPnoUp_M = ( (*muons_noUp[0])->p4() + (*muons_noUp[1])->p4() ).M();
			  hist ("dimu_noUpsi")->Fill( oPnoUp_M/1000 );
			  if ( oPnoUp_M/1000 < 9.2 )
			    {
			      hist ("X_noUpsiL")->Fill( X_M/1000 );
			    }
			  if ( oPnoUp_M/1000 > 9.7 )
			    {
			      hist ("X_noUpsiH")->Fill( X_M/1000 );
			    }
			}
		    }
		  
		  if ( upsiPairs.size() == 1 )
                    {
                      oneUpsi += 1;
		      //ANA_MSG_INFO("There are " << OS_oniPairs.size() << " in this event.");
		      
		      hist ("X_oneUpsi")->Fill( X_M/1000 );
		      for ( auto oneUp: oniPairs )
			{
			  vector muons_oneUp = oneUp.auxdata< vector<ElementLink<xAOD::MuonContainer>> >("MuonLinks");
			  double oneUPair_M = ( (*muons_oneUp[0])->p4() + (*muons_oneUp[1])->p4() ).M();
			  hist ("dimu_oneUpsi")->Fill( oneUPair_M/1000 );
			}
		    }
		  
		  if ( upsiPairs.size() == 2 )
                    {
                      twoUpsi += 1;
		      hist ("X_twoUpsi")->Fill( X_M/1000 );
                    }
		  		  		  
		  if ( upsiPairs.size() > 0 )
		    {
		      hist ("M_2muUpsi")->Fill( X_M/1000 );
		      if ( X_charge == 0 )
			{
			  //ANA_MSG_INFO( "This event has a 0 charge event");
			  hist ("M_0Charge")->Fill( X_M/1000 );
			  if ( quadMass/1000 < 50 )
			    {
			      hist ("M_finalIM")->Fill( X_M/1000 );
			    }
			  
			  if ( upsiPairs.size() == 1 )
			    {
			      hist ("M_0Charge_1upsi")->Fill( X_M/1000 );
			    }

			  if ( upsiPairs.size() == 2 )
			    {
			      hist ("M_0Charge_2upsi")->Fill( X_M/1000 );
			    }
			}

		      if ( X_charge == abs(2) )
			{
			  hist ("M_2Charge")->Fill( X_M/1000 ); 
			}
		    }
		}
	    }
	}

      if ( mu1_IDt->eta() < abs(2.5) and  mu2_IDt->eta() < abs(2.5) and  mu3_IDt->eta() < abs(2.5) and  mu4_IDt->eta() < abs(2.5) )
	{
	  hist ("M_etaCut")->Fill( (X_M)/1000 );
	}      
		     
      //hist ("total_pT")->Fill( (quad_pT)/1000 );
      /*    
      if (  ((Quads[0]).auxdata< float >("QUAD_mass") * 0.001) < 19 && ((Quads[0]).auxdata< float >("QUAD_mass") * 0.001) > 17 ) //mass cut on best chiSquared quad to be 18 GeV
	{
	  Dren666 << runNum << " " << evntNum << " " <<  endl;
	}
      
      if (run_event_set20.count(goodSet20) == 1)
	{
	  hist ("passMu20_pT")->Fill( (quad_pT)/1000 );
	}
      
      if (run_event_setM14.count(goodSet14) == 1)
	{
	  hist ("passM14_pT")->Fill( (quad_pT)/1000 );
	}

      if (run_event_setM11.count(goodSet11) == 1)
        {
          hist ("passM11_pT")->Fill( (quad_pT)/1000 );
        }

      if (run_event_set22.count(goodSet22) == 1)
        {
          hist ("passMu22_pT")->Fill( (quad_pT)/1000 );
        }
      */
      
    }
  //Write run and event numbers for events that passed specfic trigger  
  /*
  if ( m_7M14_L1->isPassed() )
    {
      M14L1evnt_list << runNum << " " << evntNum << endl;
    }
  
  if ( m_mu20_mu8noL1->isPassed() )
    {
      mu20evnt_list << runNum << " " << evntNum << endl;
    }
  
  pair<unsigned int, unsigned long long> goodSet;
  goodSet.first = runNum;
  goodSet.second = evntNum;

  for (const auto& trigg: m_allHLT_mu->getListOfTriggers() ) //looping over HLT-mu trigger list
    {
      if (run_event_set8.count(goodSet) == 1)
	{
	  //ANA_MSG_INFO( "matching run-event number set going forward");

	  if (m_trigDec->isPassed(trigg) == 1 )
	    {
	      hist ("trigFreq_800852")->Fill ( trigg.c_str(), 1 ); //fill histogram with string trigg in the form of a c-string (char string) 
	      //Dren8 << trigg.c_str() << endl; //prints entire c-string
	    }  
	}  
    }
      
  //obtaining features and dimuon combinations for di-muon triggers
  FeatureRequestDescriptor featureRequestDescriptor_M14("HLT_2mu4_L1BPH-7M14-0DR25-MU5VFMU3VF", TrigDefs::Physics, "", TrigDefs::lastFeatureOfType, "feature", -1);
  vector < TrigCompositeUtils::LinkInfo<xAOD::IParticleContainer> > features_M14 = m_trigDec->features<xAOD::IParticleContainer>(featureRequestDescriptor_M14);
  
  FeatureRequestDescriptor featureRequestDescriptor_M11("HLT_2mu4_L1BPH-7M11-25DR99-2MU3VF", TrigDefs::Physics, "", TrigDefs::lastFeatureOfType, "feature", -1);
  vector < TrigCompositeUtils::LinkInfo<xAOD::IParticleContainer> > features_M11 = m_trigDec->features<xAOD::IParticleContainer>(featureRequestDescriptor_M11);

  const TrigConf::HLTChain* chainInfo_M14 = m_trigDec->ExperimentalAndExpertMethods().getChainConfigurationDetails("HLT_2mu4_L1BPH-7M14-0DR25-MU5VFMU3VF");
  const TrigConf::HLTChain* chainInfo_M11 = m_trigDec->ExperimentalAndExpertMethods().getChainConfigurationDetails("HLT_2mu4_L1BPH-7M11-25DR99-2MU3VF");

  TrigCompositeUtils::Combinations dimu_combos14 = TrigCompositeUtils::buildCombinations("HLT_2mu4_L1BPH-7M14-0DR25-MU5VFMU3VF", features_M14, chainInfo_M14, TrigCompositeUtils::FilterType::UniqueObjects);
  TrigCompositeUtils::Combinations dimu_combos11 = TrigCompositeUtils::buildCombinations("HLT_2mu4_L1BPH-7M11-25DR99-2MU3VF", features_M11, chainInfo_M11, TrigCompositeUtils::FilterType::UniqueObjects); 
  
  bool pass_OS11 = false;
  bool pass_OS_mass11 = false;
  bool pass_OS14 = false;
  bool pass_OS_mass14 = false;

  //iterate over each combination in the di muon combinations object for 7M11 trigger
  //size_t combo_n = 0; //define the combination number in the dimu_combos
  for (const auto& combo : dimu_combos11)
    {
      //size_t mu_n = 0; //define the muon number in each combination
      
      for (const auto& lInfo : combo) 
	{
	  const xAOD::IParticle* mu = *(lInfo.link); //identify each individual muon; identical to **linkInfo.link
	  ANA_MSG_INFO( "Muon " << mu_n++ << " from combination " << combo_n << " has charge  = " << dynamic_cast<const xAOD::Muon*>(mu)->charge() );
	}
      
      const xAOD::Muon* mu1 = dynamic_cast<const xAOD::Muon*>( *(combo[0].link) );
      const xAOD::Muon* mu2 = dynamic_cast<const xAOD::Muon*>( *(combo[1].link) ); 
                  
      TLorentzVector mu1_p4 = mu1->p4();
      TLorentzVector mu2_p4 = mu2->p4();
            
      const double dimu_IM = (mu1_p4 + mu2_p4).M();
      
      if ( mu1->charge() != mu2->charge() )
	{
	  pass_OS11 = true;
	}

      if ( mu1->charge() != mu2->charge() && (dimu_IM/1000) >= 8.5 && (dimu_IM/1000) <= 10 )
        {
          pass_OS_mass11 = true;
	}

      if ( pass_OS11 )
	{
	  hist ("dimu_M11")->Fill( (dimu_IM)/1000 );
	  //break;
	}

      if ( pass_OS_mass11 )
	{
	  hist ("dimu_M11_OSU")->Fill( (dimu_IM)/1000 );
	  // break;
        }
    }
  
  for (const auto& combo: dimu_combos14)
    {
      const xAOD::Muon* mu1 = dynamic_cast<const xAOD::Muon*>( *(combo[0].link) );
      const xAOD::Muon* mu2 = dynamic_cast<const xAOD::Muon*>( *(combo[1].link) );

      TLorentzVector mu1_p4 = mu1->p4();
      TLorentzVector mu2_p4 = mu2->p4();
      
      const double dimu_IM = (mu1_p4 + mu2_p4).M();
            
      if ( mu1->charge() != mu2->charge() )
	{
	  pass_OS14 = true;
	  hist ("dimu_M14")->Fill( (dimu_IM)/1000 );
	  //break;
	}

      if ( mu1->charge() != mu2->charge() && (dimu_IM/1000) >= 8.5 && (dimu_IM/1000) <= 10 )
        {
          pass_OS_mass14 = true;
	  hist ("dimu_M14_OSU")->Fill( (dimu_IM)/1000 );
        }
    }

  if (pass_OS14)
    {
      M14_evntNum += 1;
    }

  if (pass_OS_mass14)
    {
      M11_evntNum += 1;
    }
  */

   return StatusCode::SUCCESS;
}



StatusCode MyxAODAnalysis :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they a re written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.
  /*
  //loop over x bins of histogram
  for (auto x = 0; x < 408; x++)
    {
      auto eventCnt = hist ("trigFreq_800852")->GetBinContent(x); //obtain number of entries for each bin
	  
      Dren8 << eventCnt << " " << " " << hist ("trigFreq_800852")->GetXaxis()->GetBinLabel(x) << endl;
    }
  
  TFile pFile = TFile( "trig_Eff.root", "recreate" );   
  
  if (TEfficiency::CheckConsistency( *hist("passM11_pT"), *hist("total_pT") ))
    {
      pEff_M11 = new TEfficiency( *hist("passM11_pT"), *hist("total_pT") );
      pEff_M11->SetTitle("HLT_M11 Quad pT Efficiency");
      pFile.WriteTObject(pEff_M11, "pEff_M11");
    }
  
  if (TEfficiency::CheckConsistency( *hist("passM14_pT"), *hist("total_pT") ))
    {
      pEff_M14 = new TEfficiency( *hist("passM14_pT"), *hist("total_pT") );
      pEff_M14->SetTitle("HLT_M14 Quad pT Efficiency");
      pFile.WriteTObject(pEff_M14, "pEff_M14");
    }
  
 if (TEfficiency::CheckConsistency( *hist("passMu20_pT"), *hist("total_pT") ))
    {
      pEff_mu20 = new TEfficiency( *hist("passMu20_pT"), *hist("total_pT") );
      pEff_mu20->SetTitle("HLT_mu20_mu8noL1 Quad pT Efficiency");
      pFile.WriteTObject(pEff_mu20, "pEff_mu20");
    }

 if (TEfficiency::CheckConsistency( *hist("passMu22_pT"), *hist("total_pT") ))
    {
      pEff_mu22 = new TEfficiency( *hist("passMu22_pT"), *hist("total_pT") );
      pEff_mu22->SetTitle("HLT_mu22_mu8noL1 Quad pT Efficiency");
      pFile.WriteTObject(pEff_mu22, "pEff_mu22");
    }
 
  if (TEfficiency::CheckConsistency( *hist("passM14L1_pT"), *hist("total_pT") ))
    {
      pEff_M14L1 = new TEfficiency( *hist("passM14L1_pT"), *hist("total_pT") );
      pEff_M14L1->SetTitle("L1_M14 Quad pT Efficiency");
      pFile.WriteTObject(pEff_M14L1, "pEff_M14L1");
    }
 
  pFile.Close();
  
  Dren665.close();
  Dren666.close();
  Dren8.close();
  Dset8.close();
  */

  ANA_MSG_INFO( "There are " << m_evtNr << " total events" );
  ANA_MSG_INFO("There are " << chi_1K << " pairs with chi squared > 100.");
  //ANA_MSG_INFO( "There are " << noUpsi << " events with no Upsilon candidate");
  //ANA_MSG_INFO( "There are " << oneUpsi << " events with one Upsilon candidate");
  //ANA_MSG_INFO( "There are " << twoUpsi << " events with two Upsilon candidates");
  M11evnt_list.close();
  M14evnt_list.close();
  mu22evnt_list.close();
  mu20evnt_list.close();
  M14L1evnt_list.close();
  
  hist ("noUpsi_1v2OS")->SetStats(0);
  hist ("noUpsi_1v2SS")->SetStats(0);
  
  return StatusCode::SUCCESS;
}
