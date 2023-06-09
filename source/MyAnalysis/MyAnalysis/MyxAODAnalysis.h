#ifndef MyAnalysis_MyxAODAnalysis_H
#define MyAnalysis_MyxAODAnalysis_H

#include <AnaAlgorithm/AnaAlgorithm.h>
#include <AsgTools/ToolHandle.h>
#include <TrigDecisionTool/TrigDecisionTool.h>
#include <TrigCompositeUtils/TrigCompositeUtils.h>

//#include <EventInfo/EventInfo.h>

#include <TEfficiency.h>
#include <fstream>
#include <set>

using namespace Trig;

class MyxAODAnalysis : public EL::AnaAlgorithm
{
 public:
  // this is a standard algorithm constructor
  MyxAODAnalysis (const std::string& name, ISvcLocator* pSvcLocator);

  // these are the functions inherited from Algorithm
  virtual StatusCode initialize () override;
  virtual StatusCode execute () override;
  virtual StatusCode finalize () override;

 private:
  // Configuration, and any other types of variables go here.
  //float m_cutValue;
  //TTree *m_myTree;
  //TH1 *m_myHist;

  StatusCode checkTriggerDecision(); //!< decision bit analysis
  StatusCode checkLevels(); //!< checks if levels passed are right
  StatusCode printLevels();
  StatusCode printChainConfiguration();

  ToolHandle<Trig::TrigDecisionTool> m_trigDec; //!< TDT handle
  
  //Counts of events for specfics
  int m_evtNr;
  int M11_evntNum;
  int M14_evntNum;
  int noUpsi;
  int oneUpsi;
  int twoUpsi;

  //Counts of objects with specfics
  int chi_1K;

  const ChainGroup* m_all;
  const ChainGroup* m_allL1;
  const ChainGroup* m_allHLT;
  const ChainGroup* m_allHLT_mu;
  const ChainGroup* m_2mu4_7M14;
  const ChainGroup* m_2mu4_7M11;
  const ChainGroup* m_mu22_mu8noL1;
  const ChainGroup* m_mu20_mu8noL1;
  const ChainGroup* m_7M14_L1;

  std::vector< std::string > m_chain_names;
  std::vector< std::string > m_cfg_chains;

  std::vector<unsigned int> runlist;
  std::vector<unsigned long long>  eventlist;
  std::ofstream AODren;
  std::ofstream Dren8;
  std::ofstream Dren665;
  std::ofstream Dren666;
  std::ifstream Dset8;
  std::ifstream Dset665;
  std::ifstream Dset666;

  std::fstream M11evnt_list;
  std::fstream M14evnt_list;
  std::fstream mu22evnt_list;
  std::fstream mu20evnt_list;
  std::fstream M14L1evnt_list;
  
  std::set<std::pair<unsigned int, unsigned long long>> run_event_set8;
  std::set<std::pair<unsigned int, unsigned long long>> run_event_set665;
  std::set<std::pair<unsigned int, unsigned long long>> run_event_set666;
  
  std::set<std::pair<unsigned int, unsigned long long>> run_event_setM11;
  std::set<std::pair<unsigned int, unsigned long long>> run_event_setM14;
  std::set<std::pair<unsigned int, unsigned long long>> run_event_set22;
  std::set<std::pair<unsigned int, unsigned long long>> run_event_set20;
  std::set<std::pair<unsigned int, unsigned long long>> run_event_setM14L1;

  TEfficiency* pEff_M11;
  TEfficiency* pEff_M14;
  TEfficiency* pEff_mu20;
  TEfficiency* pEff_mu22;
  TEfficiency* pEff_M14L1;
};

#endif

