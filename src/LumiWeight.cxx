#include "UHH2/UHHNtupleConverter/include/LumiWeight.h"

#include <iostream>

//==============================================================================================
LumiWeight::LumiWeight( void ){

}

//==============================================================================================
LumiWeight::~LumiWeight( void ){

}

//==============================================================================================
double LumiWeight::getLumiWeight( TString sample ){

  if( sample.Contains("BulkGrav") or sample.Contains("Qstar") or sample.Contains("Wprime") or sample.Contains("Zprime") or sample.Contains("Radion")) return 1.0;

  if( sample.Contains( "WW_"      ) ) return 118.7; //NNLO, https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
  if( sample.Contains( "WWTo4Q"      ) ) return 53.94; //NNLO, https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
  if( sample.Contains( "WZ_"      ) ) return 47.13; //NLO, https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson 
  if( sample.Contains( "ZZ_"      ) ) return 16.5; //NLO, https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
  if( sample.Contains( "ZZTo4Q"      ) ) return 8.06; //NLO, https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
  if( sample.Contains( "WminusH_HToBB_WToQQ_M125"      ) ) return 0.533*0.5824*0.6741;
  if( sample.Contains( "WplusH_HToBB_WToQQ_M125"      ) ) return 0.84*0.5824*0.6741;
  if( sample.Contains( "ZH_HToBB_ZToQQ_M125"      ) ) return 0.8839*0.5824*0.69911;
  if( sample.Contains( "TT_Mtt-700to1000"      ) ) return 76.605;
  if( sample.Contains( "TT_Mtt-1000toInf"      ) ) return 20.578;
  if( sample.Contains( "TTToHadronic"          ) ) return 380.13;
  //if( sample.Contains( "WJetsToQQ_HT"          ) ) return 95.14;
  //if( sample.Contains( "ZJetsToQQ_HT600toInf"  ) ) return 41.34;
  if( sample.Contains( "WJetsToQQ_HT400to600"          ) ) return 315.6;//EXO-18-012, http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2017/335
  if( sample.Contains( "WJetsToQQ_HT600to800"          ) ) return 68.57;//EXO-18-012, http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2017/335
  if( sample.Contains( "WJetsToQQ_HT-800toInf"          ) ) return 34.9;//EXO-18-012, http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2017/335
  if( sample.Contains( "ZJetsToQQ_HT400to600"          ) ) return 145.4;//EXO-18-012, http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2017/335
  if( sample.Contains( "ZJetsToQQ_HT600to800"          ) ) return 34.0;//EXO-18-012, http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2017/335
  if( sample.Contains( "ZJetsToQQ_HT-800toInf"          ) ) return 18.67;//EXO-18-012, http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2017/335
  if( sample.Contains( "QCD_Pt_170to300_"      ) ) return 117276.;
  if( sample.Contains( "QCD_Pt_300to470_"      ) ) return 7823.;                   
  if( sample.Contains( "QCD_Pt_470to600_"      ) ) return 648.2;                  
  if( sample.Contains( "QCD_Pt_600to800_"      ) ) return 186.9;                 
  if( sample.Contains( "QCD_Pt_800to1000_"     ) ) return 32.293;               
  if( sample.Contains( "QCD_Pt_1000to1400_"    ) ) return 9.4183; 
  if( sample.Contains( "QCD_Pt_1400to1800_"    ) ) return 0.84265; 
  if( sample.Contains( "QCD_Pt_1800to2400_"    ) ) return 0.114943; 
  if( sample.Contains( "QCD_Pt_2400to3200_"    ) ) return 0.006830; 
  if( sample.Contains( "QCD_Pt_3200toInf_"     ) ) return 0.000165445; 
  if( sample.Contains( "QCD_HT100to200"        ) ) return 2.785e07;
  if( sample.Contains( "QCD_HT200to300"        ) ) return 1717000.;
  if( sample.Contains( "QCD_HT300to500"        ) ) return 351300.0;
  if( sample.Contains( "QCD_HT500to700"        ) ) return 31630.;
  if( sample.Contains( "QCD_HT700to1000"       ) ) return 6802.;
  if( sample.Contains( "QCD_HT1000to1500"      ) ) return 1206.0;
  if( sample.Contains( "QCD_HT1500to2000"      ) ) return 120.4;
  if( sample.Contains( "QCD_HT2000toInf"       ) ) return 25.25;
  if( sample.Contains( "QCD_Pt-15to7000" ) or sample.Contains( "QCD_Pt_15to7000" )) return  2.022100000e+09*60.5387252324;
  
  
  std::cout << "Cross section not defined for this sample!!" << std::endl;
  return 0;
  
}

















