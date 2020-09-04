#include "UHH2/UHHNtupleConverter/include/UHHNtupleConverterSelections.h"

#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Particle.h"
#include "UHH2/common/include/Utils.h"

#include <stdexcept>
#include <iostream>

using namespace uhh2examples;
using namespace uhh2;
using namespace std;
bool PRINT = false;

DijetSelection::DijetSelection(float deta_max_, float mjj_min_): deta_max(deta_max_), mjj_min(mjj_min_){}
    
bool DijetSelection::passes(const Event & event){
    assert(event.jets); // if this fails, it probably means jets are not read in
    if(event.jets->size() < 2) return false;
    const auto & jet0 = event.jets->at(0);
    const auto & jet1 = event.jets->at(1);
    auto deta = fabs(jet0.eta() - jet1.eta());
    if(deta > deta_max) return false;
    if(inv_mass_safe(jet0.v4()+jet1.v4()) < mjj_min ) return false;
    return true;
}

MuonVeto::MuonVeto(const MuonId & muid_, float deltaR_min_): deltaR_min(deltaR_min_), muid(muid_){}

bool MuonVeto::passes(const Event & event){
  assert(event.jets); // if this fails, it probably means jets are not read in                                                                                                                                                                                             
  assert(event.muons); // if this fails, it probably means jets are not read in                                                                                                                                                                                               
 
  for(const auto & muon : *event.muons) {
    if(muid(muon,event)) {
      for(const auto & jet : *event.jets){
	if(deltaR(jet,muon)  < deltaR_min ) return false;
      }
    }
  }

  return true;

}

ElectronVeto::ElectronVeto(const ElectronId & eleid_, float deltaR_min_): deltaR_min(deltaR_min_), eleid(eleid_){}

bool ElectronVeto::passes(const Event & event){
  assert(event.jets); // if this fails, it probably means jets are not read in                                                                                                                                                                                             
  assert(event.electrons); // if this fails, it probably means jets are not read in                                                                                                                                                                                           

  for(const auto & electron : *event.electrons){
    if(eleid(electron,event)) {
	for(const auto & jet : *event.jets){
	  if(deltaR(jet,electron)  < deltaR_min) return false;
	}
      }
  }
  return true;
  
}



GenHbbEventSelection::GenHbbEventSelection(){}

bool GenHbbEventSelection::passes(const Event & event, Jet & jet){
  std::vector<GenParticle> genQuarks;
  int associatedQuarks=0;
  //bool Higgsboson = false;
  if(PRINT) cout << " GenHbbEventSelection" << endl;
  for(auto genp:*event.genparticles){
    if(PRINT) cout << "  genparticle pdg id = " << genp.pdgId() << endl; 
    if(abs(genp.pdgId())==25){
      if(PRINT) cout << "    I have a Higgs boson!" << endl;
      //Higgsboson = true;
      if(abs(event.genparticles->at(genp.daughter1()).pdgId())==5 && abs(event.genparticles->at(genp.daughter2()).pdgId())==5){
       genQuarks.push_back( event.genparticles->at(genp.daughter1()) ); 
       genQuarks.push_back( event.genparticles->at(genp.daughter2()) );
      }
    }
    
    if(genQuarks.size()>1)
    {
      associatedQuarks=0;
      for(unsigned int i=0; i<genQuarks.size(); ++i){
         if( deltaR(jet.v4(),genQuarks[i].v4()) < 0.8) associatedQuarks +=1;
      } 
      if(associatedQuarks != 2){
         genQuarks.clear();
         associatedQuarks=0;
         continue;
      }
      else break;
    }
  }
  
  //if(Higgsboson && genQuarks.size()>2) cout << "WARNING: found more than two daughters for Higgs boson!" << endl;
  //if(Higgsboson && genQuarks.size()<2) cout << "WARNING: found less than two daughters for Higgs boson!" << endl;
  
  if (associatedQuarks == 2) return true;
  else return false;
}


GenVqqEventSelection::GenVqqEventSelection(){}

std::tuple<bool, bool> GenVqqEventSelection::passes(const Event & event, Jet & jet){
  assert(event.genparticles); // if this fails, it probably means genparticles are not read in    
  std::vector<GenParticle> genQuarks;
  bool isZbb = false;
  int associatedQuarks=0;
  //bool vectorboson = false;
  if(PRINT) cout << " GenVqqEventSelection" << endl;
  for(auto genp:*event.genparticles){
    if(PRINT) cout << "  genparticle pdg id = " << genp.pdgId() << endl;   
    if(abs(genp.pdgId())==24 || abs(genp.pdgId())==23){
       if(PRINT) cout << "    I have a vector boson!" << endl;
       //vectorboson = true;
       int dau1 = abs(event.genparticles->at(genp.daughter1()).pdgId());
       if(PRINT) cout << "    I have one daughter" << endl;
       int dau2 = abs(event.genparticles->at(genp.daughter2()).pdgId());
       if(PRINT) cout << "    I have 2 daughters" << endl;
       if(PRINT) cout << "      daughters: " << event.genparticles->at(genp.daughter1()).pdgId() << " and " << event.genparticles->at(genp.daughter2()).pdgId() << endl;
       if( (dau1 >=1 && dau1 <=6) && (dau2 >=1 && dau2 <=6) ){
          genQuarks.push_back( event.genparticles->at(genp.daughter1()) ); 
          genQuarks.push_back( event.genparticles->at(genp.daughter2()) );
	  if( dau1==5 && dau2==5 && abs(genp.pdgId())==23 ) isZbb = true;
       }
    }

    if(genQuarks.size()>1)
    {
      associatedQuarks=0;
      for(unsigned int i=0; i<genQuarks.size(); ++i){
         if( deltaR(jet.v4(),genQuarks[i].v4()) < 0.8) associatedQuarks +=1;
      } 
      if(associatedQuarks != 2){
         genQuarks.clear();
         associatedQuarks=0;
	 isZbb=false;
         continue;
      }
      else break;
    }
  }

  if (associatedQuarks == 2) return std::make_tuple(true,isZbb);
  else return std::make_tuple(false,isZbb);
  
}   

GenTopHadrEventSelection::GenTopHadrEventSelection(){}

bool GenTopHadrEventSelection::passes(const Event & event, Jet & jet){
  assert(event.genparticles); // if this fails, it probably means genparticles are not read in    
  std::vector<GenParticle> genQuarks;
  int associatedQuarks=0;
  if(PRINT) cout << " GenTopHadrEventSelection" << endl;
  for(auto genp:*event.genparticles){
    //if(PRINT) cout << "  genparticle pdg id = " << genp.pdgId() << endl;   
    if(abs(genp.pdgId())==6){
       if(PRINT) cout << "    I have a top quark!" << endl;
       //vectorboson = true;
       int dau1 = abs(event.genparticles->at(genp.daughter1()).pdgId());
       if(PRINT) cout << "    I have one daughter" << endl;
       int dau2 = abs(event.genparticles->at(genp.daughter2()).pdgId());
       if(PRINT) cout << "    I have 2 daughters" << endl;
       if(PRINT) cout << "      daughters: " << event.genparticles->at(genp.daughter1()).pdgId() << " and " << event.genparticles->at(genp.daughter2()).pdgId()  << endl;
       
       if( dau1 == 24 && dau2 == 5 ){
        if(PRINT) cout << "CASE 1: dau1=" << dau1 << " - dau2=" << dau2;
        genQuarks.push_back( event.genparticles->at(genp.daughter2()) ); //this is the b
        int W_dau1 = abs(event.genparticles->at(event.genparticles->at(genp.daughter1()).daughter1()).pdgId());
	int W_dau2 = abs(event.genparticles->at(event.genparticles->at(genp.daughter1()).daughter2()).pdgId());
	if(PRINT) cout << " with W decay to " << W_dau1 << "," << W_dau2 << endl;
	if( (W_dau1 >=1 && W_dau1 <=6) && (W_dau2 >=1 && W_dau2 <=6) ){
	  if(PRINT) cout << "  --> FOUND HADRONIC W!" << endl;
          genQuarks.push_back( event.genparticles->at(event.genparticles->at(genp.daughter1()).daughter1()) ); 
          genQuarks.push_back( event.genparticles->at(event.genparticles->at(genp.daughter1()).daughter2()) );	
	}
       }//close case dau1=W, dau2=b
       else if( dau1 == 5 && dau2 == 24){
        if(PRINT) cout << "CASE 2: dau1=" << dau1 << " - dau2=" << dau2;
        int W_dau1 = abs(event.genparticles->at(event.genparticles->at(genp.daughter2()).daughter1()).pdgId());
	int W_dau2 = abs(event.genparticles->at(event.genparticles->at(genp.daughter2()).daughter2()).pdgId());    
	if(PRINT) cout << " with W decay to " << W_dau1 << "," << W_dau2 << endl;
	if( (W_dau1 >=1 && W_dau1 <=6) && (W_dau2 >=1 && W_dau2 <=6) ){
	  if(PRINT) cout << "  --> FOUND HADRONIC W!" << endl;
          genQuarks.push_back( event.genparticles->at(event.genparticles->at(genp.daughter2()).daughter1()) ); 
          genQuarks.push_back( event.genparticles->at(event.genparticles->at(genp.daughter2()).daughter2()) );	
	}           
       }//close case dau1=b, dau2=W
    }

    if(genQuarks.size()>1)
    {
      associatedQuarks=0;
      for(unsigned int i=0; i<genQuarks.size(); ++i){
         if( deltaR(jet.v4(),genQuarks[i].v4()) < 0.8) associatedQuarks +=1;
      } 
      if(associatedQuarks != 3){
         genQuarks.clear();
         associatedQuarks=0;
         continue;
      }
      else break;
    }
  }
  
  //if(vectorboson && genQuarks.size()>2) cout << "WARNING: found more than two daughters for Vector boson!" << endl;
  //if(vectorboson && genQuarks.size()<2) cout << "WARNING: found less than two daughters for Vector boson! Are you using an inclusive sample?" << endl;

  if (associatedQuarks == 3) return true;
  else return false;
} 

VBFjetSelection::VBFjetSelection(Context & ctx, string const & VBFjet_, float deta_min_, float mjj_min_ ):VBFjet(VBFjet_), deta_min(deta_min_), mjj_min(mjj_min_){
  h_VBFjet = ctx.get_handle<vector<Jet>>(VBFjet);
}

bool VBFjetSelection::passes(const Event & event){
  if (!event.is_valid(h_VBFjet)) {
    cerr << "VBF selections: Handle not valid!\n";
    assert(false);
  }

  vector<Jet> vbfjets = event.get(h_VBFjet);
  if(vbfjets.size() < 2) return false;
  const auto & jet0 = vbfjets.at(0);
  const auto & jet1 = vbfjets.at(1);
  auto deta = fabs(jet0.eta() - jet1.eta());
  if(deta < deta_min) return false;
  if(jet0.eta() * jet1.eta() > 0) return false;
  if(inv_mass_safe(jet0.v4()+jet1.v4()) < mjj_min ) return false;
  return true;
}


///Brute Force Decorrelation
BruteForceDecorrelation::BruteForceDecorrelation( uhh2::Context & ctx, string percentage_, string folder){

  percentage = percentage_;
  string h_name = "DeepBoosted_ZHbbvsQCD_"+percentage;
  string h_name2 = "DeepBoosted_WvsQCD_"+percentage;
  h_l1_DeepBoosted_ZHbbvsQCD= ctx.declare_event_output<float>("jj_l1_"+h_name+"_"+folder);
  h_l1_DeepBoosted_WvsQCD= ctx.declare_event_output<float>("jj_l1_"+h_name2+"_"+folder);
  h_l2_DeepBoosted_ZHbbvsQCD= ctx.declare_event_output<float>("jj_l2_"+h_name+"_"+folder);
  h_l2_DeepBoosted_WvsQCD= ctx.declare_event_output<float>("jj_l2_"+h_name2+"_"+folder);
  
  std::string filename_zh = "UHHNtupleConverter/data/DDTMap_ZHbbvsQCD/"+folder+"/myDeepBoostedMap"+percentage+".root" ;
  std::string filename_w = "UHHNtupleConverter/data/DDTMap_WvsQCD/"+folder+"/myDeepBoostedMap"+percentage+".root"; 
  infile_ZHbbvsQCD.reset(TFile::Open(locate_file(filename_zh).c_str()));
  infile_WvsQCD.reset(TFile::Open(locate_file(filename_w).c_str()));
 
  map_ZHbbvsQCD = (TH2F*)infile_ZHbbvsQCD->Get("DeepBoosted_ZHbbvsQCD_v_rho_v_pT_yx");
  map_WvsQCD = (TH2F*)infile_WvsQCD->Get("DeepBoosted_WvsQCD_v_rho_v_pT_yx");
  
}
    
bool BruteForceDecorrelation::process(Event & event, TopJet const* jet1_, TopJet const* jet2_){

  //find cut value for jet1
  pt_bin = map_ZHbbvsQCD->GetYaxis()->FindFixBin(jet1_->pt());
  if(pt_bin > map_ZHbbvsQCD->GetYaxis()->GetNbins()){
    pt_bin = map_ZHbbvsQCD->GetYaxis()->GetNbins();
  }else if(pt_bin <=0){
    pt_bin = 1;
  }
  
  x = 2 * TMath::Log(jet1_->softdropmass()/ jet1_->pt());
  x_bin = map_ZHbbvsQCD->GetXaxis()->FindFixBin(x);
  if(x_bin > map_ZHbbvsQCD->GetXaxis()->GetNbins()){
    x_bin = map_ZHbbvsQCD->GetXaxis()->GetNbins();
  }else if(x_bin <= 0){
    x_bin = 1;
  }
  
  cut_ZHbbvsQCD = map_ZHbbvsQCD->GetBinContent(x_bin,pt_bin);
  cut_WvsQCD = map_WvsQCD->GetBinContent(x_bin,pt_bin);        
  event.set(h_l1_DeepBoosted_ZHbbvsQCD,cut_ZHbbvsQCD);
  event.set(h_l1_DeepBoosted_WvsQCD,cut_WvsQCD);

  //find cut value for jet2
  pt_bin = map_ZHbbvsQCD->GetYaxis()->FindFixBin(jet2_->pt());
  if(pt_bin > map_ZHbbvsQCD->GetYaxis()->GetNbins()){
    pt_bin = map_ZHbbvsQCD->GetYaxis()->GetNbins();
  }else if(pt_bin <=0){
    pt_bin = 1;
  }
  
  x = 2 * TMath::Log(jet2_->softdropmass()/ jet2_->pt());
  x_bin = map_ZHbbvsQCD->GetXaxis()->FindFixBin(x);
  if(x_bin > map_ZHbbvsQCD->GetXaxis()->GetNbins()){
    x_bin = map_ZHbbvsQCD->GetXaxis()->GetNbins();
  }else if(x_bin <= 0){
    x_bin = 1;
  }
  
  cut_ZHbbvsQCD = map_ZHbbvsQCD->GetBinContent(x_bin,pt_bin);
  cut_WvsQCD = map_WvsQCD->GetBinContent(x_bin,pt_bin);        
  event.set(h_l2_DeepBoosted_ZHbbvsQCD,cut_ZHbbvsQCD);
  event.set(h_l2_DeepBoosted_WvsQCD,cut_WvsQCD);
    
  return true;
}






