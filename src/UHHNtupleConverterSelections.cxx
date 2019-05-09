#include "UHH2/UHHNtupleConverter/include/UHHNtupleConverterSelections.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Particle.h"
#include "UHH2/common/include/Utils.h"

#include <stdexcept>
#include <iostream>

using namespace uhh2examples;
using namespace uhh2;
using namespace std;

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

MuonVeto::MuonVeto(float deltaR_min_, const boost::optional<MuonId> & muid_): deltaR_min(deltaR_min_), muid(muid_){}

bool MuonVeto::passes(const Event & event){
  assert(event.topjets); // if this fails, it probably means jets are not read in                                                                                                                                                                                             
  assert(event.muons); // if this fails, it probably means jets are not read in                                                                                                                                                                                               
  if(muid){
    for(const auto & muons : *event.muons){
      for(const auto & topjets : *event.topjets){
	if(deltaR(topjets,muons)  < deltaR_min) return false;
	else return true;
      }
    }
  }
  return true;

}

ElectronVeto::ElectronVeto(float deltaR_min_, const boost::optional<ElectronId> & eleid_): deltaR_min(deltaR_min_), eleid(eleid_){}

bool ElectronVeto::passes(const Event & event){
  assert(event.topjets); // if this fails, it probably means jets are not read in                                                                                                                                                                                             
  assert(event.electrons); // if this fails, it probably means jets are not read in                                                                                                                                                                                           
  if(eleid){
    for(const auto & topjets : *event.topjets){
      for(const auto & electrons : *event.electrons){
	if(deltaR(topjets,electrons)  < deltaR_min) return false;
	else return true;
      }
    }
  }
  return true;
  
}



GenHbbEventSelection::GenHbbEventSelection(){}
    
bool GenHbbEventSelection::passes(const Event & event, Jet & jet){
  std::vector<GenParticle> genQuarks;
  bool Higgsboson = false;
  for(auto genp:*event.genparticles){
    if(abs(genp.pdgId())==25){
      Higgsboson = true;
      if(abs(event.genparticles->at(genp.daughter1()).pdgId())==5 && abs(event.genparticles->at(genp.daughter2()).pdgId())==5){
       genQuarks.push_back( event.genparticles->at(genp.daughter1()) ); 
       genQuarks.push_back( event.genparticles->at(genp.daughter2()) );
      }
    }
    if(genQuarks.size()>1) break;
  }
  
  if(Higgsboson && genQuarks.size()>2) cout << "WARNING: found more than two daughters for Higgs boson!" << endl;
  if(Higgsboson && genQuarks.size()<2) cout << "WARNING: found less than two daughters for Higgs boson!" << endl;
  
  int associatedQuarks=0;
  for(unsigned int i=0; i<genQuarks.size(); ++i){
   if( deltaR(jet.v4(),genQuarks[i].v4()) < 0.8) associatedQuarks +=1;
  } 
  
  if (associatedQuarks ==2) return true;
  else return false;
}

GenVqqEventSelection::GenVqqEventSelection(){}
    
bool GenVqqEventSelection::passes(const Event & event, Jet & jet){
  std::vector<GenParticle> genQuarks;
  bool vectorboson = false;
  for(auto genp:*event.genparticles){
    if(abs(genp.pdgId())==24 || abs(genp.pdgId())==23){
       vectorboson = true;
       int dau1 = abs(event.genparticles->at(genp.daughter1()).pdgId());
       int dau2 = abs(event.genparticles->at(genp.daughter2()).pdgId());
       if( (dau1 >=1 && dau1 <=6) && (dau2 >=1 && dau2 <=6) ){
        genQuarks.push_back( event.genparticles->at(genp.daughter1()) ); 
        genQuarks.push_back( event.genparticles->at(genp.daughter2()) );
       }
    }
    if(genQuarks.size()>1) break;
  }
  
  if(vectorboson && genQuarks.size()>2) cout << "WARNING: found more than two daughters for Vector boson!" << endl;
  if(vectorboson && genQuarks.size()<2) cout << "WARNING: found less than two daughters for Vector boson! Are you using an inclusive sample?" << endl;

 
  int associatedQuarks=0;
  for(unsigned int i=0; i<genQuarks.size(); ++i){
   if( deltaR(jet.v4(),genQuarks[i].v4()) < 0.8) associatedQuarks +=1;
  } 
  
  if (associatedQuarks ==2) return true;
  else return false;
}