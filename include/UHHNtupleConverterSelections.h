#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/core/include/AnalysisModule.h"

#include "TFile.h"

#include "TH2F.h"



using namespace std;
using namespace uhh2;

namespace uhh2examples {
    
/* Select events with at least two jets in which the leading two jets have deltaphi > 2.7 and the third jet pt is
 * below 20% of the average of the leading two jets, where the minimum deltaphi and
 * maximum third jet pt fraction can be changed in the constructor.
 * The jets are assumed to be sorted in pt.
 */
class DijetSelection: public uhh2::Selection {
public:
    DijetSelection(float deta_max = 2.7f, float mjj_min = 0.0f);
    virtual bool passes(const uhh2::Event & event) override;
private:
    float deta_max, mjj_min;
};


class MuonVeto: public uhh2::Selection {
public:
    MuonVeto(const MuonId & muid, float deltR_min = 0.8f  );
    virtual bool passes(const uhh2::Event & event) override;
private:
  float deltaR_min;
  const MuonId & muid;
};

class ElectronVeto: public uhh2::Selection {
public:
    ElectronVeto(const ElectronId & eleid, float deltR_min = 0.8f  );
    virtual bool passes(const uhh2::Event & event) override;
private:
  float deltaR_min;
  const ElectronId & eleid;
};




class GenHbbEventSelection: public uhh2::Selection {
public:
    GenHbbEventSelection();
    virtual bool passes(const uhh2::Event & event) override {return false;};
    bool passes(const uhh2::Event & event, Jet & jet_);
    
private:    

};

class GenVqqEventSelection: public uhh2::Selection {
public:
    GenVqqEventSelection();
    virtual bool passes(const uhh2::Event & event) override {return false;};
    std::tuple<bool, bool> passes(const uhh2::Event & event, Jet & jet_);
    
private:    

};

class GenTopHadrEventSelection: public uhh2::Selection {
public:
    GenTopHadrEventSelection();
    virtual bool passes(const uhh2::Event & event) override {return false;};
    bool passes(const uhh2::Event & event, Jet & jet_);
    
private:    

};


class VBFjetSelection: public uhh2::Selection {
public:
    VBFjetSelection(uhh2::Context & ctx, const string & VBFjet = "jetsAk4Puppi", float deta_min = 4.5f, float mjj_min = 800.0f);
    virtual bool passes(const uhh2::Event & event) override;
private:
    const string & VBFjet;
    float deta_min, mjj_min;
    Event::Handle<vector<Jet>> h_VBFjet;
  };



  class BruteForceDecorrelation: public uhh2::AnalysisModule{

  public:
    explicit BruteForceDecorrelation(uhh2::Context & ctx,string percentage_, string folder);
    //virtual bool process(uhh2::Event & event) override; 
    virtual bool process(uhh2::Event & event) override {return false;};
    bool process(uhh2::Event & event, TopJet const* jet1_, TopJet const* jet2_);

  private:
    string percentage;
    uhh2::Event::Handle<float> h_l1_DeepBoosted_ZHbbvsQCD;
    uhh2::Event::Handle<float> h_l1_DeepBoosted_WvsQCD;
    uhh2::Event::Handle<float> h_l2_DeepBoosted_ZHbbvsQCD;
    uhh2::Event::Handle<float> h_l2_DeepBoosted_WvsQCD;
    std::unique_ptr<TFile> infile_ZHbbvsQCD, infile_WvsQCD;
    TH2F* map_ZHbbvsQCD, *map_WvsQCD;
    int pt_bin,x_bin;
    double x;
    float cut_ZHbbvsQCD,cut_WvsQCD;
  };



}
