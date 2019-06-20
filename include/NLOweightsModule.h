#pragma once


/* This file is meant to integrate in this converter a script done by Daniela Schaefer.
 */
#include <iostream>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ObjectIdUtils.h"

#include "TFile.h"
#include "TH1F.h"

using namespace std;

/// saving NLOweights to the output tree
class NLOweight: public uhh2::AnalysisModule {
 public:

  explicit NLOweight(uhh2::Context & ctx,
		     const std::string & NLOWeightsFilename="",
		     const std::string & s_jj_l1_mergedVTruth="jj_l1_mergedVTruth",
		     const std::string & s_jj_l1_gen_pt="jj_l1_gen_pt",
		     const std::string & s_jj_l2_gen_pt="jj_l2_gen_pt");
  virtual ~NLOweight() {};
  virtual bool process(uhh2::Event & event) override;
  float getWeight(float gen_pt);
  float getEWCorr(float gen_pt);

 private:
  uhh2::Event::Handle<bool> jj_l1_mergedVTruth;
  uhh2::Event::Handle<float> jj_l1_gen_pt,jj_l2_gen_pt;
  std::unique_ptr<TFile> NLOWeightsFile;
  std::unique_ptr<TH1F> h_kfactor;
  std::unique_ptr<TH1F> h_ewcorr; //nominal weight
  uhh2::Event::Handle<float>  m_o_kfactor;
};
