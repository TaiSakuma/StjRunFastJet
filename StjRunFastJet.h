// -*- mode: c++;-*-
// $Id$
// Copyright (C) 2010 Tai Sakuma <sakuma@bnl.gov>
#ifndef RUNFASTJET_H
#define RUNFASTJET_H

//____________________________________________________________________________________________||
#include <TObject.h>

#include "StjFourVecList.h"
#include "StjJetList.h"

//____________________________________________________________________________________________||
namespace fastjet {
  class JetDefinition;
  class PseudoJet;
  class ClusterSequence;
  class Plugin;
}

//____________________________________________________________________________________________||
class StjRunFastJet : public TObject {

public:

  //__________________________________________________________________________________________||
  enum JetAlgorithm {
    kt_algorithm = 0,
    cambridge_algorithm = 1,
    antikt_algorithm = 2, 
    genkt_algorithm = 3, 
    cambridge_for_passive_algorithm = 11,
    genkt_for_passive_algorithm = 13, 
    ee_kt_algorithm = 50,
    ee_genkt_algorithm = 53,
    plugin_algorithm = 99
  };

  enum RecombinationScheme {
    E_scheme = 0,
    pt_scheme = 1,
    pt2_scheme = 2,
    Et_scheme = 3,
    Et2_scheme = 4,
    BIpt_scheme = 5,
    BIpt2_scheme = 6,
    external_scheme = 99
  };

  enum Strategy {
    N2MinHeapTiled   = -4, 
    N2Tiled     = -3, 
    N2PoorTiled = -2, 
    N2Plain     = -1, 
    N3Dumb      =  0, 
    Best        =  1, 
    NlnN        =  2, 
    NlnN3pi     =  3, 
    NlnN4pi     =  4,
    NlnNCam4pi   = 14,
    NlnNCam2pi2R = 13,
    NlnNCam      = 12, // 2piMultD
    plugin_strategy = 999
  };

  //__________________________________________________________________________________________||
  StjRunFastJet();
  StjRunFastJet(JetAlgorithm jet_algorithm, double R, RecombinationScheme recomb_scheme, Strategy strategy);
  StjRunFastJet(const void* plugin);

  static StjRunFastJet SISCone(double cone_radius, double overlap_threshold);
  static StjRunFastJet CDFMidpoint(double cone_radius, double overlap_threshold, double seed_threshold = 1.0, double cone_area_fraction = 1.0);
  static StjRunFastJet D0RunIICone(double cone_radius, double min_jet_Et, double split_ratio = 0.5);


  virtual ~StjRunFastJet() { }

  void Init();

  StjJetList operator()(const StjFourVecList& fourList);

  //__________________________________________________________________________________________||
private:

  double computeNeuRt(const StjFourVecList& fourList);

  std::vector<fastjet::PseudoJet> create_input_particles(const StjFourVecList& fourList);

  StjJetList create_jetList(const StjFourVecList& fourVecList, fastjet::ClusterSequence& clust_seq);

  fastjet::JetDefinition* _jetDefinition;


  ClassDef(StjRunFastJet, 1)

};
//____________________________________________________________________________________________||


#endif // RUNFASTJET_H

