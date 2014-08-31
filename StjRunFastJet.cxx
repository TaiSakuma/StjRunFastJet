// $Id$
// Copyright (C) 2010 Tai Sakuma <sakuma@bnl.gov>
#include "StjRunFastJet.h"

#include <StJetMaker/StjFourVecForJetFinder.h>
#include <StJetMaker/StjEtaToDetectorEta.h>

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/SISConePlugin.hh>
#include <fastjet/CDFMidPointPlugin.hh>
#include <fastjet/D0RunIIConePlugin.hh>

#include <TLorentzVector.h>
#include <iostream>

ClassImp(StjRunFastJet)

using namespace std;


//____________________________________________________________________________________________||
void StjRunFastJet::Init()
{

}
//____________________________________________________________________________________________||
StjRunFastJet::StjRunFastJet()
  : _jetDefinition(new fastjet::JetDefinition())
{

}
//____________________________________________________________________________________________||
StjRunFastJet::StjRunFastJet(JetAlgorithm jet_algorithm, double R, RecombinationScheme recomb_scheme, Strategy strategy)
  : _jetDefinition(new fastjet::JetDefinition((fastjet::JetAlgorithm)jet_algorithm, R, (fastjet::RecombinationScheme)recomb_scheme, (fastjet::Strategy)strategy))
{

}
//____________________________________________________________________________________________||
StjRunFastJet::StjRunFastJet(const void* plugin)
  : _jetDefinition(new fastjet::JetDefinition((const fastjet::JetDefinition::Plugin*) plugin))
{

}
//____________________________________________________________________________________________||
StjRunFastJet StjRunFastJet::SISCone(double cone_radius, double overlap_threshold)
{
  fastjet::JetDefinition::Plugin *plugin = new fastjet::SISConePlugin (cone_radius, overlap_threshold);
  return StjRunFastJet((void*)plugin);
}

//____________________________________________________________________________________________||
StjRunFastJet StjRunFastJet::CDFMidpoint(double cone_radius, double overlap_threshold, double seed_threshold, double cone_area_fraction)
{
  fastjet::JetDefinition::Plugin *plugin = new fastjet::CDFMidPointPlugin (cone_radius, overlap_threshold, seed_threshold, cone_area_fraction);
  return StjRunFastJet((void*)plugin);
}

//____________________________________________________________________________________________||
StjRunFastJet StjRunFastJet::D0RunIICone(double cone_radius, double min_jet_Et, double split_ratio)
{
  fastjet::JetDefinition::Plugin *plugin = new fastjet::D0RunIIConePlugin (cone_radius, min_jet_Et, split_ratio);
  return StjRunFastJet((void*)plugin);
}

//____________________________________________________________________________________________||
StjJetList StjRunFastJet::operator()(const StjFourVecList& fourVecList)
{
  vector<fastjet::PseudoJet> input_particles = create_input_particles(fourVecList);
  fastjet::ClusterSequence clust_seq(input_particles, *_jetDefinition);
  return create_jetList(fourVecList, clust_seq);
}
//____________________________________________________________________________________________||
std::vector<fastjet::PseudoJet> StjRunFastJet::create_input_particles(const StjFourVecList& fourVecList)
{
  vector<fastjet::PseudoJet> input_particles;
  for (unsigned int i = 0; i < fourVecList.size(); ++i)
    {
      TLorentzVector l;
      l.SetPtEtaPhiM(fourVecList[i].pt, fourVecList[i].eta, fourVecList[i].phi, fourVecList[i].m);
      fastjet::PseudoJet pseudojet(l.Px(), l.Py(), l.Pz(), l.E());
      pseudojet.set_user_index(i);
      input_particles.push_back(pseudojet); 
    }
  return input_particles;
}
//____________________________________________________________________________________________||
StjJetList StjRunFastJet::create_jetList(const StjFourVecList& fourVecList, fastjet::ClusterSequence& clust_seq)
{
  vector<fastjet::PseudoJet> jets = clust_seq.inclusive_jets(0.0);

  StjJetList jetList;
  int jetId(1);

  vector<fastjet::PseudoJet> sorted_jets = sorted_by_pt(jets);  
  for (unsigned int i = 0; i < sorted_jets.size(); i++)
    {
      StjJet jet;
      jet.jetId = jetId++;
      TLorentzVector l(sorted_jets[i].px(), sorted_jets[i].py(), sorted_jets[i].pz(), sorted_jets[i].E());
      jet.pt =  l.Pt();
      jet.eta = l.Eta();
      jet.phi = l.Phi();
      jet.m =   l.M();

      vector<fastjet::PseudoJet> constituents = clust_seq.constituents(sorted_jets[i]);
      for(vector<fastjet::PseudoJet>::const_iterator p4 = constituents.begin(); p4 != constituents.end(); ++p4)
	{
	  jet.runNumber = fourVecList[(*p4).user_index()].runNumber;
	  jet.eventId = fourVecList[(*p4).user_index()].eventId;
	  jet.vertexZ = fourVecList[(*p4).user_index()].vertexZ;
	  jet.fourVecList.push_back(fourVecList[(*p4).user_index()]);
	}
      jet.neuRt = computeNeuRt(jet.fourVecList);
      StjEtaToDetectorEta eta2deta;
      jet.detectorEta = eta2deta(jet.eta, jet.vertexZ);
      jetList.push_back(jet);
    }
  return jetList;
}
//____________________________________________________________________________________________||
double StjRunFastJet::computeNeuRt(const StjFourVecList& fourList)
{
  double totalEt = 0.0;
  double neutralEt = 0.0;
  for(StjFourVecList::const_iterator it = fourList.begin(); it != fourList.end(); ++it) {
    TLorentzVector p4;
    p4.SetPtEtaPhiM((*it).pt, (*it).eta, (*it).phi, (*it).m);
    totalEt += p4.Et();
    if((*it).type == 2) neutralEt += p4.Et();
  }
  return (totalEt) ? neutralEt/totalEt: 0.0;
}
//____________________________________________________________________________________________||

