// -*- C++ -*-
//
// Package:    GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer
// Class:      MuMuTauTauRecoAnalyzer
// 
/**\class MuMuTauTauRecoAnalyzer MuMuTauTauRecoAnalyzer.cc GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer/plugins/MuMuTauTauRecoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mengyao Shi
//         Created:  Fri, 11 Mar 2016 09:17:50 GMT
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "Tools/Common/interface/Common.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TSystem.h"
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Pruner.hh"
//trigger stuff
//#include "DataFormats/HLTReco/interface/TriggerObject.h"
//#include "DataFormats/Common/interface/TriggerResults.h"
//#include "DataFormats/HLTReco/interface/TriggerEvent.h"
//#include "HLTrigger/HLTanalyzers/interface/HLTInfo.h"
//#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
//#include "HLTrigger/HLTcore/interface/HLTConfigData.h"
//#include "FWCore/Common/interface/TriggerNames.h"
using namespace std;
using namespace edm;
using namespace reco;
using namespace fastjet;
using namespace trigger;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MuMuTauTauRecoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuMuTauTauRecoAnalyzer(const edm::ParameterSet&);
      ~MuMuTauTauRecoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override; 
 void reset(const bool);      
// ----------member data ---------------------------
      TH1F* muHadMass_;
      TH2F* InvMass2D_;
      TH1F* mumuInvMass_;
      TH1F* FourBInvMass_;
      TH1F* DR_;
      TH1F* DR23_;
      TH1F* DR13_;
  edm::EDGetTokenT<reco::PFTauRefVector> tauTag_;
  edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> muonTag1_;
  edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> muonTag2_; 
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleTag_;
  edm::EDGetTokenT<edm::ValueMap<reco::MuonRefVector>>  jetMuonMapTag_;
  std::vector<double> muHadMassBins_;
  std::vector<double> FourBInvMassBins_;
  TFile* out_;
 std::string outFileName_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuMuTauTauRecoAnalyzer::MuMuTauTauRecoAnalyzer(const edm::ParameterSet& iConfig):
  tauTag_(consumes<reco::PFTauRefVector>(iConfig.getParameter<edm::InputTag>("tauTag"))),
  muonTag1_(consumes<edm::RefVector<std::vector<reco::Muon>>>(iConfig.getParameter<edm::InputTag>("muonTag1"))),
  muonTag2_(consumes<edm::RefVector<std::vector<reco::Muon>>>(iConfig.getParameter<edm::InputTag>("muonTag2"))),
  genParticleTag_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleTag"))),
  jetMuonMapTag_(consumes<edm::ValueMap<reco::MuonRefVector> >(iConfig.getParameter<edm::InputTag>("jetMuonMapTag"))),
  muHadMassBins_(iConfig.getParameter<std::vector<double> >("muHadMassBins")),
  FourBInvMassBins_(iConfig.getParameter<std::vector<double>>("FourBInvMassBins")),
  outFileName_(iConfig.getParameter<std::string>("outFileName"))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   reset(false);
}


MuMuTauTauRecoAnalyzer::~MuMuTauTauRecoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuMuTauTauRecoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  edm::Handle<reco::PFTauRefVector> pTaus;
  iEvent.getByToken(tauTag_, pTaus);  
 
  edm::Handle<edm::RefVector<std::vector<reco::Muon>>> pMuons1;
  iEvent.getByToken(muonTag1_, pMuons1);
  
  edm::Handle<edm::RefVector<std::vector<reco::Muon>>> pMuons2;
  iEvent.getByToken(muonTag2_,pMuons2); 
 
  edm::Handle<reco::GenParticleCollection> pGenParticles;
   iEvent.getByToken(genParticleTag_, pGenParticles);

  edm::Handle<edm::ValueMap<reco::MuonRefVector> > pMuonJetMap;
  iEvent.getByToken(jetMuonMapTag_, pMuonJetMap);

  //fill an STL container
  std::vector<reco::GenParticle*> genObjPtrs; 
   for(typename std::vector<reco::GenParticle>::const_iterator iGenObj=pGenParticles->begin(); iGenObj!=pGenParticles->end(); ++iGenObj)
   {
     const unsigned int absPDGID=fabs(iGenObj->pdgId());
     if(absPDGID==13)
       genObjPtrs.push_back(const_cast<reco::GenParticle*>(&(*iGenObj)));
   }
  double PUWeight = 1.0;
  double HiggsPTWeight = 1.0;
   double tauHadPTWeight = 1.0;
  std::vector<reco::PFTauRef> pTSortedTaus;
  for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
       ++iTau) { pTSortedTaus.push_back(*iTau); }
  std::vector<reco::PFTauRef> taus = pTSortedTaus;
  std::vector<reco::PFTauRef>::const_iterator iTau = taus.begin();
  std::vector<reco::PFTauRef>::const_iterator endTau = taus.end();
  while (iTau != endTau) {
    const reco::PFJetRef& tauJetRef = (*iTau)->jetRef();
    const reco::MuonRefVector& removedMuons = (*pMuonJetMap)[tauJetRef];
    
    //make a collection of corrected old jets in |eta| < 2.4 not overlapping the W muon or tau
  int nearestGenObjKey=-1;
  std::vector<reco::MuonRef> removedMuonRefs;
  const reco::GenParticle* nearestGenObj= Common::nearestObject(*iTau, genObjPtrs, nearestGenObjKey);
   int pdgID=nearestGenObj->motherRef()->pdgId();
   if(fabs(pdgID)==511||fabs(pdgID)==521||fabs(pdgID)==411||fabs(pdgID)==441)
     DR_->Fill(reco::deltaR((*iTau)->eta(), (*iTau)->phi(), nearestGenObj->eta(), nearestGenObj->phi()));
   if(fabs(pdgID)==23)
     DR23_->Fill(reco::deltaR((*iTau)->eta(), (*iTau)->phi(), nearestGenObj->eta(), nearestGenObj->phi()));
   if(fabs(pdgID)==13)
     DR13_->Fill(reco::deltaR((*iTau)->eta(), (*iTau)->phi(), nearestGenObj->eta(), nearestGenObj->phi()));

   for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); 
	 iMuon != removedMuons.end(); ++iMuon) { removedMuonRefs.push_back(*iMuon); }
   const double muHadMass = 
      (removedMuonRefs[removedMuonRefs.size() - 1]->p4() + (*iTau)->p4()).M();

   double muHadMassForPlotting = muHadMass;
   muHadMass_->Fill(muHadMassForPlotting, PUWeight*tauHadPTWeight*HiggsPTWeight);
   for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iMuon2=pMuons2->begin();iMuon2!=pMuons2->end(); ++iMuon2)
   {
     for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iMuon1=pMuons1->begin(); iMuon1!=pMuons1->end();++iMuon1)
     {
       const double mumuInvMass =
         ((*iMuon1)->p4()+(*iMuon2)->p4()).M();
       const double FourBInvMass=((*iMuon1)->p4()+(*iMuon2)->p4()+removedMuonRefs[removedMuonRefs.size() - 1]->p4() + (*iTau)->p4()).M();
       InvMass2D_->Fill(mumuInvMass, muHadMassForPlotting);
       InvMass2D_->SetXTitle("Invariant mass of di-muon(GeV)");
       InvMass2D_->SetYTitle("Invariant mass m_{#mu+X} (GeV)");
       mumuInvMass_->Fill(mumuInvMass);
       mumuInvMass_->SetXTitle("Invariant mass of di-muon(GeV)");
       FourBInvMass_->Fill(FourBInvMass);
     }
   }
    ++iTau;
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuMuTauTauRecoAnalyzer::beginJob()
{
  out_ = new TFile(outFileName_.c_str(), "RECREATE");
  DR_=
    new TH1F("DR_nearestMu_tau", "DR between nearestMuon with mom absolute pdgID 511 521 411 441 to input tau", 100, 0, 1.0);
 DR23_=
    new TH1F("DR_nearestMu_tau", "DR between nearestMuon with mom pdgID 23 to input tau", 100, 0, 1.0);
 DR13_=
    new TH1F("DR_nearestMu_tau", "DR between nearestMuon with mom absolute pdgID 13 to input tau", 100, 0, 1.0);

  muHadMass_=
    new TH1F("H750a09 muHadMass", ";H750a09 m_{#mu+X} (GeV);", muHadMassBins_.size()-1, &muHadMassBins_[0]);
    muHadMass_->Sumw2();
  InvMass2D_=
    new TH2F("H750a09 InvMass",";H750a09 mumuInvtMass  (GeV);",14, 0, 14,14, 0, 14 );
  mumuInvMass_= 
    new TH1F("H750a09 mumuInvMass",";H750a09 mumuInvtMass", muHadMassBins_.size() - 1, &muHadMassBins_[0] ); 
  mumuInvMass_->Sumw2();
  FourBInvMass_=
    new TH1F("H750a09 four object final state Invariant Mass", "; H750a09 Four object invariant Mass", FourBInvMassBins_.size() -1, &FourBInvMassBins_[0]);
  FourBInvMass_->Sumw2();
}


// ------------ method called once each job just after ending the event loop  ------------
void 
MuMuTauTauRecoAnalyzer::endJob() 
{
 out_->cd();
  TCanvas muHadMassCanvas("muHadMassCanvas", "", 600, 600);
  TCanvas InvMass2DCanvas("InvMass2DCanvas","",600,600);
  TCanvas FourBInvMassCanvas("FourBInvMassCanvas","",600,600);
  TCanvas mumuInvMassCanvas("mumuInvMassCanvas","",600,600);
  TCanvas DRCanvas("DR_nearestMu(mom 511 521 441 411)_tau","",600,600);
  TCanvas DR23Canvas("DR_nearestMu(mom 23)_tau","",600,600);
  TCanvas DR13Canvas("DR_nearestMu(mom 13)_tau","",600,600);
  Common::draw1DHistograms(muHadMassCanvas, muHadMass_);
  Common::draw2DHistograms(InvMass2DCanvas, InvMass2D_);
  Common::draw1DHistograms(FourBInvMassCanvas, FourBInvMass_);
  Common::draw1DHistograms(mumuInvMassCanvas, mumuInvMass_);
  Common::draw1DHistograms(DRCanvas, DR_);
  Common::draw1DHistograms(DR23Canvas, DR23_);
  Common::draw1DHistograms(DR13Canvas, DR13_);
  out_->cd();
  muHadMassCanvas.Write();
  InvMass2DCanvas.Write();
  FourBInvMassCanvas.Write();
  mumuInvMassCanvas.Write();
  DRCanvas.Write();
  DR23Canvas.Write();
  DR13Canvas.Write();
  out_->Write();
  out_->Close();
}
void
MuMuTauTauRecoAnalyzer::reset(const bool doDelete)
{
  if (doDelete && (out_ != NULL)) delete out_;
  out_ = NULL;
 if (doDelete && (muHadMass_ != NULL)) delete muHadMass_;
  muHadMass_ = NULL;
  if (doDelete && (InvMass2D_ != NULL)) delete InvMass2D_;
  InvMass2D_ = NULL;
  if(doDelete && (mumuInvMass_!=NULL)) delete mumuInvMass_;
  mumuInvMass_=NULL;
  if(doDelete && (FourBInvMass_!=NULL)) delete FourBInvMass_;
  FourBInvMass_=NULL;
  if (doDelete &&(DR_!=NULL)) delete DR_;
  DR_ =NULL;
  if (doDelete &&(DR23_!=NULL)) delete DR23_;
  DR23_ =NULL;
  if (doDelete &&(DR13_!=NULL)) delete DR13_;
  DR13_ =NULL;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuMuTauTauRecoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuMuTauTauRecoAnalyzer);
