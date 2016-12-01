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
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
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
      void reset(const bool);
      virtual void endJob() override; 
// ----------member data ---------------------------
  edm::EDGetTokenT<reco::PFTauRefVector> tauTag_;
  edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> Mu1Mu2_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleTag_;
  edm::EDGetTokenT<edm::ValueMap<reco::MuonRefVector>>  jetMuonMapTag_;
  std::vector<double> muHadMassBins_;
  std::vector<double> FourBInvMassBins_;
  TFile *out_;
  TH1F *muHadMass_;
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
  Mu1Mu2_(consumes<edm::RefVector<std::vector<reco::Muon>>>(iConfig.getParameter<edm::InputTag>("Mu1Mu2"))),
  genParticleTag_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleTag"))),
  jetMuonMapTag_(consumes<edm::ValueMap<reco::MuonRefVector> >(iConfig.getParameter<edm::InputTag>("jetMuonMapTag"))),
  muHadMassBins_(iConfig.getParameter<std::vector<double> >("muHadMassBins")),
  FourBInvMassBins_(iConfig.getParameter<std::vector<double>>("FourBInvMassBins")),
  outFileName_(iConfig.getParameter<std::string>("outFileName"))
{
   //now do what ever initialization is needed
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
 
  edm::Handle<edm::RefVector<std::vector<reco::Muon>>> pMu1Mu2;
  iEvent.getByToken(Mu1Mu2_, pMu1Mu2);
  
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
  std::vector<reco::MuonRef> removedMuonRefs;

   for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); 
	 iMuon != removedMuons.end(); ++iMuon) { removedMuonRefs.push_back(*iMuon); }
   const double muHadMass = 
      (removedMuonRefs[removedMuonRefs.size() - 1]->p4() + (*iTau)->p4()).M();

   double muHadMassForPlotting = muHadMass;
   muHadMass_->Fill(muHadMassForPlotting, PUWeight*tauHadPTWeight*HiggsPTWeight);
   
   //std::vector<reco::Muon*> Mu1Mu2Ptrs;
 // reco::Muon* HighestPtMu1Mu2;
 // reco::Muon* LowestPtMu1Mu2;
 // if((*pMu1Mu2)[0]->pt()> (*pMu1Mu2)[1]->pt()){
  //  HighestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[0])));
  //  LowestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[1])));
   // }
 // else
  //  {
   // HighestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[1])));
   // LowestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[0])));
   // }

       //const double mumuInvMass =
         //(HighestPtMu1Mu2->p4()+LowestPtMu1Mu2->p4()).M();
       //const double FourBInvMass=(HighestPtMu1Mu2->p4()+LowestPtMu1Mu2->p4()+removedMuonRefs[removedMuonRefs.size() - 1]->p4() + (*iTau)->p4()).M();
    ++iTau;
  }
}

void 
MuMuTauTauRecoAnalyzer::reset(const bool doDelete)
{
  if(doDelete && (out_ !=NULL)) delete out_;
  out_=NULL;
  if(doDelete && (muHadMass_!=NULL)) delete muHadMass_;
  muHadMass_=NULL;
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuMuTauTauRecoAnalyzer::beginJob()
{
  out_= new TFile(outFileName_.c_str(),"RECREATE");
  muHadMass_=new TH1F("muHadMass", ";H750a09 m_{#mu+X} (GeV);", muHadMassBins_.size()-1, &muHadMassBins_[0]);  

}


// ------------ method called once each job just after ending the event loop  ------------
void 
MuMuTauTauRecoAnalyzer::endJob() 
{
  out_->cd();
  TCanvas muHadMassCanvas("muHadMassCanvas","",600,600);
  Common::draw1DHistograms(muHadMassCanvas, muHadMass_);
  out_->cd();
  muHadMassCanvas.Write();
  out_->Close();
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
