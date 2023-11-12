// -*- C++ -*-
//
// Package:    KShortLifetime/KShortLifetimeAnalyzer
// Class:      KShortLifetimeAnalyzer
//
/**\class KShortLifetimeAnalyzer KShortLifetimeAnalyzer.cc KShortLifetime/KShortLifetimeAnalyzer/plugins/KShortLifetimeAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Walker Sundquist <wwsundquist@ucsb.edu>
//         Created:  Tue, 29 Aug 2023 04:01:04 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <TGraph.h>
#include "TH2.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class KShortLifetimeAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit KShortLifetimeAnalyzer(const edm::ParameterSet&);
      ~KShortLifetimeAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfcands_;
      edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> kscands_;
      edm::EDGetTokenT<reco::VertexCollection> primaryVertices_;
      edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> lamcands_;
//      edm::EDGetTokenT<reco::VertexCollection> secondaryVertices_;

// Necessary for MC Truth access-----------------------------------------v
//      edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticles_;
//      edm::EDGetTokenT<reco::PackedGenParticleCollection> packedGenParticles_;
// End-------------------------------------------------------------------^

      edm::Service<TFileService> fs_;
      TTree* tree_;
      std::vector<float> kshort_mass;
      TH1D* h_kshort_mass;

      std::vector<float> pfNeutralHadron_mass;
      TH1D* h_pfNeutralHadron_mass;

      std::vector<float> dipion_mass;
      TH1D* h_dipion_mass;

      std::vector<float> kslifetimes;
      TH1D* h_kslifetimes;

      std::vector<float> ksdistances;
      TH1D* h_ksdistances;

      std::vector<float> lambda_mass;
      TH1D* h_lambda_mass;

      std::vector<float> lamdistances;
      TH1D* h_lamdistances;

      std::vector<float> lamlifetimes;
      TH1D* h_lamlifetimes;

      std::vector<float> ks_bkgd_lifetimes;
      TH1D* h_ks_bkgd_lifetimes;

      std::vector<float> ks_distVp_x;
      std::vector<float> ks_distVp_y;
      TH2F* h_ks_distVp;
//      TGraph* s_ks_distVp;
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
KShortLifetimeAnalyzer::KShortLifetimeAnalyzer(const edm::ParameterSet& iConfig)
 :
  pfcands_(consumes<pat::PackedCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("packedPFCandidates"))),
  kscands_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("slimmedKshortVertices"))),
  primaryVertices_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("offlineSlimmedPrimaryVertices"))),
  lamcands_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("slimmedLambdaVertices")))
  //secondaryVertices_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("offlineSlimmedSecondaryVertices")))
{
   //now do what ever initialization is needed
   tree_ = fs_->make<TTree>("tree", "tree");
   tree_->Branch("kshort_mass", &kshort_mass);
   h_kshort_mass = fs_->make<TH1D>("h_kshort_mass", "h_kshort_mass", 100, 0.4, 0.6);

   tree_->Branch("ks_distVp_x",&ks_distVp_x);
   tree_->Branch("ks_distVp_y",&ks_distVp_y);

   h_ks_distVp = fs_->make<TH2F>("h_ks_distVp","h_ks_distVp",30,0,10,75,0,25);

   tree_->Branch("pfNeutralHadron_mass", &pfNeutralHadron_mass);
   h_pfNeutralHadron_mass = fs_->make<TH1D>("h_pfNeutralHadron_mass", "h_pfNeutralHadron_mass", 100, 0.4, 0.6);

   tree_->Branch("dipion_mass", &dipion_mass);
   h_dipion_mass = fs_->make<TH1D>("h_dipion_mass", "h_dipion_mass", 200, 0, 2);

   tree_->Branch("kslifetimes", &kslifetimes);
   h_kslifetimes = fs_->make<TH1D>("h_kslifetimes", "h_kslifetimes", 30, 0, 8*pow(10,-10));

   tree_->Branch("ksdistances",&ksdistances);
   h_ksdistances = fs_->make<TH1D>("h_ksdistances", "h_ksdistances", 120, 0, 60);

   tree_->Branch("lambda_mass", &lambda_mass);
   h_lambda_mass = fs_->make<TH1D>("h_lambda_mass", "h_lambda_mass", 100, 1.0, 1.3);

   tree_->Branch("lamdistances", &lamdistances);
   h_lamdistances = fs_->make<TH1D>("h_lamdistances", "h_lamdistances", 20, 0, 40);

   tree_->Branch("lamlifetimes", &lamlifetimes);
   h_lamlifetimes = fs_->make<TH1D>("h_lamlifetimes", "h_lamlifetimes", 20, 0, 5.5*pow(10,-10));

   tree_->Branch("ks_bkgd_lifetimes", &ks_bkgd_lifetimes);
   h_ks_bkgd_lifetimes = fs_->make<TH1D>("h_ks_bkgd_lifetimes" , "h_ks_bkgd_lifetimes", 30, 0, 8*pow(10,-10));
}


KShortLifetimeAnalyzer::~KShortLifetimeAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

int count_infs = 0; //for some reason there are a lot of candidates that return *cycrad = inf or nan, so I'm counting them
int count_nans = 0;

int kscounter = 0;
int mctruth_ksNum = 0;
int mctruth_elseNum = 0;

//
// member functions
//

// ------------ method called for each event  ------------
void KShortLifetimeAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   kshort_mass.clear();
   dipion_mass.clear();
   kslifetimes.clear();
   ksdistances.clear();
   lambda_mass.clear();
   lamdistances.clear();
   lamlifetimes.clear();
   ks_bkgd_lifetimes.clear();

   std::cout << "Run: "<< iEvent.id().run() << " luminosity block: "<< iEvent.eventAuxiliary().luminosityBlock() <<" event: "<< iEvent.id().event() << std::endl;

   std::vector<float> pvXs;
   std::vector<float> pvYs;
   std::vector<float> pvZs;
   for(const auto& primaryVertex : iEvent.get(primaryVertices_) ) {
      pvXs.push_back(primaryVertex.position().x());
      pvYs.push_back(primaryVertex.position().y());
      pvZs.push_back(primaryVertex.position().z());
   }

   unsigned counter = 0;
   for(const auto& kscand : iEvent.get(kscands_) ) {
      std::cout << "KShort Candidate pt=" << kscand.pt() << " | mass=" << kscand.mass() << " | No. Daughters: " << kscand.numberOfDaughters() << " | vertex position (x,y,z): " << kscand.position().x() << ", " << kscand.position().y() << ", " << kscand.position().z() << std::endl;
      counter++;
      for (std::size_t i = 0; i < kscand.numberOfDaughters(); i = i + 1 ) {
         std::cout << "    Mass of daughter #" << i << " is " << kscand.daughter(i)->mass() << std::endl;
         std::cout << "    Charge of daughter #" << i << " is " << kscand.daughter(i)->charge()  << std::endl;
      }

      //mass cuts:
      float lowerlim = 0.468098; //5 sigma
      float upperlim = 0.528689;

      kshort_mass.push_back(kscand.mass());
      h_kshort_mass->Fill(kscand.mass());

      if (kscand.mass() < lowerlim or kscand.mass() > upperlim) {
         //Fill the background lifetimes with lifetime of candidates with masses outside cuts (cf. display.cxx)
         double delta_x = pow(pow(kscand.position().x()-pvXs[0],2) + pow(kscand.position().y()-pvYs[0],2) + pow(kscand.position().z()-pvZs[0],2),0.5);
         double tau = delta_x*pow(10,-2)*kscand.mass()*pow(kscand.energy(), -1)*pow(2.998e8, -1);
         ks_bkgd_lifetimes.push_back(tau);
         h_ks_bkgd_lifetimes->Fill(tau);
      } else {

         //kshort_mass.push_back(kscand.mass());
         //h_kshort_mass->Fill(kscand.mass());

         std::cout << "ID: " << kscand.pdgId() << std::endl;

         if (kscand.pdgId()==0) {
            mctruth_ksNum++;
         } else {
            mctruth_elseNum++;
         }

        //distance of No. 1 primary vertex from origin

        double pvdist = pow(pow(pvXs[0],2) + pow(pvYs[0],2) + pow(pvZs[0],2),0.5);
        std::cout << "PV DIST : " << pvdist << std::endl;

         double three_mag = pow((kscand.px()*kscand.px()) + (kscand.py()*kscand.py()) + (kscand.pz()*kscand.pz()),0.5);
         double delta_x = pow(pow(kscand.position().x()-pvXs[0],2) + pow(kscand.position().y()-pvYs[0],2) + pow(kscand.position().z()-pvZs[0],2),0.5); //had pvXs[1]?
         ksdistances.push_back(delta_x);
         h_ksdistances->Fill(delta_x);

         ks_distVp_x.push_back(three_mag);
         ks_distVp_y.push_back(delta_x);

         h_ks_distVp->Fill(three_mag, delta_x);

         double tau = delta_x*pow(10,-2)*kscand.mass()*pow(kscand.energy(), -1)*pow(2.998e8, -1);
         std::cout << "Lifetime : " << tau << std::endl;
         kslifetimes.push_back(tau);
         h_kslifetimes->Fill(tau);

      }

   }


//K-short MC truth---------------------------------------------------------------------------------------------v
// **** When uncommenting this code, but sure to also uncomment KShortLifetimeAnalyzer > Private > member data code ****

// Get generator information
//   edm::Handle<reco::GenParticleCollection> prunedGenParts;
//   iEvent.getByToken(prunedGenParticles_, prunedGenParts);
//   edm::Handle<pat::PackedGenParticleCollection> packedGenParts;
//   iEvent.getByToken(packedGenParticles_, packedGenParts);

   // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookGenParticleCandidate
   // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015
//   for(size_t i=0; i<prunedGenParts->size();++i){
//     const reco::Candidate * genParticle = &(*prunedGenParts)[i];
     // https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
//     if (abs(genParticle->pdgId()) != 310) continue;
//     std::cout << "index: "<<i<<" PdgID: " << genParticle->pdgId() << " pt " << genParticle->pt() << " eta: " << genParticle->eta() << " phi: " << genParticle->phi() << std::endl;
//     std::cout << "  status: "<<genParticle->status()<<std::endl;
//     std::cout << "  vertex position (x,y,z): "<<genParticle->position().x()<<", "<<genParticle->position().y()<<", "<<genParticle->position().z()<<std::endl;
//     std::cout << "  daughter:"<<std::endl;
//     for(size_t j=0; j<packedGenParts->size(); ++j) {
       //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection 
//       const reco::Candidate * motherInPrunedCollection = (*packedGenParts)[j].mother(0) ;
//       if(motherInPrunedCollection != nullptr && isAncestor( genParticle , motherInPrunedCollection)){
//         std::cout << "    PdgID: " << (*packedGenParts)[j].pdgId() << " pt " << (*packedGenParts)[j].pt() << " eta: " << (*packedGenParts)[j].eta() << " phi: " << (*packedGenParts)[j].phi() << std::endl;
//       }
//     }
//   }

//   for(size_t i=0; i<packedGenParts->size(); ++i) {
//    const pat::PackedGenParticle * genParticle = &(*packedGenParts)[i];
//     if (abs(genParticle->pdgId()) != 310) continue;
//     std::cout << "index: "<<i<<" PdgID: " << genParticle->pdgId() << " pt " << genParticle->pt() << " eta: " << genParticle->eta() << " phi: " << genParticle->phi() << std::endl;
//     std::cout<<  "  nDaugh: "<<genParticle->numberOfDaughters()<<std::endl;

//End----------------------------------------------------------------------------------------------------------^


   std::cout << "No. K-short Candidates: " << counter << std::endl;

   for(const auto& lamcand : iEvent.get(lamcands_) ) {
     
      if (lamcand.mass() < 1.11) {
         std::cout << "My name is Inigo Montoya!" << std::endl;
         break;
      } else if (lamcand.mass() > 1.13) {
         std::cout << "My name is Inigo Montoya!" << std::endl;
         break;
      } else {
         lambda_mass.push_back(lamcand.mass());
         h_lambda_mass->Fill(lamcand.mass());

         double three_mag = pow((lamcand.px()*lamcand.px()) + (lamcand.py()*lamcand.py()) + (lamcand.pz()*lamcand.pz()),0.5);
         double delta_x = pow(pow(lamcand.position().x()-pvXs[0],2) + pow(lamcand.position().y()-pvYs[0],2) + pow(lamcand.position().z()-pvZs[0],2),0.5); //had pvXs[1]?
         lamdistances.push_back(delta_x);
         h_lamdistances->Fill(delta_x);
         double tau = delta_x*pow(10,-2)*lamcand.mass()*pow(three_mag,-1)*pow(2.998*pow(10,8),-1);
         std::cout << "Lifetime : " << tau << std::endl;
         lamlifetimes.push_back(tau);
         h_lamlifetimes->Fill(tau);
      }

   }

   //try to find two pions
   for(const auto& pfcand1 : iEvent.get(pfcands_) ) {

      if( pfcand1.pdgId() == 310 ) {
         kscounter++;
         pfNeutralHadron_mass.push_back(pfcand1.mass());
         h_pfNeutralHadron_mass->Fill(pfcand1.mass());
      }
   }

   tree_->Fill();



#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   // if the SetupData is always needed
//   auto setup = iSetup.getData(setupToken_);
   // if need the ESHandle to check if the SetupData was there or not
//   auto pSetup = iSetup.getHandle(setupToken_);
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   std::cout << "No. infs : " << count_infs << std::endl;
   std::cout << "No. NaNs : " << count_nans << std::endl;
   std::cout << "-------------------\nNo. K-short: " << kscounter << std::endl;
   std::cout << "No. other: " << mctruth_elseNum << std::endl;
}




// ------------ method called once each job just before starting event loop  ------------
void
KShortLifetimeAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
KShortLifetimeAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
KShortLifetimeAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(KShortLifetimeAnalyzer);
