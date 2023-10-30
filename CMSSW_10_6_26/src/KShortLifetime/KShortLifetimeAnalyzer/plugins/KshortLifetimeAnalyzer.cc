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

      edm::Service<TFileService> fs_;
      TTree* tree_;
//      TTree* s_tree_;
      std::vector<float> kshort_mass;
      TH1D* h_kshort_mass;

      std::vector<float> pfNeutralHadron_mass;
      TH1D* h_pfNeutralHadron_mass;

      std::vector<float> dipion_mass;
      TH1D* h_dipion_mass;

      std::vector<float> kslifetimes;
      TH1D* h_kslifetimes;

      std::vector<float> cycrads;
      TH1D* h_cycrads;

      std::vector<float> ksdistances;
      TH1D* h_ksdistances;

//      std::vector<float> pvdistances;
//      TH1D* h_pvdistances;

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
  //tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))
//  s_tree_(fs_->make<TTree>("scatterTree","Scatter Plot Data")),
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

//   tree_->Branch("ks_distVp",&ks_distVp);
   h_ks_distVp = fs_->make<TH2F>("h_ks_distVp","h_ks_distVp",30,0,15,40,0,20);

//   s_tree_ = fs_->make<TTree>("s_tree_","s_tree_"); //one possible thing to come back to: does the "" name have to be diff?
//   s_tree_->Branch("ks_distVp",&ks_distVp);
//   s_ks_distVp = fs_->make<TGraph>; //am I allowed to do this?

   tree_->Branch("pfNeutralHadron_mass", &pfNeutralHadron_mass);
   h_pfNeutralHadron_mass = fs_->make<TH1D>("h_pfNeutralHadron_mass", "h_pfNeutralHadron_mass", 100, 0.4, 0.6);

   tree_->Branch("dipion_mass", &dipion_mass);
   h_dipion_mass = fs_->make<TH1D>("h_dipion_mass", "h_dipion_mass", 200, 0, 2);

   tree_->Branch("kslifetimes", &kslifetimes);
   h_kslifetimes = fs_->make<TH1D>("h_kslifetimes", "h_kslifetimes", 30, 0, 8*pow(10,-10));

   tree_->Branch("cycrads",&cycrads);
   h_cycrads = fs_->make<TH1D>("h_cycrads", "h_cycrads", 35, 0, 15);

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
  
//   double ks_distVp_x, ks_distVp_y; 
//   s_tree_->Branch("ks_distVp_x",&ks_distVp_x);
//   s_tree_->Branch("ks_distVp_y",&ks_distVp_y);
}


KShortLifetimeAnalyzer::~KShortLifetimeAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

float getCen_x(float p_planar_mag, float p_phiAtVtx, float r_xAtVtx, float r_yAtVtx, float mass, float p_thetaAtVtx, float p_mag, float charge) {
      float p_tAtVtx = p_planar_mag; //(since sin^2 + cos^2 = 1)
      std::cout << "p_tAtVtx : " << p_tAtVtx << std::endl;
      float p_xAtVtx = p_tAtVtx*cos(p_phiAtVtx);//convert momenta to cartesian coords [GeV/c]
      float p_yAtVtx = p_tAtVtx*sin(p_phiAtVtx); // [GeV/c]
      float B = 3.8; //strength of magnetic field in Teslas (T)
      float p_zAtVtx = p_tAtVtx*pow( tan(p_thetaAtVtx) , -1); //p_z in GeV/c
      std::cout << "p_zAtVtx : " << p_zAtVtx << std::endl;
//      float beta = p_mag*pow( pow(p_mag,2) + pow( mass*2.998*pow(10,8) ,2) , -0.5); //speed in units of v
//      std::cout << "beta : " << beta << std::endl;
      float pointer_x = p_yAtVtx*B;//Pointer -> vector to determine where center of cyclotron orbit is (equal to lorentz force vector up to some factors of charge, c, mass, and gamma)
      float pointer_y = -1*p_xAtVtx*B; //(these may be got from taking a determinant)
      float pointer_mag = pow( pow(pointer_x,2) + pow(pointer_y,2) , 0.5);
      float v_planar = p_tAtVtx*pow(p_mag,-1)*1; //planar velocity using v ~ c approx [in units of c]
      float cycrad_strangeunits = mass*v_planar*pow( charge*B ,-1); //cyclotron radius [in GeV/(c*e*T)]
      float cycrad = cycrad_strangeunits * 3.336 * pow(10,2); //cyclotron radius in centimeters, using fact that 1 GeV/(ceT) = 3.336 m.
      std::cout << "cycrad : " << cycrad << " cm" << std::endl;
      float toCenter_x = cycrad*pointer_x*pow(pointer_mag,-1); //x-vector leading from Vtx to center of cyclotron orbit (note that strange pointer units divide away)
      float toCenter_y = cycrad*pointer_y*pow(pointer_mag,-1); //y-vector leading from Vtx to center of cyclotron orbit
      float cen_x = r_xAtVtx + toCenter_x; //center x coord
      float cen_y = r_yAtVtx + toCenter_y; //center y coord
      //In the x,y-plane, all points on the particle's trajectory lie *cycrad away from the coordinates (*cen_x, *cen_y)
      return cen_x;
   }

float getCycrad(float p_planar_mag, float p_phiAtVtx, float r_xAtVtx, float r_yAtVtx, float mass, float p_thetaAtVtx, float p_mag, float charge) {
      float p_tAtVtx = p_planar_mag; //(since sin^2 + cos^2 = 1)
      std::cout << "p_tAtVtx : " << p_tAtVtx << std::endl;
      float p_xAtVtx = p_tAtVtx*cos(p_phiAtVtx);//convert momenta to cartesian coords [GeV/c]
      float p_yAtVtx = p_tAtVtx*sin(p_phiAtVtx); // [GeV/c]
      float B = 3.8; //strength of magnetic field in Teslas (T)
      float p_zAtVtx = p_tAtVtx*pow( tan(p_thetaAtVtx) , -1); //p_z in GeV/c
      std::cout << "p_zAtVtx : " << p_zAtVtx << std::endl;
//      float beta = p_mag*pow( pow(p_mag,2) + pow( mass*2.998*pow(10,8) ,2) , -0.5); //speed in units of v
//      std::cout << "beta : " << beta << std::endl;
      float pointer_x = p_yAtVtx*B;//Pointer -> vector to determine where center of cyclotron orbit is (equal to lorentz force vector up to some factors of charge, c, mass, and gamma)
      float pointer_y = -1*p_xAtVtx*B; //(these may be got from taking a determinant)
      float pointer_mag = pow( pow(pointer_x,2) + pow(pointer_y,2) , 0.5);
      float v_planar = p_tAtVtx*pow(p_mag,-1)*1; //planar velocity using v ~ c approx [in units of c]
      float cycrad_strangeunits = mass*v_planar*pow( charge*B ,-1); //cyclotron radius [in GeV/(c*e*T)]
      float cycrad = cycrad_strangeunits * 3.336 * pow(10,2); //cyclotron radius in centimeters, using fact that 1 GeV/(ceT) = 3.336 m.
      return cycrad;
   }

float getCen_y(float p_planar_mag, float p_phiAtVtx, float r_xAtVtx, float r_yAtVtx, float mass, float p_thetaAtVtx, float p_mag, float charge) {
      float p_tAtVtx = p_planar_mag; //(since sin^2 + cos^2 = 1)
      std::cout << "p_tAtVtx : " << p_tAtVtx << std::endl;
      float p_xAtVtx = p_tAtVtx*cos(p_phiAtVtx);//convert momenta to cartesian coords [GeV/c]
      float p_yAtVtx = p_tAtVtx*sin(p_phiAtVtx); // [GeV/c]
      float B = 3.8; //strength of magnetic field in Teslas (T)
      float p_zAtVtx = p_tAtVtx*pow( tan(p_thetaAtVtx) , -1); //p_z in GeV/c
      std::cout << "p_zAtVtx : " << p_zAtVtx << std::endl;
//      float beta = p_mag*pow( pow(p_mag,2) + pow( mass*2.998*pow(10,8) ,2) , -0.5); //speed in units of v
//      std::cout << "beta : " << beta << std::endl;
      float pointer_x = p_yAtVtx*B;//Pointer -> vector to determine where center of cyclotron orbit is (equal to lorentz force vector up to some factors of charge, c, mass, and gamma)
      float pointer_y = -1*p_xAtVtx*B; //(these may be got from taking a determinant)
      float pointer_mag = pow( pow(pointer_x,2) + pow(pointer_y,2) , 0.5);
      float v_planar = p_tAtVtx*pow(p_mag,-1)*1; //planar velocity using v ~ c approx [in units of c]
      float cycrad_strangeunits = mass*v_planar*pow( charge*B ,-1); //cyclotron radius [in GeV/(c*e*T)]
      float cycrad = cycrad_strangeunits * 3.336 * pow(10,2); //cyclotron radius in centimeters, using fact that 1 GeV/(ceT) = 3.336 m.
      std::cout << "cycrad : " << cycrad << " cm" << std::endl;
      float toCenter_x = cycrad*pointer_x*pow(pointer_mag,-1); //x-vector leading from Vtx to center of cyclotron orbit (note that strange pointer units divide away)
      float toCenter_y = cycrad*pointer_y*pow(pointer_mag,-1); //y-vector leading from Vtx to center of cyclotron orbit
      float cen_x = r_xAtVtx + toCenter_x; //center x coord
      float cen_y = r_yAtVtx + toCenter_y; //center y coord
      std::cout << "cen_y : " << cen_y << std::endl;
      //In the x,y-plane, all points on the particle's trajectory lie *cycrad away from the coordinates (*cen_x, *cen_y)
      return cen_y;
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
//Begin edited bit-----------------
   kshort_mass.clear();
   dipion_mass.clear();
   kslifetimes.clear();
   cycrads.clear();
   ksdistances.clear();
//   pvdistances.clear();
   lambda_mass.clear();
   lamdistances.clear();
   lamlifetimes.clear();
   ks_bkgd_lifetimes.clear();
//   ks_distVp_x.clear();
//   ks_distVp_y.clear();
//   ks_distVp.clear();

//   s_ks_distVp = fs_->make<TGraph>("ks_distVp","ks_distVp");
//   TTree *s_tree_ = new TTree("scatterTree","Scatter Plot Data");

//   Double_t ks_distVp_x, ks_distVp_y;

//   s_tree_ = fs_->make<TTree>("s_tree","s_tree");

   std::cout << "Run: "<< iEvent.id().run() << " luminosity block: "<< iEvent.eventAuxiliary().luminosityBlock() <<" event: "<< iEvent.id().event() << std::endl;

//   for(const auto& secondaryVertex : iEvent.get(secondaryVertices_) ) {
//      float svX = secondaryVertex.position().x();
//      float svY = secondaryVertex.position().y();
//      float svZ = secondaryVertex.position().z();
//      std::cout << "SV @ (" << svX << ", " << svY << ", " << svZ << ")\nNo. Assoc. tracks : " << secondaryVertex.nTracks() << std::endl;
//   }

//   std::vector<float> massList;
//   massList.clear();

//   for(const auto& pfcand : iEvent.get(pfcands_) ) {
//      std::cout << "pfcand pt=" << pfcand.pt() << std::endl;
//      float vX = pfcand.vertex().x();
//      float vY = pfcand.vertex().y();
//      float vZ = pfcand.vertex().z();
//      std::cout << "Closest approach to PV at : (" << vX << ", " << vY << ", " << vZ << ")" << std::endl;
//      int count = 0;
//      std::cout<<"Charge : "<<pfcand.charge()<<std::endl;
//      if(pfcand.charge() == 1) { std::cout<<"Got one."<<std::endl; };
//      for(const auto& pfcand2 : iEvent.get(pfcands_) ) {
//         float vX2 = pfcand2.vertex().x();
//         float vY2 = pfcand2.vertex().y();
//         float vZ2 = pfcand2.vertex().z();
//         if(vX == vX2) {
//            if(vY == vY2) {
//               if(vZ == vZ2) {
//                  std::cout << "There exist two+ tracks with same SV." << std::endl;
//                  count++;
//                  if(count == 1) {
//                    std::cout<< "-----------------------\nThere exist exactly two tracks with same SV." << std::endl;
//                     if(pfcand.charge() == 1) {
//                        std::cout<<"\nThere is a +1 particle."<<std::endl;
//                        if(pfcand2.charge() == -1) {
//                           std::cout<<"************************\nThere exist two such track with appropriate signs."<<std::endl;
//                           TLorentzVector p1(pfcand.pt(), pfcand.eta(), pfcand.phi(), pfcand.energy());
//                           TLorentzVector p2(pfcand2.pt(), pfcand2.eta(), pfcand2.phi(), pfcand2.energy());
//                           TLorentzVector pParent = p1 + p2;
//                           float mParent = pParent.M();
//                           dipion_mass.push_back(mParent);
//                           h_dipion_mass->Fill(mParent);
//                           std::cout<<"Parent mass : "<<mParent<<std::endl;
//                           std::cout<<"Got one!"<<std::endl;
//                        };
//                     } else if(pfcand.charge() == -1) {
//                        std::cout<<"\nThere is a -1 particle."<<std::endl; 
//                        if(pfcand2.charge() == 1) {
//                           std::cout<<"************************\nThere exist two such track with appropriate signs."<<std::endl;
//                           TLorentzVector p1(pfcand.pt(), pfcand.eta(), pfcand.phi(), pfcand.energy());
//                           TLorentzVector p2(pfcand2.pt(), pfcand2.eta(), pfcand2.phi(), pfcand2.energy());
//                           TLorentzVector pParent = p1 + p2;
//                           float mParent = pParent.M();
//                           dipion_mass.push_back(mParent);
//                           h_dipion_mass->Fill(mParent);
//                           std::cout<<"Parent mass : "<<mParent<<std::endl;
//                           std::cout<<"Got one!"<<std::endl;
//                        };
//                     } else {
//                        break;
//                     }
//                  }
//               }
//            }
//         }
//      }
//   }

//   float getCenter(float p_planar_mag, float p_phiAtVtx, float r_xAtVtx, float r_yAtVtx, float mass, float p_thetaAtVtx, float p_mag, float charge) {
//      float p_tAtVtx = p_planar_mag; //(since sin^2 + cos^2 = 1)
//      float p_xAtVtx = p_tAtVtx*cos(p_phiAtVtx)*pow(2.998*pow(10,8),-1);//convert momenta to cartesian coords [GeV/c]
//      float p_yAtVtx = p_tAtVtx*sin(p_phiAtVtx)*pow(2.998*pow(10,8),-1); // [GeV/c]
//      float B = 3.8; //strength of magnetic field in Teslas (T)
//      float p_zAtVtx = p_tAtVtx*pow( p_thetaAtVtx , -1);
//      float beta = p_mag*2.998*pow(10,8)*pow( pow(p_mag,2) + pow( mass*2.998*pow(10,8) ,2) , -0.5); //speed in units of c
//      float gamma = 1*pow(1 - pow(beta,2),-0.5); //lorentz factor
//      float beta_xAtVtx = p_xAtVtx*pow(mass*gamma, -1); //velocity planar components [in c]
//      float beta_yAtVtx = p_yAtVtx*pow(mass*gamma, -1); //[in c]
//      float beta_magAtVtx = pow( pow(beta_xAtVtx,2) + pow(beta_yAtVtx,2) ,0.5)  //[in c]
//      float pointer_x = beta_yAtVtx*B;//Pointer -> vector to determine where center of cyclotron orbit is (equal to lorentz force vector up to some factors of charge and c)
//      float pointer_y = -1*beta_xAtVtx*B; //(these may be got from taking a determinant)
//      float pointer_mag = pow( pow(pointer_x,2) + pow(pointer_y,2) , 0.5);
//      float cycrad_strangeunits = mass*beta_matAtVtx*pow( charge*B ,-1); //cyclotron radius [in GeV/(c*e*T)]
//      float cycrad = cycrad_strangeunits * 3.336; //cyclotron radius in meters
//      float toCenter_x = cycrad*pointer_x*pow(pointer_mag,-1); //x-vector leading from Vtx to center of cyclotron orbit (note that strange pointer units divide away)
//      float toCenter_y = cycrad*pointer_y*pow(pointer_mag,-1); //y-vector leading from Vtx to center of cyclotron orbit
//      float cen_x = r_xAtVtx + toCenter_x; //center x coord
//      float cen_y = r_yAtVtx + toCenter_y; //center y coord
      //In the x,y-plane, all points on the particle's trajectory lie *cycrad away from the coordinates (*cen_x, *cen_y)
//      std::vector<float> cenCoords; //vector that contakins x, y coords of center
//      cenCoords.push_back(cen_x);
//      cenCoords.push_back(cen_y);
//      return cenCoords;
//   }

//   for(const auto& pfcand1 : iEvent.get(pfcands_) ) {
      //Testing
//      std::cout << "phiAtVtx : " << pfcand1.phiAtVtx() << std::endl;
      //std::cout << "pxAtVtx : " << pfcand1.pxAtVtx() << std::endl;
      //std::cout << "v at vtx : " << pfcand1.vxAtVtx() << std::endl;

      //Data analysis
      //Find components of 3-momentum at .vertex() (aka closest pont on trajectory to PV)
//      float p_planar_mag1 = pow( pow(pfcand1.px(),2) + pow(pfcand1.py(),2) , 0.5); //magnitude of planar (viz. non-z) component of 3-mom [GeV/c]
//      std::cout << "p_planar_mag : " << p_planar_mag << std::endl;
//      float p_phiAtVtx1 = pfcand1.phiAtVtx(); //momentum phi component [GeV/c]
//      float r_xAtVtx1 = pfcand1.vertex().x(); //get positions [cm]
//      float r_yAtVtx1 = pfcand1.vertex().y();
//      std::cout << "r_xAtVtx : " << r_xAtVtx << std::endl;
//      std::cout << "r_yAtVtx : " << r_yAtVtx << std::endl;
//      float mass1 = pfcand1.mass(); //mass in units of GeV/c^2
//      float p_thetaAtVtx1 = 2*atan(exp(-1*pfcand1.etaAtVtx())); //[Gev/c]
//      float p_mag1 = pow( pow(p_planar_mag1,2) + pow(pfcand1.pz(),2) , 0.5); //magnitude of 3-momentum in GeV/c
//      float charge1 = pfcand1.charge(); //charge in elementary charges (e)

      //First pass: approximate the trajectory to a cylinder of cyclotron radius and infinite length in z-direction, viz. it suffices for their paths to cross in the x,y-plane at any point in time or along the z-dimension
//      if (charge1 == 1 || -1) {
//         float x1 = getCen_x( p_planar_mag1, p_phiAtVtx1, r_xAtVtx1, r_yAtVtx1, mass1, p_thetaAtVtx1, p_mag1, charge1 );
//         float y1 = getCen_y( p_planar_mag1, p_phiAtVtx1, r_xAtVtx1, r_yAtVtx1, mass1, p_thetaAtVtx1, p_mag1, charge1 );
//         float cycrad1 = getCycrad( p_planar_mag1, p_phiAtVtx1, r_xAtVtx1, r_yAtVtx1, mass1, p_thetaAtVtx1, p_mag1, charge1 );
//         if (isinf(cycrad1)) {
//            std::cout << "Infinite." << std::endl;
//            count_infs++;
//         } else if (isnan(cycrad1)) {
//            std::cout << "NaN." << std::endl;
//            count_nans++;
//         } else if (isnormal(cycrad1)) {
//            cycrads.push_back(cycrad1);
//            h_cycrads->Fill(cycrad1);
//         }
//         std::cout << "Center coordinates:  (" << x1 << ", " << y1 << ")" <<std::endl;
//      } else { std::cout << "Failed charge criteria." << std::endl; }

//      for(const auto& pfcand2 : iEvent.get(pfcands_) ) {
      //Data analysis
      //Find components of 3-momentum at .vertex() (aka closest pont on trajectory to PV)
//         float p_planar_mag2 = pow( pow(pfcand2.px(),2) + pow(pfcand2.py(),2) , 0.5); //magnitude of planar (viz. non-z) component of 3-mom [GeV/c]
//         float p_phiAtVtx2 = pfcand2.phiAtVtx(); //momentum phi component [GeV/c]
//         float r_xAtVtx2 = pfcand2.vertex().x(); //get positions [cm]
//         float r_yAtVtx2 = pfcand2.vertex().y();
//         float mass2 = pfcand2.mass(); //mass in units of GeV/c^2
//         float p_thetaAtVtx2 = 2*atan(exp(-1*pfcand2.etaAtVtx())); //[Gev/c]
//         float p_mag2 = pow( pow(p_planar_mag2,2) + pow(pfcand2.pz(),2) , 0.5); //magnitude of 3-momentum in GeV/c
//         float charge2 = pfcand2.charge(); //charge in elementary charges (e)
//
      //First pass: approximate the trajectory to a cylinder of cyclotron radius and infinite length in z-direction, viz. it suffices for their paths to cross in the x,y-plane at any point in time or along the z-dimension
//         if (charge2 == 1 || -1) {
//            float x2 = getCen_x( p_planar_mag2, p_phiAtVtx2, r_xAtVtx2, r_yAtVtx2, mass2, p_thetaAtVtx2, p_mag2, charge2 );
//            float y2 = getCen_y( p_planar_mag2, p_phiAtVtx2, r_xAtVtx2, r_yAtVtx2, mass2, p_thetaAtVtx2, p_mag2, charge2 );
//            float cycrad2 = getCycrad( p_planar_mag2, p_phiAtVtx2, r_xAtVtx2, r_yAtVtx2, mass2, p_thetaAtVtx2, p_mag2, charge2 );

//            std::cout << "Center coordinates:  (" << x2 << ", " << y2 << ")" <<std::endl;
//         } else { std::cout << "Failed charge criteria." << std::endl; }
      
//   }

//   }

   std::vector<float> pvXs;
   std::vector<float> pvYs;
   std::vector<float> pvZs;
   for(const auto& primaryVertex : iEvent.get(primaryVertices_) ) {
      pvXs.push_back(primaryVertex.position().x());
      pvYs.push_back(primaryVertex.position().y());
      pvZs.push_back(primaryVertex.position().z());
   }

//   int index = 0;

   unsigned counter = 0;
   for(const auto& kscand : iEvent.get(kscands_) ) {
      std::cout << "KShort Candidate pt=" << kscand.pt() << " | mass=" << kscand.mass() << " | No. Daughters: " << kscand.numberOfDaughters() << " | vertex position (x,y,z): " << kscand.position().x() << ", " << kscand.position().y() << ", " << kscand.position().z() << std::endl;
      counter++;
      for (std::size_t i = 0; i < kscand.numberOfDaughters(); i = i + 1 ) {
         std::cout << "    Mass of daughter #" << i << " is " << kscand.daughter(i)->mass() << std::endl;
         std::cout << "    Charge of daughter #" << i << " is " << kscand.daughter(i)->charge()  << std::endl;
      }
      //mass cuts: ignore candidates w/ masses >1 sigma from mean
      if (kscand.mass() < 0.495439) {
         break;
      } else if (kscand.mass() > 0.501482 ) {
         break;
      } else {
         kshort_mass.push_back(kscand.mass());
         h_kshort_mass->Fill(kscand.mass());

         std::cout << "ID: " << kscand.pdgId() << std::endl;

         if (kscand.pdgId()==0) {
            mctruth_ksNum++;
         } else {
            mctruth_elseNum++;
         }

        //distance of No. 1 primary vertex from origin

        double pvdist = pow(pow(pvXs[0],2) + pow(pvYs[0],2) + pow(pvZs[0],2),0.5);
        std::cout << "PV DIST : " << pvdist << std::endl;
//        pvdistances.push_back(pvdist);
//        h_pvdistances->Fill(pvdist);

         double three_mag = pow((kscand.px()*kscand.px()) + (kscand.py()*kscand.py()) + (kscand.pz()*kscand.pz()),0.5);
         double delta_x = pow(pow(kscand.position().x()-pvXs[0],2) + pow(kscand.position().y()-pvYs[0],2) + pow(kscand.position().z()-pvZs[0],2),0.5); //had pvXs[1]?
         ksdistances.push_back(delta_x);
         h_ksdistances->Fill(delta_x);

         ks_distVp_x.push_back(three_mag);
         ks_distVp_y.push_back(delta_x);

         h_ks_distVp->Fill(three_mag, delta_x);

//         double ks_distVp_x, ks_distVp_y;

//         ks_distVp_x = three_mag; s_tree_->Fill();
//         ks_distVp_y = delta_x; s_tree_->Fill();

         //add points to scatterplot of distance traveled before decay vs (3-)momentum magnitude
         //s_ks_distVp->SetPoint(index, 1.0, 1.0); //s_ks_distVp->GetN(), three_mag, delta_x

//         double tau = delta_x*pow(10,-2)*kscand.mass()*pow(three_mag,-1)*pow(2.998*pow(10,8),-1);

         double tau = delta_x*pow(10,-2)*kscand.mass()*pow(kscand.energy(), -1)*pow(2.998e8, -1);
         std::cout << "Lifetime : " << tau << std::endl;
         kslifetimes.push_back(tau);
         h_kslifetimes->Fill(tau);

         //Fill the background lifetimes with lifetime of candidates with masses more than 5 sigma of fit mean (cf. display.cxx)
         if (kscand.mass() < 0.468245 or kscand.mass() > 0.528676) {
            ks_bkgd_lifetimes.push_back(tau);
            h_ks_bkgd_lifetimes->Fill(tau);
         }
//         index++;
//	 std::cout<<"Index : "<<index<<std::endl;

      }
   }

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
//      for(const auto& pfcand2 : iEvent.get(pfcands_) ) {
//         if( pfcand1.charge() == 1)  {
//            if( pfcand2.charge() != -1 ) {
//               break;
//            } else {
//               double mSum = pfcand1.mass() + pfcand2.mass();
//               dipion_mass.push_back(mSum);
//               h_dipion_mass->Fill(mSum);
//            }
//         } else if( pfcand1.charge() == -1 ) {
//             if( pfcand2.charge() != 1 ) {
//               break;
//            } else {
//               double mSum = pfcand1.mass() + pfcand2.mass();
//               dipion_mass.push_back(mSum);
//               h_dipion_mass->Fill(mSum);
//            }           
//         }
//      }
      if( pfcand1.pdgId() == 310 ) {
         kscounter++;
         pfNeutralHadron_mass.push_back(pfcand1.mass());
         h_pfNeutralHadron_mass->Fill(pfcand1.mass());
      }
   }

//   unsigned iVertex = 0;
//   for(const auto& primaryVertex : iEvent.get(primaryVertices_) ) {
//      std::cout << "primary vertex["<<iVertex<<"]: (x,y,z): " <<primaryVertex.position().x()<<", "<<primaryVertex.position().y()<<", "<<primaryVertex.position().z()<<std::endl;
//      iVertex++;
//   }


   tree_->Fill();

//   s_tree_->Write();

//End edited bit------------------

//   for(const auto& track : iEvent.get(tracksToken_) ) {
      // do something with track parameters, e.g, plot the charge.
      // int charge = track.charge();
//      std::cout << "track pt: " << track.pt() << std::endl;
//   }

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
