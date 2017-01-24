// -*- C++ -*-
//
// Package:    EgammaWork/ElectronNtupler
// Class:      SimpleElectronNtupler
// 
/**\class SimpleElectronNtupler SimpleElectronNtupler.cc EgammaWork/ElectronNtupler/plugins/SimpleElectronNtupler.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Ilya Kravchenko
//         Created:  Thu, 10 Jul 2014 09:54:13 GMT
//         A.Juodagalvis modified to derive smearing variance, Sept 20, 2016
//
//


// system include files
#include <memory>
#include <vector>
#include <regex>

#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "Math/VectorUtil.h"
#include "TClonesArray.h"

#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefHolder.h"
#include "DataFormats/Common/interface/RefVectorHolder.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "/afs/cern.ch/work/r/rchawla/public/utils.h"

#include "RecoEgamma/EgammaIsolationAlgos/interface/HcalPFClusterIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EcalPFClusterIsolation.h"

// To produce smearing variance
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "EgammaAnalysis/ElectronTools/interface/SimpleElectron.h"
#include "EgammaAnalysis/ElectronTools/interface/EpCombinationTool.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h"
#include "EgammaAnalysis/ElectronTools/src/EnergyScaleCorrection_class.cc"
#include <TRandom3.h>
#include <TH1F.h>

//
// class declaration
//

class SimpleElectronNtupler : public edm::EDAnalyzer {
  public:
    explicit SimpleElectronNtupler(const edm::ParameterSet&);
    ~SimpleElectronNtupler();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    enum ElectronMatchType {UNMATCHED = 0, 
      TRUE_PROMPT_ELECTRON, 
      TRUE_ELECTRON_FROM_TAU,
      TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    void printCutFlowResult(vid::CutFlowResult &cutflow);

    // ----------member data ---------------------------
    // Data members that are the same for AOD and miniAOD
    edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales>triggerPrescale_; 
    edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupToken_;
    edm::EDGetTokenT<double> rhoToken_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    edm::EDGetTokenT<pat::METCollection> metToken_;

    // AOD case data members
    edm::EDGetToken electronsToken_;
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
    edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;

    // MiniAOD case data members
    edm::EDGetToken muonsMiniAODToken_;
    edm::EDGetToken electronsMiniAODToken_;
    edm::EDGetToken photonsMiniAODToken_;
    edm::EDGetTokenT<reco::VertexCollection> vtxMiniAODToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;
    edm::EDGetTokenT<reco::ConversionCollection> conversionsMiniAODToken_;

    // ID decisions objects
    edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;
    //edm::EDGetTokenT<edm::ValueMap<bool> > eletrigMVAlooseMapToken_;
    //edm::EDGetTokenT<edm::ValueMap<bool> > eletrigMVAtightMapToken_;    
    //edm::EDGetTokenT<edm::ValueMap<bool> > elenontrigMVAlooseMapToken_;
    //edm::EDGetTokenT<edm::ValueMap<bool> > elenontrigMVAtightMapToken_;

    edm::EDGetTokenT<GenEventInfoProduct> genToken_;

    // One example of full information about the cut flow
    edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleMediumIdFullInfoMapToken_;

    // Verbose output for ID
    bool verboseIdFlag_;

    //
    edm::Service<TFileService> fs;
    TTree *electronTree_;

    //edm::LumiReWeighting LumiWeights_;

    // If MC
    bool misMC;
    bool misSIG;


    // Weights for MC@NLO
    double theWeight;

    // General
    Double_t RunNo, EvtNo, Lumi, Bunch;

    // Vars for PVs
    Int_t pvNTracks;

    // Vars for pile-up
    Int_t nPUTrue;    // true pile-up
    Int_t nPU;        // generated pile-up
    Int_t nPV;        // number of reconsrtucted primary vertices
    Float_t rho;      // the rho variable
    //Double_t PUWeight_;

    // Trigger
    std::regex ele23_WPLoose;
    std::regex ele27_WP85;
    std::regex isoMu20;
    std::regex ele17_ele12;
    std::regex mu8_ele17;
    std::regex mu17_ele12;
    std::regex mu8_ele23;
    
    bool Ele23_WPLoose;
    bool Ele27_WP85;
    bool IsoMu20;
    bool Ele17_Ele12;
    bool Mu8_Ele17;
    bool Mu17_Ele12;
    bool Mu8_Ele23;

    // tau variables
    Int_t tauFlag;
    Int_t nGenTaus;

    // all electron variables
    Int_t nGenElectrons;

    //std::vector<Int_t>   genPhoton_Id_;
    //std::vector<Int_t>   genlepMother_Id_;
    //std::vector<Bool_t>  fromHProcessFinalState_;
    //std::vector<Bool_t>  fromHProcessDecayed_;

    std::vector<Float_t> genPhoton_Px_;
    std::vector<Float_t> genPhoton_Py_;
    std::vector<Float_t> genPhoton_Pz_;
    std::vector<Float_t> genPhoton_Pt_;
    std::vector<Float_t> genPhoton_Eta_;
    std::vector<Float_t> genPhoton_Rap_;
    std::vector<Float_t> genPhoton_Phi_;
    std::vector<Float_t> genPhoton_En_;

    std::vector<Float_t> genPostFSR_Px_;
    std::vector<Float_t> genPostFSR_Py_;
    std::vector<Float_t> genPostFSR_Pz_;
    std::vector<Float_t> genPostFSR_Pt_;
    std::vector<Float_t> genPostFSR_Eta_;
    std::vector<Float_t> genPostFSR_Rap_;
    std::vector<Float_t> genPostFSR_Phi_;
    std::vector<Float_t> genPostFSR_En_;

    std::vector<Float_t> genPreFSR_Px_;
    std::vector<Float_t> genPreFSR_Py_;
    std::vector<Float_t> genPreFSR_Pz_;
    std::vector<Float_t> genPreFSR_Pt_;
    std::vector<Float_t> genPreFSR_Eta_;
    std::vector<Float_t> genPreFSR_Rap_;
    std::vector<Float_t> genPreFSR_Phi_;
    std::vector<Float_t> genPreFSR_En_;

    Int_t nElectrons;

    std::vector<Float_t> ptElec_;
    std::vector<Float_t> etaElec_;
    std::vector<Float_t> rapElec_;
    std::vector<Float_t> phiElec_;
    std::vector<Float_t> energyElec_;
    std::vector<Float_t> massElec_;
    std::vector<Float_t> chargeElec_;

    std::vector<Float_t> error_scale_;
    std::vector<Float_t> sigma_;
    std::vector<Float_t> sigma1_;
    std::vector<Float_t> sigma2_;
    std::vector<Float_t> sigma3_;
    std::vector<Float_t> sigma4_;
    
    std::vector<Float_t> enSC_;
    std::vector<Float_t> preEnSC_;
    std::vector<Float_t> rawEnSC_;
    std::vector<Float_t> etSC_;
    std::vector<Float_t> etaSC_;
    std::vector<Float_t> phiSC_;

    std::vector<Float_t> full5x5_sigmaIetaIeta_;
    std::vector<Float_t> E1x5_;
    std::vector<Float_t> E2x5_;
    std::vector<Float_t> E5x5_;
    std::vector<Float_t> hOverE_;
    std::vector<Float_t> etaScWidth_;
    std::vector<Float_t> phiScWidth_;
    std::vector<Float_t> r9_;

    std::vector<Float_t> dEtaIn_;
    std::vector<Float_t> dPhiIn_;
    std::vector<Float_t> isoChargedHadrons_;
    std::vector<Float_t> isoNeutralHadrons_;
    std::vector<Float_t> isoPhotons_;
    std::vector<Float_t> isoChargedFromPU_;
    std::vector<Float_t> isoDeltaBeta_;
    std::vector<Float_t> isoRho_;
    std::vector<Float_t> ooEmooP_;
    std::vector<Float_t> d0_;
    std::vector<Float_t> dz_;
    std::vector<Int_t>   expectedMissingInnerHits_;
    std::vector<Int_t>   passConversionVeto_;
    std::vector<Float_t> brem_;

    std::vector<Float_t> eleInBarrel_;
    std::vector<Float_t> eleInEndcap_;

    std::vector<Int_t> passVetoId_;
    std::vector<Int_t> passLooseId_;
    std::vector<Int_t> passMediumId_;
    std::vector<Int_t> passTightId_;
    std::vector<Int_t> passHEEPId_;
    std::vector<Int_t> eleEcalDrivenSeed_;

    std::vector<Int_t> isPassMedium_NoPt_;
    std::vector<Int_t> isPassMedium_NoScEta_; 
    std::vector<Int_t> isPassMedium_NoDEta_;
    std::vector<Int_t> isPassMedium_NoDPhi_;
    std::vector<Int_t> isPassMedium_NoSigmaEtaEta_;
    std::vector<Int_t> isPassMedium_NoHOverE_;
    std::vector<Int_t> isPassMedium_NoDxy_;
    std::vector<Int_t> isPassMedium_NoDz_;
    std::vector<Int_t> isPassMedium_NoEInvP_;
    std::vector<Int_t> isPassMedium_NoPFIso_;
    std::vector<Int_t> isPassMedium_NoConVeto_;
    std::vector<Int_t> isPassMedium_NoMissHits_;

    // Trigger objects
    std::vector<double> pt_Ele23;
    std::vector<double> eta_Ele23;
    std::vector<double> phi_Ele23;

    // Triggers for Fake-Rate method
    std::regex Photon_22;
    std::regex Photon_30;
    std::regex Photon_36;
    std::regex Photon_50;
    std::regex Photon_75;
    std::regex Photon_90;
    std::regex Photon_120;
    std::regex Photon_175;

    int singlePhoton;
    int prescalePhoton;

    bool photon22;
    bool photon30;
    bool photon36;
    bool photon50;
    bool photon75;
    bool photon90;
    bool photon120;
    bool photon175;

    std::vector<double> etSPhoHLT_;
    std::vector<double> ptSPhoHLT_;
    std::vector<double> etaSPhoHLT_;
    std::vector<double> phiSPhoHLT_;

    // all muon variables
    Int_t nMuons;

    std::vector<bool>   isLooseMuon_;
    std::vector<bool>   isTightMuon_;
    std::vector<bool>   isHighPtMuon_;

    std::vector<double> ptMuon_;
    std::vector<double> etaMuon_;
    std::vector<double> phiMuon_;
    std::vector<double> energyMuon_;
    std::vector<double> chargeMuon_;

    std::vector<double> isoChargedHadronPfR04Muon_;
    std::vector<double> isoNeutralHadronPfR04Muon_;
    std::vector<double> isoGammaPfR04Muon_;
    std::vector<double> isoChargedFromPUMuon_;
    std::vector<double> isoPFMuon_;
    std::vector<double> isoTrkMuon_;

    // MET variables
    std::vector<double> metPt_;
    std::vector<double> metPhi_;
    std::vector<double> metSumEt_;

    // photon variables
    Int_t nPhotons;

    std::vector<Float_t> ptPhoton_;
    std::vector<Float_t> etaPhoton_;
    std::vector<Float_t> phiPhoton_;

    double DeltaR(const pat::Electron& e, std::vector<pat::TriggerObjectStandAlone> object);

};

SimpleElectronNtupler::SimpleElectronNtupler(const edm::ParameterSet& iConfig):

  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  eleHEEPIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap"))),
  eleMediumIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleMediumIdFullInfoMap"))),
  verboseIdFlag_(iConfig.getParameter<bool>("eleIdVerbose"))
{

  ele23_WPLoose = "(HLT_Ele23_WPLoose_Gsf_v)(.*)";
  ele27_WP85    = "(HLT_Ele27_WP85_Gsf_v)(.*)";
  isoMu20       = "(HLT_IsoMu20_v)(.*)";
  ele17_ele12   = "(HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v)(.*)";
  mu8_ele17     = "(HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v)(.*)";
  mu17_ele12    = "(HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v)(.*)";
  mu8_ele23     = "(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v)(.*)";

  // triggers for Fake-Rate method
  Photon_22  = "(HLT_Photon22_v)(.*)";
  Photon_30  = "(HLT_Photon30_v)(.*)";
  Photon_36  = "(HLT_Photon36_v)(.*)";
  Photon_50  = "(HLT_Photon50_v)(.*)";
  Photon_75  = "(HLT_Photon75_v)(.*)";
  Photon_90  = "(HLT_Photon90_v)(.*)";
  Photon_120 = "(HLT_Photon120_v)(.*)";
  Photon_175 = "(HLT_Photon175_v)(.*)";

  misMC = iConfig.getUntrackedParameter<bool>("isMC");
  misSIG  = iConfig.getUntrackedParameter<bool>("isSIG");

  // initialize the random number generator

  // Prepare tokens for all input collections and objects

  // Trigger
  triggerToken_  = mayConsume<edm::TriggerResults>
    (iConfig.getParameter<edm::InputTag>
     ("trigger"));

  triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>
    (iConfig.getParameter<edm::InputTag>
     ("objects"));

  triggerPrescale_ = consumes<pat::PackedTriggerPrescales>
    (iConfig.getParameter<edm::InputTag>
     ("prescale"));

  // Universal tokens for AOD and miniAOD  
  pileupToken_ = consumes<edm::View<PileupSummaryInfo> >
    (iConfig.getParameter <edm::InputTag>
     ("pileup"));

  rhoToken_    = consumes<double> 
    (iConfig.getParameter <edm::InputTag>
     ("rho"));

  beamSpotToken_    = consumes<reco::BeamSpot> 
    (iConfig.getParameter <edm::InputTag>
     ("beamSpot"));

  // AOD tokens
  electronsToken_    = mayConsume<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
     ("electrons"));

  vtxToken_          = mayConsume<reco::VertexCollection>
    (iConfig.getParameter<edm::InputTag>
     ("vertices"));

  genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticles"));

  conversionsToken_ = mayConsume< reco::ConversionCollection >
    (iConfig.getParameter<edm::InputTag>
     ("conversions"));

  // MiniAOD tokens
  // For electrons, use the fact that pat::Electron can be cast into GsfElectron
  muonsMiniAODToken_        = mayConsume<edm::View<pat::Muon> >
    (iConfig.getParameter<edm::InputTag>
     ("muonsMiniAOD"));

  electronsMiniAODToken_    = mayConsume<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
     ("electronsMiniAOD"));

  photonsMiniAODToken_ = mayConsume<edm::View<reco::Photon> >
    (iConfig.getParameter<edm::InputTag>
     ("photonsMiniAOD"));

  vtxMiniAODToken_          = mayConsume<reco::VertexCollection>
    (iConfig.getParameter<edm::InputTag>
     ("verticesMiniAOD"));

  genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticlesMiniAOD"));

  conversionsMiniAODToken_ = mayConsume< reco::ConversionCollection >
    (iConfig.getParameter<edm::InputTag>
     ("conversionsMiniAOD"));

  metToken_ = mayConsume<pat::METCollection>
    (iConfig.getParameter<edm::InputTag>
     ("metsMiniAOD"));

  genToken_ = mayConsume<GenEventInfoProduct>
  (iConfig.getParameter<edm::InputTag>("eventWeight"));

  //weightSrcToken_ = iC.consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("eventWeight"));

  //
  // Set up the ntuple structure
  //
  edm::Service<TFileService> fs;
  electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");

  electronTree_->Branch("RunNo", &RunNo, "RunNo/D");
  electronTree_->Branch("EvtNo", &EvtNo, "EvtNo/D");
  electronTree_->Branch("Lumi", &Lumi, "Lumi/D");
  electronTree_->Branch("Bunch", &Bunch, "Bunch/D");

  electronTree_->Branch("theWeight", &theWeight, "theWeight/D");

  electronTree_->Branch("pvNTracks"    ,  &pvNTracks , "pvNTracks/I");

  electronTree_->Branch("nPV"        ,  &nPV     , "nPV/I");
  electronTree_->Branch("nPU"        ,  &nPU     , "nPU/I");
  electronTree_->Branch("nPUTrue"    ,  &nPUTrue , "nPUTrue/I");
  electronTree_->Branch("rho"        ,  &rho , "rho/F");
  //electronTree_->Branch("PUWeight", &PUWeight_, "PUWeight/D");

  electronTree_->Branch("Ele23_WPLoose" ,  &Ele23_WPLoose);
  electronTree_->Branch("Ele27_WP85"    ,  &Ele27_WP85);
  electronTree_->Branch("IsoMu20"       ,  &IsoMu20);
  electronTree_->Branch("Ele17_Ele12"   ,  &Ele17_Ele12);
  electronTree_->Branch("Mu8_Ele17"     ,  &Mu8_Ele17);
  electronTree_->Branch("Mu17_Ele12"    ,  &Mu17_Ele12);
  electronTree_->Branch("Mu8_Ele23"     ,  &Mu8_Ele23);

  electronTree_->Branch("etSPhoHLT", &etSPhoHLT_);
  electronTree_->Branch("ptSPhoHLT", &ptSPhoHLT_);
  electronTree_->Branch("etaSPhoHLT", &etaSPhoHLT_);
  electronTree_->Branch("phiSPhoHLT", &phiSPhoHLT_);
  electronTree_->Branch("singlePhoton", &singlePhoton);
  electronTree_->Branch("prescalePhoton", &prescalePhoton);

  electronTree_->Branch("pt_Ele23"    ,  &pt_Ele23    );
  electronTree_->Branch("eta_Ele23"   ,  &eta_Ele23   );
  electronTree_->Branch("phi_Ele23"   ,  &phi_Ele23   );

  electronTree_->Branch("nEle"    ,  &nElectrons , "nEle/I");
  electronTree_->Branch("nGenEle"    ,  &nGenElectrons , "nGenEle/I");
  electronTree_->Branch("nGenTau"    ,  &nGenTaus , "nGenTau/I");
  electronTree_->Branch("tauFlag", &tauFlag, "tauFlag/I");

  //electronTree_->Branch("genlep_Id", &genlep_Id_);
  //electronTree_->Branch("genlepMother_Id", &genlepMother_Id_);
  //electronTree_->Branch("fromHProcessFinalState", &fromHProcessFinalState_);
  //electronTree_->Branch("fromHProcessDecayed", &fromHProcessDecayed_);

  electronTree_->Branch("genPhoton_Px"    ,  &genPhoton_Px_    );
  electronTree_->Branch("genPhoton_Py"    ,  &genPhoton_Py_    );
  electronTree_->Branch("genPhoton_Pz"    ,  &genPhoton_Pz_    );
  electronTree_->Branch("genPhoton_Pt"    ,  &genPhoton_Pt_    );
  electronTree_->Branch("genPhoton_Eta"   ,  &genPhoton_Eta_   );
  electronTree_->Branch("genPhoton_Rap"   ,  &genPhoton_Rap_   );
  electronTree_->Branch("genPhoton_Phi"   ,  &genPhoton_Phi_   );
  electronTree_->Branch("genPhoton_En"    ,  &genPhoton_En_    );

  electronTree_->Branch("genPostFSR_Px"   ,  &genPostFSR_Px_    );
  electronTree_->Branch("genPostFSR_Py"   ,  &genPostFSR_Py_    );
  electronTree_->Branch("genPostFSR_Pz"   ,  &genPostFSR_Pz_    );
  electronTree_->Branch("genPostFSR_Pt"   ,  &genPostFSR_Pt_    );
  electronTree_->Branch("genPostFSR_Eta"  ,  &genPostFSR_Eta_   );
  electronTree_->Branch("genPostFSR_Rap"  ,  &genPostFSR_Rap_   );
  electronTree_->Branch("genPostFSR_Phi"  ,  &genPostFSR_Phi_   );
  electronTree_->Branch("genPostFSR_En"   ,  &genPostFSR_En_    );

  electronTree_->Branch("genPreFSR_Px"   ,  &genPreFSR_Px_    );
  electronTree_->Branch("genPreFSR_Py"   ,  &genPreFSR_Py_    );
  electronTree_->Branch("genPreFSR_Pz"   ,  &genPreFSR_Pz_    );
  electronTree_->Branch("genPreFSR_Pt"   ,  &genPreFSR_Pt_    );
  electronTree_->Branch("genPreFSR_Eta"  ,  &genPreFSR_Eta_   );
  electronTree_->Branch("genPreFSR_Rap"  ,  &genPreFSR_Rap_   );
  electronTree_->Branch("genPreFSR_Phi"  ,  &genPreFSR_Phi_   );
  electronTree_->Branch("genPreFSR_En"   ,  &genPreFSR_En_    );

  electronTree_->Branch("ptElec"    ,  &ptElec_    );
  electronTree_->Branch("etaElec"   ,  &etaElec_    );
  electronTree_->Branch("rapElec"   ,  &rapElec_    );
  electronTree_->Branch("phiElec"   ,  &phiElec_    );
  electronTree_->Branch("energyElec",  &energyElec_    );
  electronTree_->Branch("massElec"  ,  &massElec_    );
  electronTree_->Branch("chargeElec",  &chargeElec_    );
  electronTree_->Branch("enSC"    ,  &enSC_ );
  electronTree_->Branch("preEnSC" ,  &preEnSC_ );
  electronTree_->Branch("rawEnSC" ,  &rawEnSC_ );
  electronTree_->Branch("etSC"  ,  &etSC_ );
  electronTree_->Branch("etaSC" ,  &etaSC_ );
  electronTree_->Branch("phiSC" ,  &phiSC_ );
  electronTree_->Branch("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta_);
  electronTree_->Branch("E1x5" ,  &E1x5_ );
  electronTree_->Branch("E2x5" ,  &E2x5_ );
  electronTree_->Branch("E5x5" ,  &E5x5_ );
  electronTree_->Branch("hOverE",  &hOverE_);
  electronTree_->Branch("etaScWidth" ,  &etaScWidth_ );
  electronTree_->Branch("phiScWidth" ,  &phiScWidth_ );
  electronTree_->Branch("r9",  &r9_);
  
  electronTree_->Branch("error_scale", &error_scale_);
  electronTree_->Branch("sigma", &sigma_);
  electronTree_->Branch("sigma1", &sigma1_);
  electronTree_->Branch("sigma2", &sigma2_);
  electronTree_->Branch("sigma3", &sigma3_);
  electronTree_->Branch("sigma4", &sigma4_);

  electronTree_->Branch("dEtaIn",  &dEtaIn_);
  electronTree_->Branch("dPhiIn",  &dPhiIn_);
  electronTree_->Branch("isoChargedHadrons"      , &isoChargedHadrons_);
  electronTree_->Branch("isoNeutralHadrons"      , &isoNeutralHadrons_);
  electronTree_->Branch("isoPhotons"             , &isoPhotons_);
  electronTree_->Branch("isoChargedFromPU"       , &isoChargedFromPU_);
  electronTree_->Branch("isoDeltaBeta"     , &isoDeltaBeta_);
  electronTree_->Branch("isoRho"     , &isoRho_);
  electronTree_->Branch("ooEmooP", &ooEmooP_);
  electronTree_->Branch("d0"     , &d0_);
  electronTree_->Branch("dz"     , &dz_);
  electronTree_->Branch("expectedMissingInnerHits", &expectedMissingInnerHits_);
  electronTree_->Branch("passConversionVeto", &passConversionVeto_);
  electronTree_->Branch("brem",  &brem_);

  electronTree_->Branch("passVetoId"  ,  &passVetoId_ );
  electronTree_->Branch("passLooseId"  ,  &passLooseId_ );
  electronTree_->Branch("passMediumId" ,  &passMediumId_ );
  electronTree_->Branch("passTightId"  ,  &passTightId_ );
  electronTree_->Branch("passHEEPId" , &passHEEPId_ );

  electronTree_->Branch("isPassMedium_NoPt", &isPassMedium_NoPt_);
  electronTree_->Branch("isPassMedium_NoScEta", &isPassMedium_NoScEta_);
  electronTree_->Branch("isPassMedium_NoDEta", &isPassMedium_NoDEta_);
  electronTree_->Branch("isPassMedium_NoDPhi", &isPassMedium_NoDPhi_);
  electronTree_->Branch("isPassMedium_NoSigmaEtaEta", &isPassMedium_NoSigmaEtaEta_);
  electronTree_->Branch("isPassMedium_NoHOverE", &isPassMedium_NoHOverE_);
  electronTree_->Branch("isPassMedium_NoDxy", &isPassMedium_NoDxy_);
  electronTree_->Branch("isPassMedium_NoDz", &isPassMedium_NoDz_);
  electronTree_->Branch("isPassMedium_NoEInvP", &isPassMedium_NoEInvP_);
  electronTree_->Branch("isPassMedium_NoPFIso", &isPassMedium_NoPFIso_);
  electronTree_->Branch("isPassMedium_NoConVeto", &isPassMedium_NoConVeto_);
  electronTree_->Branch("isPassMedium_NoMissHits", &isPassMedium_NoMissHits_);

  electronTree_->Branch("eleEcalDrivenSeed", &eleEcalDrivenSeed_);
  electronTree_->Branch("eleInBarrel", &eleInBarrel_);
  electronTree_->Branch("eleInEndcap", &eleInEndcap_);

  electronTree_->Branch("nMuons", &nMuons);
  electronTree_->Branch("isLooseMuon", &isLooseMuon_);
  electronTree_->Branch("isTightMuon", &isTightMuon_);
  electronTree_->Branch("isHighPtMuon", &isHighPtMuon_);

  electronTree_->Branch("ptMuon", &ptMuon_);
  electronTree_->Branch("etaMuon", &etaMuon_);
  electronTree_->Branch("phiMuon", &phiMuon_);
  electronTree_->Branch("energyMuon", &energyMuon_);
  electronTree_->Branch("chargeMuon", &chargeMuon_);

  electronTree_->Branch("isoChargedHadronPfR04Muon", &isoChargedHadronPfR04Muon_);
  electronTree_->Branch("isoNeutralHadronPfR04Muon", &isoNeutralHadronPfR04Muon_);
  electronTree_->Branch("isoGammaPfR04Muon", &isoGammaPfR04Muon_);
  electronTree_->Branch("isoChargedFromPUMuon", &isoChargedFromPUMuon_);
  electronTree_->Branch("isoPFMuon", &isoPFMuon_);
  electronTree_->Branch("isoTrkMuon", &isoTrkMuon_);

  electronTree_->Branch("metPt", &metPt_);
  electronTree_->Branch("metPhi", &metPhi_);
  electronTree_->Branch("metSumEt", &metSumEt_);

  electronTree_->Branch("nPhotons",  &nPhotons , "nPhotons/I");
  electronTree_->Branch("ptPhoton"  ,  &ptPhoton_    );
  electronTree_->Branch("etaPhoton" ,  &etaPhoton_ );
  electronTree_->Branch("phiPhoton" ,  &phiPhoton_ );
}

SimpleElectronNtupler::~SimpleElectronNtupler()
{
}

// member functions

// ------------ method called for each event  ------------
  void
SimpleElectronNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace std;
  using namespace edm;
  using namespace reco;

  RunNo = iEvent.id().run();
  EvtNo = iEvent.id().event();
  Lumi  = iEvent.luminosityBlock();
  Bunch = iEvent.bunchCrossing();

  //cout<<"1"<<endl;

  if(misMC){

    edm::Handle<GenEventInfoProduct> genEvtInfo;
    //edm::Handle<edm::View<GenEventInfoProduct> > genEvtInfo;   
    //iEvent.getByToken("generator", genEvtInfo);
    iEvent.getByToken(genToken_,genEvtInfo);

    if (genEvtInfo.isValid()) {
      theWeight = genEvtInfo->weight();
      //cout << theWeight << endl;

    }
  }

  // Get Pileup info
  if(misMC){
    Handle<edm::View<PileupSummaryInfo> > pileupHandle;
    iEvent.getByToken(pileupToken_, pileupHandle);
    for( auto & puInfoElement : *pileupHandle){
      if( puInfoElement.getBunchCrossing() == 0 ){
	nPU     = puInfoElement.getPU_NumInteractions();
	nPUTrue = puInfoElement.getTrueNumInteractions();
      }
    }
  }

  //cout<<"2"<<endl;

  // Get Triggers
  Handle<edm::TriggerResults> triggerHandle;
  iEvent.getByToken(triggerToken_, triggerHandle);

  Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescale_, triggerPrescales);

  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerHandle);

  // Filter Labels
  std::string SEFilter("hltEle23WPLooseGsfTrackIsoFilter");

  std::string photon22Filter("hltEG22HEFilter");
  std::string photon30Filter("hltEG30HEFilter");
  std::string photon36Filter("hltEG36HEFilter");
  std::string photon50Filter("hltEG50HEFilter");
  std::string photon75Filter("hltEG75HEFilter");
  std::string photon90Filter("hltEG90HEFilter");
  std::string photon120Filter("hltEG120HEFilter");
  std::string photon175Filter("hltEG175HEFilter");

  /*
     bool trigResult = false;
     for (unsigned int i=0; i<triggerHandle->size(); i++)
     {
     std::string trigName = triggerNames.triggerName(i);
     trigResult = triggerHandle->accept(i);
  //cout<<"Name of Trigger = "<<trigName<<"   Trigger Result = "<<trigResult<<"   Trigger Number = "<<i<<endl;
  //if(i==42) cout<<"trig result Ele23: "<<trigResult<<endl;
  }
  */

  singlePhoton = 0;
  prescalePhoton = 0;

  photon22 = false; photon30 = false; photon36 = false; photon50 = false; photon75 = false; photon90 = false; photon120 = false; photon175 = false;

  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    obj.unpackPathNames(triggerNames);
    //cout<<"Trigger: "<<obj.unpackPathNames(triggerNames)<<endl;
    for (unsigned j = 0; j < obj.filterLabels().size(); ++j){

      if((SEFilter.compare(obj.filterLabels()[j]))==0){
	//cout<<"pt: "<<obj.pt()<<endl;
	pt_Ele23.push_back(obj.pt());
	eta_Ele23.push_back(obj.eta());
	phi_Ele23.push_back(obj.phi());
      }

      if((photon22Filter.compare(obj.filterLabels()[j]))==0 || (photon30Filter.compare(obj.filterLabels()[j]))==0 || (photon36Filter.compare(obj.filterLabels()[j]))==0 || (photon50Filter.compare(obj.filterLabels()[j]))==0 || (photon75Filter.compare(obj.filterLabels()[j]))==0 || (photon90Filter.compare(obj.filterLabels()[j]))==0 || (photon120Filter.compare(obj.filterLabels()[j]))==0 || (photon175Filter.compare(obj.filterLabels()[j]))==0){

	if(obj.et() >= 22. && obj.et() < 30.){
	  if((photon22Filter.compare(obj.filterLabels()[j]))==0){
	    photon22 = true;
	    etSPhoHLT_.push_back(obj.et());
	    ptSPhoHLT_.push_back(obj.pt());
	    etaSPhoHLT_.push_back(obj.eta());
	    phiSPhoHLT_.push_back(obj.phi());
	  }
	}

	if(obj.et() >= 30. && obj.et() < 36.){
	  if((photon30Filter.compare(obj.filterLabels()[j]))==0){
	    photon30 = true;
	    etSPhoHLT_.push_back(obj.et()); 
	    ptSPhoHLT_.push_back(obj.pt());
	    etaSPhoHLT_.push_back(obj.eta());
	    phiSPhoHLT_.push_back(obj.phi());
	  }
	}

	if(obj.et() >=36. && obj.et() < 50.){
	  if((photon36Filter.compare(obj.filterLabels()[j]))==0){
	    photon36 = true;
	    etSPhoHLT_.push_back(obj.et());
	    ptSPhoHLT_.push_back(obj.pt());
	    etaSPhoHLT_.push_back(obj.eta());
	    phiSPhoHLT_.push_back(obj.phi());
	  }
	}

	if(obj.et() >= 50. && obj.et() < 75.){
	  if((photon50Filter.compare(obj.filterLabels()[j]))==0){
	    photon50 = true;
	    etSPhoHLT_.push_back(obj.et());
	    ptSPhoHLT_.push_back(obj.pt());
	    etaSPhoHLT_.push_back(obj.eta());
	    phiSPhoHLT_.push_back(obj.phi());
	  }
	}

	if(obj.et() >= 75. && obj.et() < 90.){
	  if((photon75Filter.compare(obj.filterLabels()[j]))==0){
	    photon75 = true;
	    etSPhoHLT_.push_back(obj.et());
	    ptSPhoHLT_.push_back(obj.pt());
	    etaSPhoHLT_.push_back(obj.eta());
	    phiSPhoHLT_.push_back(obj.phi());
	  }
	}

	if(obj.et() >= 90. && obj.et() < 120.){
	  if((photon90Filter.compare(obj.filterLabels()[j]))==0){
	    photon90 = true;
	    etSPhoHLT_.push_back(obj.et());
	    ptSPhoHLT_.push_back(obj.pt());
	    etaSPhoHLT_.push_back(obj.eta());
	    phiSPhoHLT_.push_back(obj.phi());
	  }
	}

	if(obj.et() >= 120. && obj.et() < 175.){
	  if((photon120Filter.compare(obj.filterLabels()[j]))==0){
	    photon120 = true;
	    etSPhoHLT_.push_back(obj.et());
	    ptSPhoHLT_.push_back(obj.pt());
	    etaSPhoHLT_.push_back(obj.eta());
	    phiSPhoHLT_.push_back(obj.phi());
	  }
	}

	if(obj.et() >= 175.){
	  if((photon175Filter.compare(obj.filterLabels()[j]))==0){
	    photon175 = true;
	    etSPhoHLT_.push_back(obj.et());
	    ptSPhoHLT_.push_back(obj.pt());
	    etaSPhoHLT_.push_back(obj.eta());
	    phiSPhoHLT_.push_back(obj.phi());
	  }
	}

      } // OR of triggers
    } // obj.filterLabels().size()
  } // triggerObjects

  for (unsigned int i=0; i<triggerHandle->size(); i++)
  {

    if(std::regex_match(triggerNames.triggerName(i),ele23_WPLoose)) Ele23_WPLoose = triggerHandle->accept(i);
    if(std::regex_match(triggerNames.triggerName(i),ele27_WP85)) Ele27_WP85       = triggerHandle->accept(i);
    if(std::regex_match(triggerNames.triggerName(i),isoMu20)) IsoMu20             = triggerHandle->accept(i);
    if(std::regex_match(triggerNames.triggerName(i),ele17_ele12)) Ele17_Ele12     = triggerHandle->accept(i);
    if(std::regex_match(triggerNames.triggerName(i),mu8_ele17)) Mu8_Ele17         = triggerHandle->accept(i);
    if(std::regex_match(triggerNames.triggerName(i),mu17_ele12)) Mu17_Ele12       = triggerHandle->accept(i);
    if(std::regex_match(triggerNames.triggerName(i),mu8_ele23)) Mu8_Ele23         = triggerHandle->accept(i);

    //if(i==42) cout<<"Ele23: "<<Ele23_WPLoose<<endl;

    if(photon22){
      if(std::regex_match(triggerNames.triggerName(i),Photon_22)) {
	singlePhoton = triggerHandle->accept(i);
	prescalePhoton = triggerPrescales->getPrescaleForIndex(i);
      }
    }

    if(photon30){
      if(std::regex_match(triggerNames.triggerName(i),Photon_30)) {
	singlePhoton = triggerHandle->accept(i);
	prescalePhoton = triggerPrescales->getPrescaleForIndex(i);
      }
    }

    if(photon36){
      if(std::regex_match(triggerNames.triggerName(i),Photon_36)) {
	singlePhoton = triggerHandle->accept(i);
	prescalePhoton = triggerPrescales->getPrescaleForIndex(i);
      }
    }

    if(photon50){
      if(std::regex_match(triggerNames.triggerName(i),Photon_50)) {
	singlePhoton = triggerHandle->accept(i);
	prescalePhoton = triggerPrescales->getPrescaleForIndex(i);
      }
    }

    if(photon75){
      if(std::regex_match(triggerNames.triggerName(i),Photon_75)) {
	singlePhoton = triggerHandle->accept(i);
	prescalePhoton = triggerPrescales->getPrescaleForIndex(i);
      }
    }

    if(photon90){
      if(std::regex_match(triggerNames.triggerName(i),Photon_90)) {
	singlePhoton = triggerHandle->accept(i);
	prescalePhoton = triggerPrescales->getPrescaleForIndex(i);
      }
    }

    if(photon120){
      if(std::regex_match(triggerNames.triggerName(i),Photon_120)) {
	singlePhoton = triggerHandle->accept(i);
	prescalePhoton = triggerPrescales->getPrescaleForIndex(i);
      }
    }

    if(photon175){
      if(std::regex_match(triggerNames.triggerName(i),Photon_175)) {
	singlePhoton = triggerHandle->accept(i);
	prescalePhoton = triggerPrescales->getPrescaleForIndex(i);
      }
    }

  }

  //cout<<"singlePhoton: "<<singlePhoton<<"   "<<"prescalePhoton: "<<prescalePhoton<<endl;
  //cout<<"3"<<endl;

  // Get rho value
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho = *rhoH;

  // Get the beam spot
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_,theBeamSpot);  

  // We use exactly the same handle for AOD and miniAOD formats since pat::Electron objects can be recast as reco::GsfElectron objects.
  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  bool isAOD = true;
  iEvent.getByToken(electronsToken_, electrons);
  if( !electrons.isValid() ){
    isAOD = false;
    iEvent.getByToken(electronsMiniAODToken_,electrons);
  }


  // Get the MC collection
  Handle<edm::View<reco::GenParticle> > genParticles;
  if( isAOD )
    iEvent.getByToken(genParticlesToken_,genParticles);
  else
    iEvent.getByToken(genParticlesMiniAODToken_,genParticles);

  if(misMC && misSIG){

    nGenElectrons = 0;
    nGenTaus = 0;
    tauFlag  = 0;

    for(size_t i = 0; i < genParticles->size(); ++i){
      const GenParticle &genlep = (*genParticles)[i];

      //if(abs(genlep.pdgId())==22) genlepMother_Id_.push_back(genlep.mother(0)->pdgId());

      //fromHProcessFinalState_.push_back(genlep.fromHardProcessFinalState());
      //fromHProcessDecayed_.push_back(genlep.fromHardProcessDecayed());

      nGenElectrons++;


      // gen Photon
      if(abs(genlep.pdgId())==22){

	genPhoton_Px_.push_back(genlep.px());
	genPhoton_Py_.push_back(genlep.py());
	genPhoton_Pz_.push_back(genlep.pz());
	genPhoton_Pt_.push_back(genlep.pt());
	genPhoton_Eta_.push_back(genlep.eta());
	genPhoton_Rap_.push_back(genlep.rapidity());
	genPhoton_Phi_.push_back(genlep.phi());
	genPhoton_En_.push_back(genlep.energy());
      }

      // gen Post FSR
      if(abs(genlep.pdgId())==11 && genlep.fromHardProcessFinalState()==1){

	genPostFSR_Px_.push_back(genlep.px());
	genPostFSR_Py_.push_back(genlep.py());
	genPostFSR_Pz_.push_back(genlep.pz());
	genPostFSR_Pt_.push_back(genlep.pt());
	genPostFSR_Eta_.push_back(genlep.eta());
	genPostFSR_Rap_.push_back(genlep.rapidity());
	genPostFSR_Phi_.push_back(genlep.phi());
	genPostFSR_En_.push_back(genlep.energy());
      }

      // gen info with HardProcess Flag

      if(abs(genlep.pdgId())==11 && genlep.isHardProcess()==1){

	genPreFSR_Px_.push_back(genlep.px());
	genPreFSR_Py_.push_back(genlep.py());
	genPreFSR_Pz_.push_back(genlep.pz());
	genPreFSR_Pt_.push_back(genlep.pt());
	genPreFSR_Eta_.push_back(genlep.eta());
	genPreFSR_Rap_.push_back(genlep.rapidity());
	genPreFSR_Phi_.push_back(genlep.phi());
	genPreFSR_En_.push_back(genlep.energy());
      }
      //  cout << "******************************" << endl;

      // Separating the taus coming from Z decay
      if(abs(genlep.pdgId())==15 && genlep.fromHardProcessDecayed()==1) nGenTaus++;
    }

    // Select the events containing 2 taus from hard-process
    if(nGenTaus == 2) tauFlag = 1;

  } // isMC && isSIG

  // Get PV
  edm::Handle<reco::VertexCollection> vertices;
  if( isAOD )
    iEvent.getByToken(vtxToken_, vertices);
  else
    iEvent.getByToken(vtxMiniAODToken_, vertices);

  if (vertices->empty()) return; // skip the event if no PV found
  nPV = vertices->size();
  int firstGoodVertexIdx = 0;

  VertexCollection::const_iterator firstGoodVertex = vertices->end();
  for (VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    firstGoodVertex = vtx;
    break;
  }

  if (firstGoodVertex==vertices->end())
    return; // skip event if there are no good PVs

  // Seems always zero. Not stored in miniAOD...?
  pvNTracks = firstGoodVertex->nTracks();


  // Get the conversions collection
  edm::Handle<reco::ConversionCollection> conversions;
  if(isAOD)
    iEvent.getByToken(conversionsToken_, conversions);
  else
    iEvent.getByToken(conversionsMiniAODToken_, conversions);

  // Get the electron ID data from the event stream.
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  edm::Handle<edm::ValueMap<bool> > heep_id_decisions;


  iEvent.getByToken(eleVetoIdMapToken_ ,veto_id_decisions);
  iEvent.getByToken(eleLooseIdMapToken_ ,loose_id_decisions);
  iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(eleTightIdMapToken_ ,tight_id_decisions);
  iEvent.getByToken(eleHEEPIdMapToken_ ,heep_id_decisions);

  // Full cut flow info for one of the working points:
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > medium_id_cutflow;
  iEvent.getByToken(eleMediumIdFullInfoMapToken_,medium_id_cutflow);

  nElectrons = 0;

  EnergyScaleCorrection_class eScaler("EgammaAnalysis/ElectronTools/data/ScalesSmearings/76X_16DecRereco_2015_Etunc");
  //EnergyScaleCorrection_class eScaler("EgammaAnalysis/ElectronTools/data/76X_16DecRereco_2015");
  float err_scale(0.), sig(0.), sig1(0.), sig2(0.), sig3(0.), sig4(0.);

  cout<<"Event No.: "<<EvtNo<<endl;
  
  // Loop over electrons
  for (size_t i = 0; i < electrons->size(); ++i){
    const auto el = electrons->ptrAt(i);

    // Kinematics

    if(el->pt() > 9. && el->eta() <= 2.5){

      nElectrons++;
      ptElec_.push_back( el->pt() );
      etaElec_.push_back( el->eta() );
      rapElec_.push_back( el->rapidity() );
      phiElec_.push_back( el->phi() );
      energyElec_.push_back( el->energy() );
      massElec_.push_back( el->mass() );
      chargeElec_.push_back( el->charge() );

      err_scale = eScaler.ScaleCorrectionUncertainty(RunNo, el->isEB(), el->r9(), el->superCluster()->eta(), el->et());
      
      cout<<"electron = "<<i<<" in barrel = "<<el->isEB()<<" r9 = "<<el->r9()<<" eta = "<<el->eta()<<" eta SC = "<<el->superCluster()->eta()<<" et = "<<el->et()<<" pt = "<<el->pt()<<" energy = "<<el->energy()<<endl;
      
      sig  = eScaler.getSmearingSigma(RunNo, el->isEB(), el->r9(), el->superCluster()->eta(), el->et(), 0, 0);
      sig1 = eScaler.getSmearingSigma(RunNo, el->isEB(), el->r9(), el->superCluster()->eta(), el->et(), 1, 0);
      sig2 = eScaler.getSmearingSigma(RunNo, el->isEB(), el->r9(), el->superCluster()->eta(), el->et(), -1, 0);
      sig3 = eScaler.getSmearingSigma(RunNo, el->isEB(), el->r9(), el->superCluster()->eta(), el->et(), 0, 1);
      sig4 = eScaler.getSmearingSigma(RunNo, el->isEB(), el->r9(), el->superCluster()->eta(), el->et(), 0, -1);
      
      error_scale_.push_back(err_scale);
      sigma_.push_back(sig);
      sigma1_.push_back(sig1);
      sigma2_.push_back(sig2);
      sigma3_.push_back(sig3);
      sigma4_.push_back(sig4);

      //cout<<"electron = "<<i<<endl;
      cout<<"error_scale = "<<err_scale<<endl;
      cout<<"sigma1 = "<<sig<<"   sigma2 = "<<sig1<<"   sigma3 = "<<sig2<<"   sigma4 = "<<sig3<<"   sigma5 = "<<sig4<<endl;
      cout<<""<<endl;

      double R = sqrt(el->superCluster()->x()*el->superCluster()->x() + el->superCluster()->y()*el->superCluster()->y() +el->superCluster()->z()*el->superCluster()->z());
      double Rt = sqrt(el->superCluster()->x()*el->superCluster()->x() + el->superCluster()->y()*el->superCluster()->y());

      enSC_.push_back(el->superCluster()->energy());
      preEnSC_.push_back(el->superCluster()->preshowerEnergy());
      rawEnSC_.push_back(el->superCluster()->rawEnergy());
      etSC_.push_back( (el->superCluster()->energy())*(Rt/R) );
      etaSC_.push_back( el->superCluster()->eta() );
      phiSC_.push_back( el->superCluster()->phi() );

      // ECAL
      full5x5_sigmaIetaIeta_.push_back( el->full5x5_sigmaIetaIeta() );
      E1x5_.push_back(el->e1x5());
      E2x5_.push_back(el->e2x5Max());
      E5x5_.push_back(el->e5x5());
      hOverE_.push_back( el->hcalOverEcal() );
      etaScWidth_.push_back(el->superCluster()->etaWidth());
      phiScWidth_.push_back(el->superCluster()->phiWidth());
      r9_.push_back(el->r9());

      // ECAL + Track
      dEtaIn_.push_back( el->deltaEtaSuperClusterTrackAtVtx() );
      dPhiIn_.push_back( el->deltaPhiSuperClusterTrackAtVtx() );
      // |1/E-1/p| = |1/E - EoverPinner/E| is computed below. The if protects against ecalEnergy == inf or zero
      // (always the case for miniAOD for electrons <5 GeV)
      if( el->ecalEnergy() == 0 ){
	//printf("Electron energy is zero!\n");
	ooEmooP_.push_back( 1e30 );
      }else if( !std::isfinite(el->ecalEnergy())){
	//printf("Electron energy is not finite!\n");
	ooEmooP_.push_back( 1e30 );
      }else{
	ooEmooP_.push_back( fabs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy() ) );
      }

      // Isolation
      GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();

      // Compute isolation with delta beta correction for PU
      isoChargedHadrons_.push_back( pfIso.sumChargedHadronPt );
      isoNeutralHadrons_.push_back( pfIso.sumNeutralHadronEt );
      isoPhotons_.push_back( pfIso.sumPhotonEt );
      isoChargedFromPU_.push_back( pfIso.sumPUPt );


      float abseta = fabs(el->superCluster()->eta());

      // The effective areas constants file in the local release or default CMSSW, whichever is found
      edm::FileInPath eaConstantsFile("RecoEgamma/ElectronIdentification/data/PHYS14/effAreaElectrons_cone03_pfNeuHadronsAndPhotons.txt");
      EffectiveAreas effectiveAreas(eaConstantsFile.fullPath());
      float eA = effectiveAreas.getEffectiveArea(abseta);

      isoDeltaBeta_.push_back((pfIso.sumChargedHadronPt + max<float>( 0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt))/(el->pt()));
      isoRho_.push_back((pfIso.sumChargedHadronPt + max<float>( 0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * eA))/(el->pt()));

      // Track - Impact Parameter, Conversion rejection, Converted
      reco::GsfTrackRef theTrack = el->gsfTrack();
      d0_.push_back( (-1) * theTrack->dxy(firstGoodVertex->position() ) );
      dz_.push_back( theTrack->dz( firstGoodVertex->position() ) );
      expectedMissingInnerHits_.push_back(el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) );
      bool passConvVeto = !ConversionTools::hasMatchedConversion(*el, conversions, theBeamSpot->position());
      passConversionVeto_.push_back( (int) passConvVeto );
      brem_.push_back(el->fbrem());


      // Electron ID
      bool isPassVeto  = (*veto_id_decisions)[el];

      bool isPassLoose  = (*loose_id_decisions)[el];
      bool isPassMedium = (*medium_id_decisions)[el];
      bool isPassTight  = (*tight_id_decisions)[el];
      bool isPassHEEP = (*heep_id_decisions)[el];

      passVetoId_.push_back  ( (int)isPassVeto  );
      passLooseId_.push_back ( (int)isPassLoose );
      passMediumId_.push_back( (int)isPassMedium);
      passTightId_.push_back ( (int)isPassTight );
      passHEEPId_.push_back ( (int)isPassHEEP );

      if( verboseIdFlag_ ) {
	vid::CutFlowResult fullCutFlowData = (*medium_id_cutflow)[el];

	// Full printout
	printf("\nDEBUG CutFlow, full info for cand with pt=%f:\n", el->pt());
	printCutFlowResult(fullCutFlowData);

	// Example of how to find the ID decision with one cut removed, this could be needed for N-1 studies.
	const int cutIndexToMask = 2; 
	// Here we masked the cut by cut index, but you can also do it by cut name string.
	vid::CutFlowResult maskedCutFlowData = fullCutFlowData.getCutFlowResultMasking(cutIndexToMask);
	printf("DEBUG CutFlow, the result with cut %s masked out\n", maskedCutFlowData.getNameAtIndex(cutIndexToMask).c_str());
	printCutFlowResult(maskedCutFlowData);
      }

      vid::CutFlowResult mediumID_Pt          = (*medium_id_cutflow)[el].getCutFlowResultMasking("MinPtCut_0");
      vid::CutFlowResult mediumID_ScEta       = (*medium_id_cutflow)[el].getCutFlowResultMasking("GsfEleSCEtaMultiRangeCut_0");
      vid::CutFlowResult mediumID_DEta        = (*medium_id_cutflow)[el].getCutFlowResultMasking("GsfEleDEtaInCut_0");
      vid::CutFlowResult mediumID_DPhi        = (*medium_id_cutflow)[el].getCutFlowResultMasking("GsfEleDPhiInCut_0");
      vid::CutFlowResult mediumID_SigmaEtaEta = (*medium_id_cutflow)[el].getCutFlowResultMasking("GsfEleFull5x5SigmaIEtaIEtaCut_0");
      vid::CutFlowResult mediumID_HOverE      = (*medium_id_cutflow)[el].getCutFlowResultMasking("GsfEleHadronicOverEMCut_0");
      vid::CutFlowResult mediumID_Dxy         = (*medium_id_cutflow)[el].getCutFlowResultMasking("GsfEleDxyCut_0");
      vid::CutFlowResult mediumID_Dz          = (*medium_id_cutflow)[el].getCutFlowResultMasking("GsfEleDzCut_0");
      vid::CutFlowResult mediumID_EInvP       = (*medium_id_cutflow)[el].getCutFlowResultMasking("GsfEleEInverseMinusPInverseCut_0");
      vid::CutFlowResult mediumID_PFIso       = (*medium_id_cutflow)[el].getCutFlowResultMasking("GsfEleEffAreaPFIsoCut_0");
      vid::CutFlowResult mediumID_ConVeto     = (*medium_id_cutflow)[el].getCutFlowResultMasking("GsfEleConversionVetoCut_0");
      vid::CutFlowResult mediumID_MissHits    = (*medium_id_cutflow)[el].getCutFlowResultMasking("GsfEleMissingHitsCut_0");

      isPassMedium_NoPt_.push_back(mediumID_Pt.cutFlowPassed());
      isPassMedium_NoScEta_.push_back(mediumID_ScEta.cutFlowPassed());
      isPassMedium_NoDEta_.push_back(mediumID_DEta.cutFlowPassed());
      //cout<<"isPassMedium_NoDEta: "<<mediumID_DEta.cutFlowPassed()<<endl;

      isPassMedium_NoDPhi_.push_back(mediumID_DPhi.cutFlowPassed());
      isPassMedium_NoSigmaEtaEta_.push_back(mediumID_SigmaEtaEta.cutFlowPassed());
      isPassMedium_NoHOverE_.push_back(mediumID_HOverE.cutFlowPassed());
      isPassMedium_NoDxy_.push_back(mediumID_Dxy.cutFlowPassed());
      isPassMedium_NoDz_.push_back(mediumID_Dz.cutFlowPassed());
      isPassMedium_NoEInvP_.push_back(mediumID_EInvP.cutFlowPassed());
      isPassMedium_NoPFIso_.push_back(mediumID_PFIso.cutFlowPassed());
      isPassMedium_NoConVeto_.push_back(mediumID_ConVeto.cutFlowPassed());
      isPassMedium_NoMissHits_.push_back(mediumID_MissHits.cutFlowPassed());

      eleInBarrel_.push_back(el->isEB());
      eleInEndcap_.push_back(el->isEE());

      // ECAL driven
      eleEcalDrivenSeed_.push_back(el->ecalDrivenSeed());

    }
  }


  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(muonsMiniAODToken_,muons);

  nMuons = 0;

  for(unsigned i=0; i < muons->size();++i ) {
    const auto mu = muons->ptrAt(i);

    if(mu->pt() > 10. && mu->eta() < 2.4){
      nMuons++;

      // Kinematics
      ptMuon_.push_back(mu->pt());
      etaMuon_.push_back(mu->eta());
      phiMuon_.push_back(mu->phi());
      energyMuon_.push_back(mu->energy());
      chargeMuon_.push_back(mu->charge());

      // isLoose
      bool isLooseId = muons->at(i).isLooseMuon();
      isLooseMuon_.push_back(isLooseId);

      // isTight
      bool isTightId = muons->at(i).isTightMuon(vertices->at(0));
      isTightMuon_.push_back(isTightId);

      // is High Pt
      bool isHighPt_Id = muons->at(i).isHighPtMuon(vertices->at(0));
      isHighPtMuon_.push_back(isHighPt_Id);

      // pf isolation 
      isoChargedHadronPfR04Muon_.push_back(mu->pfIsolationR04().sumChargedHadronPt);
      isoNeutralHadronPfR04Muon_.push_back(mu->pfIsolationR04().sumNeutralHadronEt);
      isoGammaPfR04Muon_.push_back(mu->pfIsolationR04().sumPhotonEt);
      isoChargedFromPUMuon_.push_back(mu->pfIsolationR04().sumPUPt);

      isoPFMuon_.push_back((mu->pfIsolationR04().sumChargedHadronPt + max<float>(0.0, mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - 0.5 * (mu->pfIsolationR04().sumPUPt)))/(mu->pt()));

      // tracker isolation
      isoTrkMuon_.push_back((mu->isolationR03().sumPt)/mu->pt());

    }
  }

  // MET
  edm::Handle<pat::METCollection> metHandle;
  iEvent.getByToken(metToken_, metHandle);

  const pat::MET &met = metHandle->front();
  metPt_.push_back(met.pt());
  metPhi_.push_back(met.phi());
  metSumEt_.push_back(met.sumEt());

  edm::Handle<edm::View<reco::Photon> > photons;
  iEvent.getByToken(photonsMiniAODToken_,photons);

  nPhotons = 0;

  // Loop over photons
  for (size_t i = 0; i < photons->size(); ++i){
    const auto pho = photons->ptrAt(i);

    // Kinematics
    if(pho->pt() < 15.){
      nPhotons++;

      // Kinematics
      ptPhoton_  .push_back(pho->pt());
      etaPhoton_ .push_back(pho->superCluster()->eta());
      phiPhoton_ .push_back(pho->superCluster()->phi());
    }
  }


  // Save this electron's info
  electronTree_->Fill();

  // Clear vectors
  pt_Ele23.clear();
  eta_Ele23.clear();
  phi_Ele23.clear();

  etSPhoHLT_.clear();
  ptSPhoHLT_.clear();
  etaSPhoHLT_.clear();
  phiSPhoHLT_.clear();

  if(misMC && misSIG){

    //genlep_Id_.clear();
    //genlepMother_Id_.clear();
    //fromHProcessFinalState_.clear();
    //fromHProcessDecayed_.clear();

    genPhoton_Pt_.clear();
    genPhoton_Px_.clear();
    genPhoton_Py_.clear();
    genPhoton_Pz_.clear();
    genPhoton_Eta_.clear();
    genPhoton_Rap_.clear();
    genPhoton_Phi_.clear();
    genPhoton_En_.clear();

    genPostFSR_Px_.clear();
    genPostFSR_Py_.clear();
    genPostFSR_Pz_.clear();
    genPostFSR_Pt_.clear();
    genPostFSR_Eta_.clear();
    genPostFSR_Rap_.clear();
    genPostFSR_Phi_.clear();
    genPostFSR_En_.clear();

    genPreFSR_Px_.clear();
    genPreFSR_Py_.clear();
    genPreFSR_Pz_.clear();
    genPreFSR_Pt_.clear();
    genPreFSR_Eta_.clear();
    genPreFSR_Rap_.clear();
    genPreFSR_Phi_.clear();
    genPreFSR_En_.clear();

  }

  ptElec_.clear();
  etaElec_.clear();
  rapElec_.clear();
  phiElec_.clear();
  energyElec_.clear();
  massElec_.clear();
  chargeElec_.clear();

  error_scale_.clear();
  sigma_.clear();
  sigma1_.clear();
  sigma2_.clear();
  sigma3_.clear();
  sigma4_.clear();

  enSC_.clear();
  preEnSC_.clear();
  rawEnSC_.clear();
  etSC_.clear();
  etaSC_.clear();
  phiSC_.clear();
  dEtaIn_.clear();
  dPhiIn_.clear();
  full5x5_sigmaIetaIeta_.clear();
  E1x5_.clear();
  E2x5_.clear();
  E5x5_.clear();
  hOverE_.clear();
  etaScWidth_.clear();
  phiScWidth_.clear();
  r9_.clear();

  isoChargedHadrons_.clear();
  isoNeutralHadrons_.clear();
  isoPhotons_.clear();
  isoChargedFromPU_.clear();
  isoDeltaBeta_.clear();
  isoRho_.clear();
  ooEmooP_.clear();
  d0_.clear();
  dz_.clear();
  expectedMissingInnerHits_.clear();
  passConversionVeto_.clear();
  brem_.clear();

  passVetoId_.clear();
  passLooseId_.clear();
  passMediumId_.clear();
  passTightId_.clear();
  passHEEPId_ .clear();

  isPassMedium_NoPt_.clear();
  isPassMedium_NoScEta_.clear();
  isPassMedium_NoDEta_.clear();
  isPassMedium_NoDPhi_.clear();
  isPassMedium_NoSigmaEtaEta_.clear();
  isPassMedium_NoHOverE_.clear();
  isPassMedium_NoDxy_.clear();
  isPassMedium_NoDz_.clear();
  isPassMedium_NoEInvP_.clear();
  isPassMedium_NoPFIso_.clear();
  isPassMedium_NoConVeto_.clear();
  isPassMedium_NoMissHits_.clear();

  eleEcalDrivenSeed_.clear();

  eleInBarrel_.clear();
  eleInEndcap_.clear();

  isLooseMuon_.clear();
  isTightMuon_.clear();
  isHighPtMuon_.clear();

  ptMuon_.clear();
  etaMuon_.clear();
  phiMuon_.clear();
  energyMuon_.clear();
  chargeMuon_.clear();

  isoChargedHadronPfR04Muon_.clear();
  isoNeutralHadronPfR04Muon_.clear();
  isoGammaPfR04Muon_.clear();
  isoChargedFromPUMuon_.clear();
  isoPFMuon_.clear();
  isoTrkMuon_.clear();

  metPt_.clear();
  metPhi_.clear();
  metSumEt_.clear();

  ptPhoton_.clear();
  etaPhoton_.clear();
  phiPhoton_.clear();

}
// ------------ method called once each job just before starting event loop  ------------
  void 
SimpleElectronNtupler::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
  void 
SimpleElectronNtupler::endJob() 
{

}

// ------------ method called when starting to processes a run  ------------
/*
   void 
   SimpleElectronNtupler::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void 
   SimpleElectronNtupler::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void 
   SimpleElectronNtupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void 
   SimpleElectronNtupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimpleElectronNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void SimpleElectronNtupler::printCutFlowResult(vid::CutFlowResult &cutflow){

  printf("    CutFlow name= %s    decision is %d\n", 
      cutflow.cutFlowName().c_str(),
      (int) cutflow.cutFlowPassed());
  int ncuts = cutflow.cutFlowSize();
  printf(" Index                               cut name              isMasked    value-cut-upon     pass?\n");
  for(int icut = 0; icut<ncuts; icut++){
    printf("  %2d      %50s    %d        %f          %d\n", icut,
	cutflow.getNameAtIndex(icut).c_str(),
	(int)cutflow.isCutMasked(icut),
	cutflow.getValueCutUpon(icut),
	(int)cutflow.getCutResultByIndex(icut));
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleElectronNtupler);
