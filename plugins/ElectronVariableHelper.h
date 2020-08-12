#ifndef _ELECTRONVARIABLEHELPER_H
#define _ELECTRONVARIABLEHELPER_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
//#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "CondFormats/L1TObjects/interface/L1CaloGeometry.h"
#include "CondFormats/DataRecord/interface/L1CaloGeometryRecord.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include <DataFormats/PatCandidates/interface/Electron.h>

#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "EgammaAnalysis/TnPTreeProducer/plugins/WriteValueMap.h"
#include "isolations.h"

#include "TMath.h"

typedef edm::View<reco::Candidate> CandView;

template <class T>
class ElectronVariableHelper : public edm::EDProducer {
 public:
  explicit ElectronVariableHelper(const edm::ParameterSet & iConfig);
  virtual ~ElectronVariableHelper() ;

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;

private:
  edm::EDGetTokenT<std::vector<T> > probesToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<BXVector<l1t::EGamma> > l1EGToken_;
  edm::EDGetTokenT<CandView> pfCandToken_;
  edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandidatesToken_;
  edm::EDGetTokenT<double> rhoToken_;
};

template<class T>
ElectronVariableHelper<T>::ElectronVariableHelper(const edm::ParameterSet & iConfig) :
  probesToken_(consumes<std::vector<T> >(iConfig.getParameter<edm::InputTag>("probes"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"))),
  l1EGToken_(consumes<BXVector<l1t::EGamma> >(iConfig.getParameter<edm::InputTag>("l1EGColl"))),
  conversionsToken_(consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"))),
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  pfCandidatesToken_(consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))) {

  produces<edm::ValueMap<float> >("dz");
  produces<edm::ValueMap<float> >("dxy");
  produces<edm::ValueMap<float> >("sip");
  produces<edm::ValueMap<float> >("missinghits");
  produces<edm::ValueMap<float> >("gsfhits");
  produces<edm::ValueMap<float> >("l1e");
  produces<edm::ValueMap<float> >("l1et");
  produces<edm::ValueMap<float> >("l1eta");
  produces<edm::ValueMap<float> >("l1phi");
  produces<edm::ValueMap<float> >("pfPt");
  produces<edm::ValueMap<float> >("convVtxFitProb");
  produces<edm::ValueMap<float> >("kfhits");
  produces<edm::ValueMap<float> >("kfchi2");
  produces<edm::ValueMap<float> >("ioemiop");
  produces<edm::ValueMap<float> >("5x5circularity");
  produces<edm::ValueMap<float> >("pfLeptonIsolation");

  produces<edm::ValueMap<float> >("isPassConversionVeto");
  produces<edm::ValueMap<float> >("full5x5sigmaIetaIeta");
  produces<edm::ValueMap<float> >("hoeCutValueForLooseEB");
  produces<edm::ValueMap<float> >("hoeCutValueForLooseEE");

  produces<edm::ValueMap<float> >("hoeCutValue");
  produces<edm::ValueMap<float> >("emhadIsoCutValue");
  produces<edm::ValueMap<float> >("isPassHEEPV70");
  produces<edm::ValueMap<float> >("isPassHEEPV70For2018");

  if( iConfig.existsAs<edm::InputTag>("pfCandColl") ) {
    pfCandToken_ = consumes<CandView>(iConfig.getParameter<edm::InputTag>("pfCandColl"));
  }

}

template<class T>
ElectronVariableHelper<T>::~ElectronVariableHelper()
{}


template<class T>
void ElectronVariableHelper<T>::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // read input
  edm::Handle<std::vector<T> > probes;
  edm::Handle<reco::VertexCollection> vtxH;

  iEvent.getByToken(probesToken_, probes);
  iEvent.getByToken(vtxToken_, vtxH);
  const reco::VertexRef vtx(vtxH, 0);

  edm::Handle<BXVector<l1t::EGamma> > l1Cands;
  iEvent.getByToken(l1EGToken_, l1Cands);

  edm::Handle<reco::ConversionCollection> conversions;
  iEvent.getByToken(conversionsToken_, conversions);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpotToken_, beamSpotHandle);
  const reco::BeamSpot* beamSpot = &*(beamSpotHandle.product());

  edm::Handle<CandView> pfCands;
  if( !pfCandToken_.isUninitialized() ) iEvent.getByToken(pfCandToken_,pfCands);

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_, rhoHandle);
  double Rho = *rhoHandle;

  // prepare vector for output
  std::vector<float> dzVals;
  std::vector<float> dxyVals;
  std::vector<float> sipVals;
  std::vector<float> mhVals;

  std::vector<float> l1EVals;
  std::vector<float> l1EtVals;
  std::vector<float> l1EtaVals;
  std::vector<float> l1PhiVals;
  std::vector<float> pfPtVals;
  std::vector<float> convVtxFitProbVals;
  std::vector<float> kfhitsVals;
  std::vector<float> kfchi2Vals;
  std::vector<float> ioemiopVals;
  std::vector<float> ocVals;

  std::vector<float> gsfhVals;

  std::vector<float> isPassConversionVetoVals;
  std::vector<float> full5x5sigmaIetaIetaVals;
  std::vector<float> hoeCutValueForLooseEBVals;
  std::vector<float> hoeCutValueForLooseEEVals;

  std::vector<float> hoeCutValueVals;
  std::vector<float> emhadIsoCutValueVals;
  std::vector<float> isPassHEEPV70Vals;
  std::vector<float> isPassHEEPV70For2018Vals;

  typename std::vector<T>::const_iterator probe, endprobes = probes->end();

  for (probe = probes->begin(); probe != endprobes; ++probe) {

    //---Clone the pat::Electron
    pat::Electron l((pat::Electron)*probe);

    dzVals.push_back(probe->gsfTrack()->dz(vtx->position()));
    dxyVals.push_back(probe->gsfTrack()->dxy(vtx->position()));

    // SIP
    float IP      = fabs(l.dB(pat::Electron::PV3D));
    float IPError = l.edB(pat::Electron::PV3D);
    sipVals.push_back(IP/IPError);

    mhVals.push_back(float(probe->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)));
    gsfhVals.push_back(float(probe->gsfTrack()->hitPattern().trackerLayersWithMeasurement()));
    float l1e = 999999.;    
    float l1et = 999999.;
    float l1eta = 999999.;
    float l1phi = 999999.;
    float pfpt = 999999.;
    float dRmin = 0.3;

    for (std::vector<l1t::EGamma>::const_iterator l1Cand = l1Cands->begin(0); l1Cand != l1Cands->end(0); ++l1Cand) {

      float dR = deltaR(l1Cand->eta(), l1Cand->phi() , probe->superCluster()->eta(), probe->superCluster()->phi());
      if (dR < dRmin) {
        dRmin = dR;
        l1e = l1Cand->energy();
        l1et = l1Cand->et();
        l1eta = l1Cand->eta();
        l1phi = l1Cand->phi();
      }
    }
    if( pfCands.isValid() )
    for( size_t ipf = 0; ipf < pfCands->size(); ++ipf ) {
        auto pfcand = pfCands->ptrAt(ipf);
    if( abs( pfcand->pdgId() ) != 11 ) continue;
    float dR = deltaR(pfcand->eta(), pfcand->phi() , probe->eta(), probe->phi());
    if( dR < 0.0001 ) pfpt = pfcand->pt();
    }

    l1EVals.push_back(l1e);
    l1EtVals.push_back(l1et);
    l1EtaVals.push_back(l1eta);
    l1PhiVals.push_back(l1phi);
    pfPtVals.push_back(pfpt);

    // Conversion vertex fit
    reco::ConversionRef convRef = ConversionTools::matchedConversion(*probe, conversions, beamSpot->position());

    float convVtxFitProb = -1.;
    if(!convRef.isNull()) {
        const reco::Vertex &vtx = convRef.get()->conversionVertex();
        if (vtx.isValid()) {
            convVtxFitProb = TMath::Prob( vtx.chi2(),  vtx.ndof());
        }
    }

    convVtxFitProbVals.push_back(convVtxFitProb);

    // kf track related variables
    bool validKf=false;
    reco::TrackRef trackRef = probe->closestCtfTrackRef();
    validKf = trackRef.isAvailable();
    validKf &= trackRef.isNonnull();
    float kfchi2 = validKf ? trackRef->normalizedChi2() : 0 ; //ielectron->track()->normalizedChi2() : 0 ;
    float kfhits = validKf ? trackRef->hitPattern().trackerLayersWithMeasurement() : -1. ;

    kfchi2Vals.push_back(kfchi2);
    kfhitsVals.push_back(kfhits);

    // 5x5circularity
    float oc = probe->full5x5_e5x5() != 0. ? 1. - (probe->full5x5_e1x5() / probe->full5x5_e5x5()) : -1.;
    ocVals.push_back(oc);

    // 1/E - 1/p
    float ele_pin_mode  = probe->trackMomentumAtVtx().R();
    float ele_ecalE     = probe->ecalEnergy();
    float ele_IoEmIop   = -1;
    if(ele_ecalE != 0 || ele_pin_mode != 0) {
        ele_IoEmIop = 1.0 / ele_ecalE - (1.0 / ele_pin_mode);
    }

    ioemiopVals.push_back(ele_IoEmIop);

    // For Z'toNN
    float scEta = probe->superCluster()->eta();
    float scEnergy = probe->superCluster()->energy();
    float et = probe->et();
    float hoe = probe->hadronicOverEm();
    float ecaliso = probe->dr03EcalRecHitSumEt();
    float hcaliso = probe->dr03HcalDepth1TowerSumEt();

    isPassConversionVetoVals.push_back( probe->passConversionVeto() );
    full5x5sigmaIetaIetaVals.push_back( probe->full5x5_sigmaIetaIeta() );
    float hoeCutValueForLooseEB = 0.05 + 1.16/scEnergy + 0.0324*Rho/scEnergy;
    float hoeCutValueForLooseEE = 0.0441 + 2.54/scEnergy + 0.183*Rho/scEnergy;
    hoeCutValueForLooseEBVals.push_back(hoeCutValueForLooseEB);
    hoeCutValueForLooseEEVals.push_back(hoeCutValueForLooseEE);

    float hoeCutValue = (-0.4 + 0.4*fabs(scEta)) * Rho / scEnergy + 0.05;
    float emhadIsoCutValue = 2.5 + (0.15 + 0.07*fabs(scEta)) * Rho;
    if(et > 50.) emhadIsoCutValue = 2.5 + 0.03*(et - 50.) + (0.15 + 0.07*fabs(scEta)) * Rho;
    // Original cut values
    /*float hoeCutValue = 5/scEnergy + 0.05;
    float emhadIsoCutValue = 2.5 + 0.28*Rho;
    if(et > 50.) emhadIsoCutValue = 2.5 + 0.03*(et - 50.) + 0.28*Rho;*/
    hoeCutValueVals.push_back(hoeCutValue);
    emhadIsoCutValueVals.push_back(emhadIsoCutValue);

    /*int IDcutBitLoose = probe->userInt("cutBasedElectronID-Fall17-94X-V2-loose"); // Only valid for 2018 dataset
    float isPassLooseNoIso = 0;
    // 1023 - 2^7 = 895 
    if((IDcutBitLoose&895) == 895) isPassLooseNoIso = 1; 
    isPassLooseNoIsoVals.push_back(isPassLooseNoIso);*/

    float isPassHEEPV70 = probe->electronID("heepElectronID-HEEPV70");
    isPassHEEPV70Vals.push_back(isPassHEEPV70);

    int IDcutBit = probe->userInt("heepElectronID-HEEPV70");
    float isPassHEEPV70For2018 = 0;
    if(fabs(scEta) < 1.566){
      if(isPassHEEPV70 == 1) isPassHEEPV70For2018 = 1;
    }
    else{  // 4095 - 2^6 - 2^8 = 3775
      if(((IDcutBit&3775) == 3775) && (hoe < hoeCutValue) && (ecaliso + hcaliso < emhadIsoCutValue)) isPassHEEPV70For2018 = 1;
    }
    isPassHEEPV70For2018Vals.push_back(isPassHEEPV70For2018);

  }

  // PF lepton isolations
  edm::Handle<pat::PackedCandidateCollection> pfCandidates;
  iEvent.getByToken(pfCandidatesToken_, pfCandidates);
  auto pfLeptonIsolations = computePfLeptonIsolations(*probes, *pfCandidates);
  for(unsigned int i = 0; i < probes->size(); ++i) {
      pfLeptonIsolations[i] /= (*probes)[i].pt();
  }


  // convert into ValueMap and store
  writeValueMap(iEvent, probes, dzVals, "dz");
  writeValueMap(iEvent, probes, dxyVals, "dxy");
  writeValueMap(iEvent, probes, sipVals, "sip");
  writeValueMap(iEvent, probes, mhVals, "missinghits");
  writeValueMap(iEvent, probes, gsfhVals, "gsfhits");
  writeValueMap(iEvent, probes, l1EVals, "l1e");
  writeValueMap(iEvent, probes, l1EtVals, "l1et");
  writeValueMap(iEvent, probes, l1EtaVals, "l1eta");
  writeValueMap(iEvent, probes, l1PhiVals, "l1phi");
  writeValueMap(iEvent, probes, pfPtVals, "pfPt");
  writeValueMap(iEvent, probes, convVtxFitProbVals, "convVtxFitProb");
  writeValueMap(iEvent, probes, kfhitsVals, "kfhits");
  writeValueMap(iEvent, probes, kfchi2Vals, "kfchi2");
  writeValueMap(iEvent, probes, ioemiopVals, "ioemiop");
  writeValueMap(iEvent, probes, ocVals, "5x5circularity");
  writeValueMap(iEvent, probes, pfLeptonIsolations, "pfLeptonIsolation");

  writeValueMap(iEvent, probes, isPassConversionVetoVals, "isPassConversionVeto");
  writeValueMap(iEvent, probes, full5x5sigmaIetaIetaVals, "full5x5sigmaIetaIeta");
  writeValueMap(iEvent, probes, hoeCutValueForLooseEBVals, "hoeCutValueForLooseEB");
  writeValueMap(iEvent, probes, hoeCutValueForLooseEEVals, "hoeCutValueForLooseEE");

  writeValueMap(iEvent, probes, hoeCutValueVals, "hoeCutValue");
  writeValueMap(iEvent, probes, emhadIsoCutValueVals, "emhadIsoCutValue");
  writeValueMap(iEvent, probes, isPassHEEPV70Vals, "isPassHEEPV70");
  writeValueMap(iEvent, probes, isPassHEEPV70For2018Vals, "isPassHEEPV70For2018");
}

#endif
