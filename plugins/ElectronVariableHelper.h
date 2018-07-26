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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"


typedef edm::View<reco::Candidate> CandView;

template <class T>
class ElectronVariableHelper : public edm::EDProducer {
 public:
  explicit ElectronVariableHelper(const edm::ParameterSet & iConfig);
  virtual ~ElectronVariableHelper() ;
  
  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup) override;
  virtual float getEffArea(float scEta);
 
private:
  edm::EDGetTokenT<std::vector<T> > probesToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<BXVector<l1t::EGamma> > l1EGTkn;
  edm::EDGetTokenT<CandView> pfCandToken_;
  edm::EDGetTokenT<double> rhoLabel_;
//  edm::EDGetToken electronsMiniAODToken_;
//  edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
};

template<class T>
ElectronVariableHelper<T>::ElectronVariableHelper(const edm::ParameterSet & iConfig) :
  probesToken_(consumes<std::vector<T> >(iConfig.getParameter<edm::InputTag>("probes"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"))),
  l1EGTkn(consumes<BXVector<l1t::EGamma> >(iConfig.getParameter<edm::InputTag>("l1EGColl"))),
  rhoLabel_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoLabel"))) {
//  mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))) {

  produces<edm::ValueMap<float> >("chi2");
  produces<edm::ValueMap<float> >("dz");
  produces<edm::ValueMap<float> >("dxy");
  produces<edm::ValueMap<float> >("missinghits");
  produces<edm::ValueMap<float> >("l1e");
  produces<edm::ValueMap<float> >("l1et");
  produces<edm::ValueMap<float> >("l1eta");
  produces<edm::ValueMap<float> >("l1phi");
  produces<edm::ValueMap<float> >("pfPt");
  produces<edm::ValueMap<float> >("passConversionVeto");
  produces<edm::ValueMap<float> >("SIP3D");
  produces<edm::ValueMap<float> >("miniIso");
  produces<edm::ValueMap<float> >("hoeLooseBarrel");
  produces<edm::ValueMap<float> >("hoeLooseEndcap");
//  produces<edm::ValueMap<float> >("mva");

  if( iConfig.existsAs<edm::InputTag>("pfCandColl") ) {
    pfCandToken_ = consumes<CandView>(iConfig.getParameter<edm::InputTag>("pfCandColl"));
  }

//  electronsMiniAODToken_ = mayConsume<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electronsMiniAOD"));
//  electronsMiniAODToken_ = mayConsume<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronsMiniAOD"));
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
  iEvent.getByToken(l1EGTkn, l1Cands);
  
  edm::Handle<CandView> pfCands;
  if( !pfCandToken_.isUninitialized() ) iEvent.getByToken(pfCandToken_,pfCands);

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoLabel_, rhoHandle);
  double rho = std::max(*(rhoHandle.product()), 0.0);
 
/*  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByToken(electronsMiniAODToken_,electrons);

  edm::Handle<edm::ValueMap<float> > mvaValues;
  iEvent.getByToken(mvaValuesMapToken_,mvaValues);*/
 
  // prepare vector for output
  std::vector<float> chi2Vals;
  std::vector<float> dzVals;
  std::vector<float> dxyVals;
  std::vector<float> mhVals;
  std::vector<float> l1EVals;
  std::vector<float> l1EtVals;
  std::vector<float> l1EtaVals;
  std::vector<float> l1PhiVals;
  std::vector<float> pfPtVals;
  std::vector<float> passConversionVetoVals;
  std::vector<float> SIP3DVals;
  std::vector<float> miniIsoVals;
  std::vector<float> hoeLooseBarrelVals;
  std::vector<float> hoeLooseEndcapVals;
//  std::vector<float> mvaVals;

  typename std::vector<T>::const_iterator probe, endprobes = probes->end();

//  int j = 0;
  for (probe = probes->begin(); probe != endprobes; ++probe) {
    
    chi2Vals.push_back(probe->gsfTrack()->normalizedChi2());
    dzVals.push_back(probe->gsfTrack()->dz(vtx->position()));
    dxyVals.push_back(probe->gsfTrack()->dxy(vtx->position()));
    mhVals.push_back(float(probe->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)));

//    const auto el = electrons->ptrAt(j);

    passConversionVetoVals.push_back( probe->passConversionVeto() );
    float ip3d = probe->dB(pat::Electron::PV3D);
    float ip3derr = probe->edB(pat::Electron::PV3D);
    float sip3d = std::abs(ip3d/ip3derr);
    SIP3DVals.push_back(sip3d);

    float scEta = probe->superCluster()->eta();
    float scEnergy = probe->superCluster()->energy();
    float ecalpt = probe->ecalDrivenMomentum().pt();

    float chg = probe->miniPFIsolation().chargedHadronIso();
    float neu = probe->miniPFIsolation().neutralHadronIso();
    float pho = probe->miniPFIsolation().photonIso();
    float ea = getEffArea(scEta);
    float x0 = 10.0, x1 = 50.0, x2 = 200.0;
    float cone = x0/std::min(std::max(ecalpt, x1), x2);
    ea *= std::pow(cone/0.3, 2);
//    float x3 = 0.0;
    float minirelIso = (chg + std::max(neu + pho - rho*ea, 0.0))/ecalpt;
    miniIsoVals.push_back(minirelIso);
    float hoe_loose_barrel = 0.05 + 1.12/scEnergy + 0.0368*rho/scEnergy;
    float hoe_loose_endcap = 0.0414 + 0.5/scEnergy + 0.201*rho/scEnergy;
    hoeLooseBarrelVals.push_back(hoe_loose_barrel);
    hoeLooseEndcapVals.push_back(hoe_loose_endcap);  

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

//    float mvaval = (*mvaValues)[el];
//    float mvaval = el->userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV1Values");    
//    float mvaval = el->userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values");
//    mvaVals.push_back(mvaval);
//    ++j; 
  }

  
  // convert into ValueMap and store
  std::unique_ptr<edm::ValueMap<float> > chi2ValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler chi2Filler(*chi2ValMap);
  chi2Filler.insert(probes, chi2Vals.begin(), chi2Vals.end());
  chi2Filler.fill();
  iEvent.put(std::move(chi2ValMap), "chi2");

  std::unique_ptr<edm::ValueMap<float> > dzValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler dzFiller(*dzValMap);
  dzFiller.insert(probes, dzVals.begin(), dzVals.end());
  dzFiller.fill();
  iEvent.put(std::move(dzValMap), "dz");

  std::unique_ptr<edm::ValueMap<float> > dxyValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler dxyFiller(*dxyValMap);
  dxyFiller.insert(probes, dxyVals.begin(), dxyVals.end());
  dxyFiller.fill();
  iEvent.put(std::move(dxyValMap), "dxy");

  std::unique_ptr<edm::ValueMap<float> > mhValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler mhFiller(*mhValMap);
  mhFiller.insert(probes, mhVals.begin(), mhVals.end());
  mhFiller.fill();
  iEvent.put(std::move(mhValMap), "missinghits");

  std::unique_ptr<edm::ValueMap<float> > l1EValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler l1EFill(*l1EValMap);
  l1EFill.insert(probes, l1EVals.begin(), l1EVals.end());
  l1EFill.fill();
  iEvent.put(std::move(l1EValMap), "l1e");

  std::unique_ptr<edm::ValueMap<float> > l1EtValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler l1EtFill(*l1EtValMap);
  l1EtFill.insert(probes, l1EtVals.begin(), l1EtVals.end());
  l1EtFill.fill();
  iEvent.put(std::move(l1EtValMap), "l1et");

  std::unique_ptr<edm::ValueMap<float> > l1EtaValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler l1EtaFill(*l1EtaValMap);
  l1EtaFill.insert(probes, l1EtaVals.begin(), l1EtaVals.end());
  l1EtaFill.fill();
  iEvent.put(std::move(l1EtaValMap), "l1eta");

  std::unique_ptr<edm::ValueMap<float> > l1PhiValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler l1PhiFill(*l1PhiValMap);
  l1PhiFill.insert(probes, l1PhiVals.begin(), l1PhiVals.end());
  l1PhiFill.fill();
  iEvent.put(std::move(l1PhiValMap), "l1phi");

  std::unique_ptr<edm::ValueMap<float> > pfPtValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler pfPtFill(*pfPtValMap);
  pfPtFill.insert(probes, pfPtVals.begin(), pfPtVals.end());
  pfPtFill.fill();
  iEvent.put(std::move(pfPtValMap), "pfPt");

  std::unique_ptr<edm::ValueMap<float> > passConversionVetoValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler passConversionVetoFiller(*passConversionVetoValMap);
  passConversionVetoFiller.insert(probes, passConversionVetoVals.begin(), passConversionVetoVals.end());
  passConversionVetoFiller.fill();
  iEvent.put(std::move(passConversionVetoValMap), "passConversionVeto");  

  std::unique_ptr<edm::ValueMap<float> > SIP3DValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler SIP3DFill(*SIP3DValMap);
  SIP3DFill.insert(probes, SIP3DVals.begin(), SIP3DVals.end());
  SIP3DFill.fill();
  iEvent.put(std::move(SIP3DValMap), "SIP3D");

  std::unique_ptr<edm::ValueMap<float> > miniIsoValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler miniIsoFill(*miniIsoValMap);
  miniIsoFill.insert(probes, miniIsoVals.begin(), miniIsoVals.end());
  miniIsoFill.fill();
  iEvent.put(std::move(miniIsoValMap), "miniIso");

  std::unique_ptr<edm::ValueMap<float> > hoeLooseBarrelValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler hoeLooseBarrelFill(*hoeLooseBarrelValMap);
  hoeLooseBarrelFill.insert(probes, hoeLooseBarrelVals.begin(), hoeLooseBarrelVals.end());
  hoeLooseBarrelFill.fill();
  iEvent.put(std::move(hoeLooseBarrelValMap), "hoeLooseBarrel");

  std::unique_ptr<edm::ValueMap<float> > hoeLooseEndcapValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler hoeLooseEndcapFill(*hoeLooseEndcapValMap);
  hoeLooseEndcapFill.insert(probes, hoeLooseEndcapVals.begin(), hoeLooseEndcapVals.end());
  hoeLooseEndcapFill.fill();
  iEvent.put(std::move(hoeLooseEndcapValMap), "hoeLooseEndcap");

/*  std::unique_ptr<edm::ValueMap<float> > mvaValMap(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler mvaFiller(*mvaValMap);
  mvaFiller.insert(probes, mvaVals.begin(), mvaVals.end());
  mvaFiller.fill();
  iEvent.put(std::move(mvaValMap), "mva");*/
}

template<class T>
float ElectronVariableHelper<T>::getEffArea(float scEta){
  float absEta = std::abs(scEta);
  if (0.0000 >= absEta && absEta < 1.0000) return 0.1566;
  if (1.0000 >= absEta && absEta < 1.4790) return 0.1626;
  if (1.4790 >= absEta && absEta < 2.0000) return 0.1073;
  if (2.0000 >= absEta && absEta < 2.2000) return 0.0854;
  if (2.2000 >= absEta && absEta < 2.3000) return 0.1051;
  if (2.3000 >= absEta && absEta < 2.4000) return 0.1204;
  if (2.4000 >= absEta && absEta < 5.0000) return 0.1524;
  return 0;
}

#endif
