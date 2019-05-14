#ifndef TAUTAGANDPROBEFILTER_H
#define TAUTAGANDPROBEFILTER_H

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>

#include <iostream>
#include <utility>
#include <vector>

using namespace edm;
using namespace std;
// using namespace reco;


class TauTagAndProbeFilter : public edm::EDFilter {

    public:
        TauTagAndProbeFilter(const edm::ParameterSet &);
        ~TauTagAndProbeFilter();

    private:
        bool filter(edm::Event &, edm::EventSetup const&);

        float ComputeMT(math::XYZTLorentzVector visP4, const pat::MET& met);
        
        EDGetTokenT<edm::View<reco::GenParticleRefVector> > _genparticlesTag;
        EDGetTokenT<pat::TauRefVector>   _tausTag;
        EDGetTokenT<pat::MuonRefVector>  _muonsTag;
        EDGetTokenT<pat::METCollection>  _metTag;
        bool _useMassCuts;
        EDGetTokenT<edm::View<reco::GsfElectron> >  _electronsTag;
        edm::EDGetTokenT<edm::ValueMap<bool> > _eleLooseIdMapTag;
        bool _electronVeto;
        EDGetTokenT<pat::JetRefVector>  _bjetsTag;
};

TauTagAndProbeFilter::TauTagAndProbeFilter(const edm::ParameterSet & iConfig) :

_genparticlesTag (consumes<edm::View<reco::GenParticleRefVector> > (iConfig.getParameter<InputTag>("Genparticles"))),
_tausTag  (consumes<pat::TauRefVector>  ( iConfig.getParameter<InputTag>("taus"))),
_muonsTag (consumes<pat::MuonRefVector> (iConfig.getParameter<InputTag>("muons"))),
_metTag   (consumes<pat::METCollection> (iConfig.getParameter<InputTag>("met"))),
_electronsTag (consumes<edm::View<reco::GsfElectron> > (iConfig.getParameter<edm::InputTag>("electrons"))),
_eleLooseIdMapTag  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
_bjetsTag  (consumes<pat::JetRefVector>  (iConfig.getParameter<InputTag>("bjets")))
{
    produces <pat::TauRefVector>  (); // probe
    produces <pat::MuonRefVector> (); // tag
    _useMassCuts = iConfig.getParameter<bool>("useMassCuts");
    _electronVeto = iConfig.getParameter<bool>("eleVeto");   
}

TauTagAndProbeFilter::~TauTagAndProbeFilter()
{}

bool TauTagAndProbeFilter::filter(edm::Event & iEvent, edm::EventSetup const& iSetup)
{
    std::unique_ptr<reco::GenParticleRefVector> genTau ( new reco::GenParticleRefVector );
    std::unique_ptr<pat::MuonRefVector> resultMuon ( new pat::MuonRefVector );
    std::unique_ptr<pat::TauRefVector>  resultTau  ( new pat::TauRefVector  );  

    // Veto events with loose electrons
    if(_electronVeto){
      Handle<edm::View<reco::GsfElectron> > electrons;
      iEvent.getByToken(_electronsTag, electrons);
      Handle<edm::ValueMap<bool> > loose_id_decisions;
      iEvent.getByToken(_eleLooseIdMapTag, loose_id_decisions);
      
      for(unsigned int i = 0; i< electrons->size(); ++i){
	
	const auto ele = electrons->ptrAt(i);
	int isLooseID = (*loose_id_decisions)[ele];
	if(isLooseID && ele->p4().Pt()>10 && fabs(ele->p4().Eta())<2.5)
	  return false;
	
      }
      
    }
    //----------------------   search for the gentaus probe in the event --------------------
  
    edm::Handle<reco::GenParticleRefVector> GenParticles; 
    iEvent.getByToken (_genparticlesTag,GenParticles);
    const reco::GenParticleRef gen;
    const reco::GenParticleRef genmu;
    const reco::GenParticleRef gentau;
    for(unsigned int igm=0;igm<GenParticles->size();igm++){
      const reco::GenParticleRef gen=(*GenParticles)[igm];
      if(abs(gen->pdgId())==13)
	const reco::GenParticleRef genmu=(*GenParticles)[igm];
    }
    
    unsigned int gtau_ix=0;
    for(unsigned int ig=0;ig<GenParticles->size();ig++){
      if(abs(gen->pdgId())==15)
	const reco::GenParticleRef gentau=(*GenParticles)[ig];
    }
    math::XYZTLorentzVector pSum = genmu->p4() + gentau->p4();
    for(unsigned int igt=0;igt<GenParticles->size();igt++){
      if(abs(gen->pdgId())==15){
      
	if (_useMassCuts && (pSum.mass() <= 40 || pSum.mass() >= 80)) continue; // visible mass in (40, 80)
	// float dphi = gentau.phi()-genmuon.phi();
	// const double PI  =3.141592653589793238463;
	// if(dphi>PI)
	// 	dphi=2*PI - dphi;
	// float dR=sqrt((gentau.eta()-genmuon.eta())*(gentau.eta()-genmuon.eta())+dphi*dphi);
	if (deltaR(*gentau, *genmu) < 0.5) continue;
	gtau_ix = igt;
      }
    }
    reco::GenParticleRef gtaus;
    if(abs(gen->pdgId())==15)
      gtaus=(*GenParticles)[gtau_ix];
    
    
    
    // ---------------------   search for the tag in the event --------------------
    Handle<pat::MuonRefVector> muonHandle;
    iEvent.getByToken (_muonsTag, muonHandle);

    const pat::MuonRef mu = (*muonHandle)[0] ;

    //---------------------   get the met for mt computation etc. -----------------
    Handle<pat::METCollection> metHandle;
    iEvent.getByToken (_metTag, metHandle);
    const pat::MET& met = (*metHandle)[0];

    float mt = ComputeMT (mu->p4(), met);

    if (mt >= 30 && _useMassCuts) return false; // reject W+jets


    // ------------------- get Taus -------------------------------
    Handle<pat::TauRefVector> tauHandle;
    iEvent.getByToken (_tausTag, tauHandle);
    if (tauHandle->size() < 1) return false;

    vector<pair<float, int>> tausIdxPtVec;
    for (uint itau = 0; itau < tauHandle->size(); ++itau)
    {
        const pat::TauRef tau = (*tauHandle)[itau] ;
        math::XYZTLorentzVector pSum = mu->p4() + tau->p4();
        if (_useMassCuts && (pSum.mass() <= 40 || pSum.mass() >= 80)) continue; // visible mass in (40, 80)
        if (deltaR(*tau, *mu) < 0.5) continue;

        // min iso
        float isoMVA = tau->tauID("byIsolationMVArun2v1DBoldDMwLTraw");
        tausIdxPtVec.push_back(make_pair(isoMVA, itau));
    }


    pat::TauRef tau;

    if (tausIdxPtVec.size() == 0) return false; //No tau found
    if (tausIdxPtVec.size() > 1) sort (tausIdxPtVec.begin(), tausIdxPtVec.end()); //Sort if multiple taus
    int tauIdx = tausIdxPtVec.back().second; // min iso --> max MVA score
    tau = (*tauHandle)[tauIdx];


    // ----------------- b-jets veto ---------------------
    Handle<pat::JetRefVector> bjetHandle;
    iEvent.getByToken (_bjetsTag, bjetHandle);

    for(unsigned int ijet = 0; ijet < bjetHandle->size(); ijet++){

      const pat::JetRef bjet = (*bjetHandle)[ijet];
      if( deltaR(*mu,*bjet)>0.5 && deltaR(*tau,*bjet)>0.5 ) return false;
      
    }
   
    resultTau->push_back (tau);
    resultMuon->push_back (mu);
    genTau->push_back(gtaus);
    iEvent.put(std::move(resultMuon));
    iEvent.put(std::move(resultTau));
    iEvent.put(std::move(genTau));
    return true;
}

float TauTagAndProbeFilter::ComputeMT (math::XYZTLorentzVector visP4, const pat::MET& met)
{
  math::XYZTLorentzVector METP4 (met.pt()*TMath::Cos(met.phi()), met.pt()*TMath::Sin(met.phi()), 0, met.pt());
  float scalSum = met.pt() + visP4.pt();

  math::XYZTLorentzVector vecSum (visP4);
  vecSum += METP4;
  float vecSumPt = vecSum.pt();
  return sqrt (scalSum*scalSum - vecSumPt*vecSumPt);
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(TauTagAndProbeFilter);

#endif
