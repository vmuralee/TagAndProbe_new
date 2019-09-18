/*! Common simple types for analysis purposes.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#pragma once

#include "EnumNameMap.h"

namespace analysis {

using LorentzVectorXYZ = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>;
using LorentzVectorM = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;
using LorentzVectorE = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>>;
using LorentzVector = LorentzVectorE;

enum class DiscriminatorWP { VVVLoose = 0, VVLoose = 1, VLoose = 2, Loose = 3, Medium = 4, Tight = 5,
                             VTight = 6, VVTight = 7, VVVTight = 8 };
ENUM_NAMES(DiscriminatorWP) = {
    { DiscriminatorWP::VVVLoose, "VVVLoose" }, { DiscriminatorWP::VVLoose, "VVLoose" },
    { DiscriminatorWP::VLoose, "VLoose" }, { DiscriminatorWP::Loose, "Loose" }, { DiscriminatorWP::Medium, "Medium" },
    { DiscriminatorWP::Tight, "Tight" }, { DiscriminatorWP::VTight, "VTight" }, { DiscriminatorWP::VVTight, "VVTight" },
    { DiscriminatorWP::VVVTight, "VVVTight" }
};
const EnumNameMap<DiscriminatorWP> __DiscriminatorWP_short_names("ShortWPNames", {
    { DiscriminatorWP::VVVLoose, "VVVL" }, { DiscriminatorWP::VVLoose, "VVL" }, { DiscriminatorWP::VLoose, "VL" },
    { DiscriminatorWP::Loose, "L" }, { DiscriminatorWP::Medium, "M" }, { DiscriminatorWP::Tight, "T" },
    { DiscriminatorWP::VTight, "VT" }, { DiscriminatorWP::VVTight, "VVT" }, { DiscriminatorWP::VVVTight, "VVVT" }
});

enum class Period { Run2011, Run2012, Run2015, Run2016, Run2017, Run2018 };
ENUM_NAMES(Period) = {
    { Period::Run2011, "Run2011" }, { Period::Run2012, "Run2012" }, { Period::Run2015, "Run2015" },
    { Period::Run2016, "Run2016" }, { Period::Run2017, "Run2017" }, { Period::Run2018, "Run2018" },
};

enum class GenLeptonMatch { Electron = 1, Muon = 2, TauElectron = 3,  TauMuon = 4, Tau = 5, NoMatch = 6 };
ENUM_NAMES(GenLeptonMatch) = {
    { GenLeptonMatch::Electron, "gen_electron" },
    { GenLeptonMatch::Muon, "gen_muon" },
    { GenLeptonMatch::TauElectron, "gen_electron_from_tau" },
    { GenLeptonMatch::TauMuon, "gen_muon_from_tau" },
    { GenLeptonMatch::Tau, "gen_tau" },
    { GenLeptonMatch::NoMatch, "no_gen_match" }
};

enum class GenQcdMatch { NoMatch = 0, Down = 1, Up = 2, Strange = 3, Charm = 4, Bottom = 5, Top = 6, Gluon = 21 };
ENUM_NAMES(GenQcdMatch) = {
    { GenQcdMatch::NoMatch, "no_gen_match" },
    { GenQcdMatch::Down, "gen_down_quark" },
    { GenQcdMatch::Up, "gen_up_quark" },
    { GenQcdMatch::Strange, "gen_strange_quark" },
    { GenQcdMatch::Charm, "gen_charm_quark" },
    { GenQcdMatch::Bottom, "gen_bottom_quark" },
    { GenQcdMatch::Top, "gen_top_quark" },
    { GenQcdMatch::Gluon, "gen_gluon" }
};

enum class TauType { e = 0, mu = 1, tau = 2, jet = 3 };
ENUM_NAMES(TauType) = {
    { TauType::e, "e" }, { TauType::mu, "mu" }, { TauType::tau, "tau" }, { TauType::jet, "jet" }
};

inline constexpr TauType GenMatchToTauType(GenLeptonMatch gen_match)
{
    if(gen_match == GenLeptonMatch::Electron || gen_match == GenLeptonMatch::TauElectron) return TauType::e;
    if(gen_match == GenLeptonMatch::Muon || gen_match == GenLeptonMatch::TauMuon) return TauType::mu;
    if(gen_match == GenLeptonMatch::Tau) return TauType::tau;
    return TauType::jet;
}

} // namespace analysis
