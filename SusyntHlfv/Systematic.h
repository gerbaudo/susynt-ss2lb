// emacs -*- C++ -*-
#ifndef SUSYNTHLVF_SYSTEMATICS_H
#define SUSYNTHLVF_SYSTEMATICS_H

#include "SusyNtuple/SusyDefs.h"
#include <string>

namespace hlfv {

///   List of possible systematic uncertainties for the hlfv study
/**
   Essentially a copy of what we used for susy::wh, with a few additional
   helper functions to handle names and to convert between
   hlfv::Systematic <-> SusyNtSys

   Note to self: enclose enum in a struct to avoid collisions; see
   explanation on
   <A HREF="http://stackoverflow.com/questions/7090130/enum-in-a-namespace">stackoveflow</A>.

   davide.gerbaudo@gmail.com
   Jun 2014
*/

struct Systematic {
    enum Value {
        CENTRAL             ///< Central value
        ,STAT               ///< Statistical Uncertainty (no tree/weight needed for Limit setting)
        ,JESUP              ///< Positive Jet energy scale
        ,JESDOWN            ///< Negative Jet energy scale
        ,JER                ///< Jet energy resolution
        ,EESLOWUP           ///< Positive shift in electron energy scale (LOW)
        ,EESLOWDOWN         ///< Negative shift in electron energy scale (LOW)
        ,EESMATUP           ///< Positive shift in electron energy scale (Mat)
        ,EESMATDOWN         ///< Negative shift in electron energy scale (Mat)
        ,EESPSUP            ///< Positive shift in electron energy scale (PS)
        ,EESPSDOWN          ///< Negative shift in electron energy scale (PS)
        ,EESZUP             ///< Positive shift in electron energy scale (Zee)
        ,EESZDOWN           ///< Negative shift in electron energy scale (Zee)
        ,EERUP              ///< Positive shift in electron energy resolution
        ,EERDOWN            ///< Negative shift in electron energy resolution
        ,MESUP              ///< Positive shift in muon energy scale
        ,MESDOWN            ///< Negative shift in muon energy scale
        ,MIDUP              ///< Positive shift in muon ID resolution
        ,MIDDOWN            ///< Negative shift in muon ID resolution
        ,MMSUP              ///< Positive shift in muon MS resolution
        ,MMSDOWN            ///< Negative shift in muon MS resolution
        ,ESFUP              ///< Positive shift in electron efficiency
        ,ESFDOWN            ///< Negative shift in electron efficiency
        ,MEFFUP             ///< Positive shift in muon efficiency
        ,MEFFDOWN           ///< Negative shift in muon efficiency
        ,ETRIGREWUP         ///< Positive shift in electron trigger weights
        ,ETRIGREWDOWN       ///< Negative shift in electron trigger weights
        ,MTRIGREWUP         ///< Positive shift in muon trigger weights
        ,MTRIGREWDOWN       ///< Negative shift in muon trigger weights
        ,BJETUP             ///< Positive shift in btag scale factor
        ,BJETDOWN           ///< Negative shift in btag scale factor
        ,CJETUP             ///< Positive shift in ctag scale factor
        ,CJETDOWN           ///< Negative shift in ctag scale factor
        ,BMISTAGUP          ///< Positive shift in ltag scale factor
        ,BMISTAGDOWN        ///< Negative shift in ltag scale factor
        ,SCALESTUP          ///< Positive shift in MET soft term scale
        ,SCALESTDOWN        ///< Negative shift in MET soft term scale
        ,RESOST             ///< MET soft term resolution
        ,GEN                ///< Uncertainty due to generator (Alpgen versus Sherpa,...) and patronshower/hadronization (Herwig versus Pythia) (possibly includes PDF choice as well, i.e. CTEQ versus HERA,...)
        ,GENUP              ///< Positive shift due to generator parameter variations (ktfac, qfac, ptmin, Iqopt, scale); includes also shift in ISR/FSR
        ,GENDOWN            ///< Negative shift due to generator parameter variations (ktfac, qfac, ptmin, Iqopt, scale) ;includes also  shift in ISR/FSR
        ,PDFERRUP           ///< Positive shift due to PDF errorset
        ,PDFERRDOWN         ///< Negative shift due to PDF errorset
        ,BKGMETHODUP        ///< Positive shift due to background estimation method (chargeflip, jetveto, fakeweights,...)
        ,BKGMETHODDOWN      ///< Negative shift due to background estimation method (chargeflip, jetveto, fakeweights,...)
        ,XSUP               ///< Positive shift in theoretical cross-section uncertainty
        ,XSDOWN             ///< Negative shift in theoretical cross-section uncertainty
        ,LUMI               ///< Uncertainty due to luminosity (no tree/weight needed for Limit setting)
        ,PILEUPUP           ///< Positive shift for mu
        ,PILEUPDOWN         ///< Negative shift for mu
        ,JVFUP              ///< Positive shift in Jet Vertex Fraction
        ,JVFDOWN            ///< Negative shift in Jet Vertex Fraction
        ,TESUP              ///< Positive shift in TES variation
        ,TESDOWN            ///< Negative shift in TES variation
        ,TIDSFUP            ///< Positive shift in tauID scale factor
        ,TIDSFDOWN          ///< Negative shift in tauID scale factor
        ,TEVSFUP            ///< Positive shift in tauEVeto scale factor
        ,TEVSFDOWN          ///< Negative shift in tauEVeto scale factor
        ,TTRIGSFUP          ///< Positive shift in tau trigger scale factor
        ,TTRIGSFDOWN        ///< Negative shift in tau trigger scale factor
        ,TFAKESFUP          ///< Positive shift in fake tau scale factors
        ,TFAKESFDOWN        ///< Negative shift in fake tau scale factors
        ,FakeTauBGSyst      ///< Uncertainty on fake tau background estimation
    }; // Value
};

const std::string SystematicNames[] = {
    "CENTRAL"
    ,"STAT"
    ,"JESUP"
    ,"JESDOWN"
    ,"JER"
    ,"EESLOWUP"
    ,"EESLOWDOWN"
    ,"EESMATUP"
    ,"EESMATDOWN"
    ,"EESPSUP"
    ,"EESPSDOWN"
    ,"EESZUP"
    ,"EESZDOWN"
    ,"EERUP"
    ,"EERDOWN"
    ,"MESUP"
    ,"MESDOWN"
    ,"MIDUP"
    ,"MIDDOWN"
    ,"MMSUP"
    ,"MMSDOWN"
    ,"ESFUP"
    ,"ESFDOWN"
    ,"MEFFUP"
    ,"MEFFDOWN"
    ,"ETRIGREWUP"
    ,"ETRIGREWDOWN"
    ,"MTRIGREWUP"
    ,"MTRIGREWDOWN"
    ,"BJETUP"
    ,"BJETDOWN"
    ,"CJETUP"
    ,"CJETDOWN"
    ,"BMISTAGUP"
    ,"BMISTAGDOWN"
    ,"SCALESTUP"
    ,"SCALESTDOWN"
    ,"RESOST"
    ,"GEN"
    ,"GENUP"
    ,"GENDOWN"
    ,"PDFERRUP"
    ,"PDFERRDOWN"
    ,"BKGMETHODUP"
    ,"BKGMETHODDOWN"
    ,"XSUP"
    ,"XSDOWN"
    ,"LUMI"
    ,"PILEUPUP"
    ,"PILEUPDOWN"
    ,"JVFUP"
    ,"JVFDOWN"
    ,"TESUP"
    ,"TESDOWN"
    ,"TIDSFUP"
    ,"TIDSFDOWN"
    ,"TEVSFUP"
    ,"TEVSFDOWN"
    ,"TTRIGSFUP"
    ,"TTRIGSFDOWN"
    ,"TFAKESFUP"
    ,"TFAKESFDOWN"
    ,"FakeTauBGSyst"
};

/// whether s is a valid enum value
inline bool isValid(const Systematic::Value &s) { return s>=Systematic::CENTRAL && s<=Systematic::FakeTauBGSyst; }
/// whether s is a valid enum value
inline bool isValid(const SusyNtSys &s)  { return s>=NtSys_NOM  && s< NtSys_N; }
/// string representation
inline std::string syst2str(const Systematic::Value &s) { return isValid(s) ? SystematicNames[s] : "unknown"; }
/// string representation
inline std::string syst2str(const SusyNtSys &s)  { return isValid(s) ? SusyNtSystNames[s] : "unknown"; }
/// convert Systematic to SusyNtSys
SusyNtSys sys2ntsys(const Systematic::Value &s);
/// convert SusyNtSys to Systematic; assert(false) if invalid
Systematic::Value ntsys2sys(const SusyNtSys &s);
/// convert Systematic to SusyDef::BTagSys; return nominal if not a btag-related sys
BTagSys sys2ntbsys(const Systematic::Value &s);

} // hlfv
#endif
