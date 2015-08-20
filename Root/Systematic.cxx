#include "susynt-ss2lb/Systematic.h"

#include <cassert>

using hlfv::Systematic;

//-----------------------------------------
bool hlfv::isTriggerSyst(const Systematic::Value &s)
{
    return (s==Systematic::ETRIGREWUP   ||
            s==Systematic::ETRIGREWDOWN ||
            s==Systematic::MTRIGREWUP   ||
            s==Systematic::MTRIGREWDOWN ||
            s==Systematic::TTRIGSFUP    ||
            s==Systematic::TTRIGSFDOWN  );
}
//-----------------------------------------
Systematic::Value hlfv::ntsys2sys(const SusyNtSys &s)
{
    Systematic::Value r = Systematic::CENTRAL;
    switch(s) {
    case NtSys_NOM             :  r =  Systematic::CENTRAL      ; break;
    case NtSys_EES_Z_UP        :  r =  Systematic::EESZUP       ; break;
    case NtSys_EES_Z_DN        :  r =  Systematic::EESZDOWN     ; break;
    case NtSys_EES_MAT_UP      :  r =  Systematic::EESMATUP     ; break;
    case NtSys_EES_MAT_DN      :  r =  Systematic::EESMATDOWN   ; break;
    case NtSys_EES_PS_UP       :  r =  Systematic::EESPSUP      ; break;
    case NtSys_EES_PS_DN       :  r =  Systematic::EESPSDOWN    ; break; 
    case NtSys_EES_LOW_UP      :  r =  Systematic::EESLOWUP     ; break;
    case NtSys_EES_LOW_DN      :  r =  Systematic::EESLOWDOWN   ; break;
    case NtSys_EER_UP          :  r =  Systematic::EERUP        ; break;
    case NtSys_EER_DN          :  r =  Systematic::EERDOWN      ; break; 
    case NtSys_MS_UP           :  r =  Systematic::MESUP        ; break;
    case NtSys_MS_DN           :  r =  Systematic::MESDOWN      ; break;
    case NtSys_ID_UP           :  r =  Systematic::MIDUP        ; break;
    case NtSys_ID_DN           :  r =  Systematic::MIDDOWN      ; break;
    case NtSys_JES_UP          :  r =  Systematic::JESUP        ; break;
    case NtSys_JES_DN          :  r =  Systematic::JESDOWN      ; break;
    case NtSys_JER             :  r =  Systematic::JER          ; break;
    case NtSys_SCALEST_UP      :  r =  Systematic::SCALESTUP    ; break;
    case NtSys_SCALEST_DN      :  r =  Systematic::SCALESTDOWN  ; break;
    case NtSys_RESOST          :  r =  Systematic::RESOST       ; break;
    case NtSys_TRIGSF_EL_UP    :  r =  Systematic::ETRIGREWUP   ; break;
    case NtSys_TRIGSF_EL_DN    :  r =  Systematic::ETRIGREWDOWN ; break;
    case NtSys_TRIGSF_MU_UP    :  r =  Systematic::MTRIGREWUP   ; break;
    case NtSys_TRIGSF_MU_DN    :  r =  Systematic::MTRIGREWDOWN ; break;
    case NtSys_TES_UP          :  r =  Systematic::TESUP        ; break;
    case NtSys_TES_DN          :  r =  Systematic::TESDOWN      ; break;
    case NtSys_JVF_UP          :  r =  Systematic::JVFUP        ; break;
    case NtSys_JVF_DN          :  r =  Systematic::JVFDOWN      ; break;        
    case NtSys_N               : assert(false)         ; break; // perhaps throw an exception instead
        // no default, so that the compiler will warn us of un-handled cases
    }
    return r;
}
//-----------------------------------------
SusyNtSys hlfv::sys2ntsys(const Systematic::Value &s)
{
    // Here we assume that NtSys_N <= N(Systematic); in general this
    // is true because we have other syst (eg. fakes) in addition to
    // the SusyNt ones.
    // Also note that this implementation relies completely on
    // ntsys2sys (where we do an explicit value check). This allows to
    // have, for example:
    // NtSys_EES_MAT_UP=4    and    Systematic::EESMATUP=14.
    int r=0;
    while(r<NtSys_N) { if(s==ntsys2sys(static_cast<SusyNtSys>(r))) break; else r++;}
    return static_cast<SusyNtSys>(r);
}
//-----------------------------------------
BTagSys hlfv::sys2ntbsys(const Systematic::Value &sys)
{
    BTagSys bsys(sys==Systematic::BJETUP      ? BTag_BJet_UP :
                 sys==Systematic::BJETDOWN    ? BTag_BJet_DN :
                 sys==Systematic::CJETUP      ? BTag_CJet_UP :
                 sys==Systematic::CJETDOWN    ? BTag_CJet_DN :
                 sys==Systematic::BMISTAGUP   ? BTag_LJet_UP :
                 sys==Systematic::BMISTAGDOWN ? BTag_LJet_DN :
                 BTag_NOM);
    return bsys;
}
//-----------------------------------------

