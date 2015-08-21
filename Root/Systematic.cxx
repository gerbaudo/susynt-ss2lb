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
Systematic::Value hlfv::ntsys2sys(const Susy::NtSys::SusyNtSys &s)
{
    Systematic::Value r = Systematic::CENTRAL;
    switch(s) {
    case Susy::NtSys::NOM             :  r =  Systematic::CENTRAL      ; break;
    case Susy::NtSys::TRIGSF_EL_UP    :  r =  Systematic::ETRIGREWUP   ; break;
    case Susy::NtSys::TRIGSF_EL_DN    :  r =  Systematic::ETRIGREWDOWN ; break;
    case Susy::NtSys::TRIGSF_MU_UP    :  r =  Systematic::MTRIGREWUP   ; break;
    case Susy::NtSys::TRIGSF_MU_DN    :  r =  Systematic::MTRIGREWDOWN ; break;
    case Susy::NtSys::SYS_UNKNOWN     :  assert(false); break; // perhaps throw an exception instead
    // no default, so that the compiler will warn us of un-handled cases
    default: assert(false); // temporarily (2015-08-21) suppress warnings
    }
    return r;
}
//-----------------------------------------
Susy::NtSys::SusyNtSys hlfv::sys2ntsys(const Systematic::Value &s)
{
    // Here we assume that NtSys_N <= N(Systematic); in general this
    // is true because we have other syst (eg. fakes) in addition to
    // the SusyNt ones.
    // Also note that this implementation relies completely on
    // ntsys2sys (where we do an explicit value check). This allows to
    // have, for example:
    // NtSys_EES_MAT_UP=4    and    Systematic::EESMATUP=14.
    int r=0;
    while(r<Susy::NtSys::SYS_UNKNOWN) { if(s==ntsys2sys(static_cast<Susy::NtSys::SusyNtSys>(r))) break; else r++;}
    return static_cast<Susy::NtSys::SusyNtSys>(r);
}
//-----------------------------------------

