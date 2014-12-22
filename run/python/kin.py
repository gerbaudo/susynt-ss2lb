
# kinematic utilities
#
# davide.gerbaudo@gmail.com
# Jan 2014, original version from SusyTest0
# Oct 2014, modified version for SusyntHlfv

import array
import math
import os

from rootUtils import importRoot
r = importRoot()
try:
    rootcoredir = os.environ['ROOTCOREDIR']
    r.gROOT.LoadMacro(rootcoredir+'/scripts/load_packages.C+')
    r.load_packages()
except KeyError :
    print "undefined ROOTCOREDIR: the functions involving susy::wh::FourMom and mt2 will not work"

fabs, cos, sin, pi, sqrt = math.fabs, math.cos, math.sin, math.pi, math.sqrt

def phi_mpi_pi(phi) :
    pi = math.pi
    while phi < -pi : phi += pi
    while phi > +pi : phi -= pi
    return phi

tlv = r.TLorentzVector
def FourMom2TLorentzVector(fm) :
    l = tlv()
    l.SetPxPyPzE(fm.px, fm.py, fm.pz, fm.E)
    return l
def addTlv(l) :
    if not hasattr(l, 'p4') : l.p4 = FourMom2TLorentzVector(l)
    return l
def computeMt(lep, met) :
    return sqrt(2.0 * lep.Pt() * met.Et() *(1.0-cos(lep.DeltaPhi(met))))
def computeHt(met, leptsJets=[]) :
    return sum(o.Pt() for o in leptsJets) + met.Et()
def computeMetRel(met, leptsJets=[]) :
    minDphi = min([0.5*pi]+[fabs(met.DeltaPhi(o)) for o in leptsJets])
    return met.Et()*sin(minDphi)

def getDilepType(fmLep0, fmLep1) :
    "Given two susy::wh::FourMom, return ee/em/mm"
    def FourMom2LepType(fm) : return 'e' if fm.isEl else 'm' if fm.isMu else None
    dilepType = ''.join(sorted(FourMom2LepType(l) for l in [fmLep0, fmLep1])) # sort -> em instead of me
    return dilepType

def computeMt2(a, b, met, zeroMass, lspMass) :
    mt2 = r.mt2_bisect.mt2()
    pa    = array.array( 'd', [0.0 if zeroMass else a.M(), a.Px(), a.Py() ] )
    pb    = array.array( 'd', [0.0 if zeroMass else b.M(), b.Px(), b.Py() ] )
    pmiss = array.array( 'd', [0.0, met.Px(), met.Py() ] )
    mt2.set_momenta(pa, pb, pmiss)
    mt2.set_mn(lspMass)
    return mt2.get_mt2()

def computeCollinearMassLepTau(l0, l1, met):
    sqrt, cos, cosh = math.sqrt, math.cos, math.cosh
    return sqrt(2.0 * l0.Pt() *
                (l1.Pt() + met.Et()) *
                (cosh(l0.Eta() - l1.Eta()) - cos(l0.DeltaPhi(l1))))

def computeRazor(l0, l1, met, throw_on_neg_sqrt=False):
    """
    razor variables from hep-ph/1310.4827.
    Inputs are TLorentzVector objects
    """
    metlv = met
    l0    = l0
    l1    = l1
    # lab frame
    vBETA_z = (l0+l1).Vect()*r.Double(1./(l0.E()+l1.E()))
    vBETA_z.SetX(0.0)
    vBETA_z.SetY(0.0)
    l0.Boost(-vBETA_z)
    l1.Boost(-vBETA_z)
    pT_CM = (l0+l1).Vect() + metlv.Vect()
    pT_CM.SetZ(0.0)
    ll = l0+l1
    SHATR = sqrt( 2.*(ll.E()*ll.E() - ll.Vect().Dot(pT_CM)
                      + ll.E()*sqrt( ll.E()*ll.E() + pT_CM.Mag2() - 2.*ll.Vect().Dot(pT_CM) )))
    vBETA_T_CMtoR = pT_CM * r.Double(1./sqrt(pT_CM.Mag2() + SHATR*SHATR))
    l0.Boost(-vBETA_T_CMtoR)
    l1.Boost(-vBETA_T_CMtoR)
    ll.Boost(-vBETA_T_CMtoR)
    # R-frame
    dphi_LL_vBETA_T = fabs((ll.Vect()).DeltaPhi(vBETA_T_CMtoR))
    dphi_L1_L2 = fabs(l0.Vect().DeltaPhi(l1.Vect()))
    vBETA_R = (l0.Vect() - l1.Vect())*r.Double(1./(l0.E()+l1.E()))
    try:
        gamma_R = 1./sqrt(1.-vBETA_R.Mag2())
    except ValueError:
        if throw_on_neg_sqrt:
            print 'trying to sqrt a negative value: ',(1.-vBETA_R.Mag2())
            raise
    dphi_vBETA_R_vBETA_T = fabs(vBETA_R.DeltaPhi(vBETA_T_CMtoR))
    l0.Boost(-vBETA_R)
    l1.Boost(vBETA_R)
    # R+1 frame
    MDELTAR = 2.*l0.E()
    costhetaRp1 = l0.Vect().Dot(vBETA_R)/(l0.Vect().Mag()*vBETA_R.Mag())
    return dphi_LL_vBETA_T, MDELTAR
