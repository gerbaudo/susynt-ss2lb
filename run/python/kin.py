
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

def massBestZcandidate(l0, l1, otherLeps) :
    "Inputs are susy::wh::FourMom; return 0.0 if there is no candidate"
    l0, l1, otherLeps = addTlv(l0), addTlv(l1), [addTlv(l) for l in otherLeps]
    def lepFlavor(l) : return 'el' if l.isEl else 'mu' if l.isMu else 'other'
    def lepIsSeparatedFromOther(l, otherLeps=[], minDr=0.05) : return all(l.p4.DeltaR(ol.p4)>minDr for ol in otherLeps)
    def lepPairIsZcand(la, lb) :
        laFl, lbFl = lepFlavor(la), lepFlavor(lb)
        elOrMu = laFl in ['el','mu']
        sameFlavor = laFl==lbFl
        oppCharge  = la.charge*lb.charge < 0.0
        return elOrMu and sameFlavor and oppCharge
    def mll((la, lb)) : return (la.p4 + lb.p4).M()
    def deltaMZ0((la, lb)) : return abs(mll((la, lb)) - 91.2)
    otherLeps = filter(lambda l : lepIsSeparatedFromOther(l, [l0, l1]), otherLeps)
    validPairs = [(lh, ls) for lh in [l0, l1] for ls in otherLeps if lepPairIsZcand(lh, ls)]
    pairsSortedByBestM = sorted(validPairs, key=deltaMZ0)
    return mll(pairsSortedByBestM[0]) if len(pairsSortedByBestM) else 0.0
def thirdLepZcandidateIsInWindow(l0, l1, otherLeps, windowHalfwidth=20.0) :
    mZ0 = 91.2
    return fabs(massBestZcandidate(l0, l1, otherLeps) - mZ0) < windowHalfwidth

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
def computeMt2j(l0, l1, j0, j1, met, zeroMass=False, lspMass=0.0) :
    "As described in CMS-SUS-13-017"
    mt2_00 = computeMt2(l0+j0, l1+j1, met, zeroMass, lspMass)
    mt2_01 = computeMt2(l1+j0, l0+j1, met, zeroMass, lspMass)
    return min([mt2_00, mt2_01])
def computeMljj(l0, l1, j0, j1) :
    "Todo: extened combinatorics to N_j>2; good enough for now (we have few cases with >=3j)"
    jj = j0+j1
    dr0, dr1 = jj.DeltaR(l0), jj.DeltaR(l1)
    return (jj+l0).M() if dr0<dr1 else (jj+l1).M()
def computeMlj(l0, l1, j) :
    dr0, dr1 = j.DeltaR(l0), j.DeltaR(l1)
    return (j+l0).M() if dr0<dr1 else (j+l1).M()
