
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
