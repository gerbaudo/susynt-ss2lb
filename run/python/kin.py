
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
    twopi = 2.0*math.pi
    while phi < -pi : phi += twopi
    while phi > +pi : phi -= twopi
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
    metlv = tlv(met) # don't boost the input TLV, make a copy
    l0    = tlv(l0)
    l1    = tlv(l1)
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
#___________________________________________________________
def selection_formulas():
    "return a dictionary with the criteria for several selection regions"
    pt_req = 'l0_pt>35.0 and l1_pt>12.0 '
    pt_req += ' and abs(l0_eta)<2.4 and abs(l1_eta)<2.4' # tmp to be dropped, fixed upstream in ee8b767
    common_req_sr = (pt_req+
                     ' and n_cl_jets==0 and n_b_jets==0'+
                     ' and dphi_l1_met<0.7'+
                     ' and dphi_l0_l1>2.3 '+
                     ' and dphi_l0_met>2.5'+
                     ' and (l0_pt-l1_pt)>7.0')
    common_req_vr = (pt_req+
                     ' and n_jets==0'+
                     ' and (dphi_l1_met>0.7 or dphi_l0_l1<2.3 or dphi_l0_met<2.5 or (l0_pt-l1_pt)<7.0)')
    formulas = {
        'sr_emu' : 'l0_is_el and l1_is_mu and '+common_req_sr,
        'sr_mue' : 'l0_is_mu and l1_is_el and '+common_req_sr,
        'vr_emu' : 'l0_is_el and l1_is_mu and '+common_req_vr,
        'vr_mue' : 'l0_is_mu and l1_is_el and '+common_req_vr,
        'sr_emu_mue' : '(is_emu or is_mue) and '+common_req_sr,
        'sr_ee' : 'l0_is_el and l1_is_el and '+common_req_sr,
        'sr_mumu' : 'l0_is_mu and l1_is_mu and '+common_req_sr,
        'vr_emu_mue' : '(is_emu or is_mue) and '+common_req_vr,
        'vr_ee' : 'l0_is_el and l1_is_el and '+common_req_vr,
        'vr_mumu' : 'l0_is_mu and l1_is_mu and '+common_req_vr,
        # 'pre_emu' : 'is_emu and '+pt_req,
        # 'pre_mue' : 'is_mue and '+pt_req,
        # 'pre_emu_mue' : '(is_emu or is_mue) and '+pt_req,
        }
    os_expr = '(is_opp_sign)'
    # ss_expr = '(is_same_sign or is_qflippable)'
    ss_expr = '(is_same_sign)'
    formulas = dict([(k+'_'+ssos, v+' and '+ssos_expr)
                     for k, v in formulas.iteritems()
                     for ssos, ssos_expr in [('ss', ss_expr), ('os', os_expr)]])
    # # symmetric selection
    # pt_sym_req = 'l0_pt>20.0 and l1_pt>20.0'
    # for lf, lf_expr in [('emu', 'is_emu'), ('mue', 'is_mue'), ('emu_mue', '(is_emu or is_mue)')]:
    #     for ssos, ssos_expr in [('ss', ss_expr), ('os', os_expr)]:
    #         formulas['sym_'+lf+'_'+ssos] = pt_sym_req+' and '+lf_expr+' and '+ssos_expr
    # validation region used by Matt in the 2L paper, see sec6.4 ATL-COM-PHYS-2012-1808
    # formulas_vrss_btag = 'num_b_jets==1 and et_miss_rel>50.0 and abs(m_ll-91.2)>10.0 if is_ee else True) and ((m_ll<90.0 or m_ll>120) if is_mumu else True)'
    # formulas['vrss_btag'] = formulas_vrss_btag

    formulas['vr_emu_razor_ss'] = '(is_emu or is_mue) and mdr>20.0 and '+ss_expr
    formulas['vr_ee_razor_ss'] = 'is_ee and mdr>20.0 and '+ss_expr
    formulas['vr_mumu_razor_ss'] = 'is_mumu and mdr>20.0 and '+ss_expr
    # real extraction formulas; trig_match is already applied when making the ntuples
    formulas['ext_mumu_zpeak_probe_mu'] = ('is_mumu and is_opp_sign and abs(m_ll-91.2)<10.0'
                                           ' and tag_is_mu and tag_is_tight and probe_is_mu')
    formulas['ext_ee_zpeak_probe_el'] = ('is_ee and is_opp_sign and abs(m_ll-91.2)<10.0'
                                         ' and tag_is_el and tag_is_tight and probe_is_el')
    # fake extraction formulas; trig_match is already applied when making the ntuples
    formulas['ext_mumu_ss_probe_mu'] = 'is_mumu and is_same_sign and tag_is_mu and tag_is_tight and probe_is_mu'
    formulas['ext_emu_mue_ss_probe_el'] = '(is_emu or is_mue) and is_same_sign and tag_is_mu and tag_is_tight and probe_is_el'
    formulas['ext_mumu_ss'] = 'is_mumu and is_same_sign'
    formulas['ext_emu_mue_ss'] = '(is_emu or is_mue) and is_same_sign'
    formulas['ext_emu_pt0_40_ss']= 'is_emu and is_same_sign and l0_pt>40.0'
    formulas['ext_mue_pt0_40_ss']= 'is_mue and is_same_sign and l0_pt>40.0'
    formulas['ext_mumu_pt0_40_ss']= 'is_mumu and is_same_sign and l0_pt>40.0'
    # dbg low pt
    formulas['sr_mue_os_low_pt1_15'] = (formulas['sr_mue_os']+' and l1_pt<15.0')
    formulas['sr_mue_os_hi_pt1_15'] = (formulas['sr_mue_os']+' and l1_pt>15.0')
    formulas['sr_emu_os_vbf'] = (formulas['sr_emu_os'].replace('n_jets==0', 'n_jets==2')+' and deta_jj>3.5 and n_b_jets==0')
    formulas['sr_mue_os_vbf'] = (formulas['sr_mue_os'].replace('n_jets==0', 'n_jets==2')+' and deta_jj>3.5 and n_b_jets==0')
    formulas['sr_emu_ss_vbf'] = (formulas['sr_emu_ss'].replace('n_jets==0', 'n_jets==2')+' and deta_jj>3.5 and n_b_jets==0')
    formulas['sr_mue_ss_vbf'] = (formulas['sr_mue_ss'].replace('n_jets==0', 'n_jets==2')+' and deta_jj>3.5 and n_b_jets==0')
    # formulas for Aielet
    common_req_cr = ('is_opp_sign and'
                     ' l0_pt>45.0 and l1_pt>12.0 and'
                     ' dphi_l1_met < 0.7 and '
                     ' (dphi_l0_l1 < 2.3 or dphi_l0_met < 2.5 ) and'
                     ' n_jets==0')
    common_req_cr40 = common_req_cr.replace('l0_pt>45.0', 'l0_pt>40.0')
    common_req_crttbar = ('is_opp_sign and'
                          ' l0_pt>45.0 and l1_pt>12.0 and'
                          ' dphi_l1_met < 0.7 and '
                          ' (dphi_l0_l1 < 2.3 or dphi_l0_met < 2.5 ) and'
                          ' n_cl_jets==0 and'
                          ' n_b_jets>=1 and'
                          ' n_f_jets==0')

    formulas['cr_emu_os'      ] = ('l0_is_el and l1_is_mu and '+common_req_cr)
    formulas['cr_mue_os'      ] = ('l0_is_mu and l1_is_el and '+common_req_cr)
    formulas['cr_40_emu_os'   ] = ('l0_is_el and l1_is_mu and '+common_req_cr40)
    formulas['cr_40_mue_os'   ] = ('l0_is_mu and l1_is_el and '+common_req_cr40)
    formulas['cr_ttbar_emu_os'] = ('l0_is_el and l1_is_mu and '+common_req_crttbar)
    formulas['cr_ttbar_mue_os'] = ('l0_is_mu and l1_is_el and '+common_req_crttbar)

    common_req_sr_jets = ('is_opp_sign'+
                          ' and l0_pt>35.0 and l1_pt>12.0'+
                          ' and n_cl_jets>0 and n_b_jets==0'+
                          ' and dphi_l1_met<0.5'+
                          ' and dphi_l0_l1>1.0 '+
                          ' and dphi_l0_met>1.0'+
                          ' and (l0_pt-l1_pt)>1.0')
    formulas['sr_emu_os_jets'] = ('is_emu and '+common_req_sr_jets)
    formulas['sr_mue_os_jets'] = ('is_mue and '+common_req_sr_jets) # todo: update also cr*jets, vr*jets
    formulas['sr_emu_ss_jets'] = formulas['sr_emu_os_jets'].replace('is_opp_sign', 'is_same_sign')
    formulas['sr_mue_ss_jets'] = formulas['sr_mue_os_jets'].replace('is_opp_sign', 'is_same_sign')

    formulas['cr_mue_os_jets'] = (formulas['cr_mue_os'].replace('n_jets==0', 'n_jets>0')+' and n_bf_jets==0')
    formulas['cr_emu_os_jets'] = (formulas['cr_emu_os'].replace('n_jets==0', 'n_jets>0')+' and n_bf_jets==0')
    formulas['vr_mue_os_jets'] = (formulas['vr_mue_os'].replace('n_jets==0', 'n_jets>0')+' and n_bf_jets==0')
    formulas['vr_emu_os_jets'] = (formulas['vr_emu_os'].replace('n_jets==0', 'n_jets>0')+' and n_bf_jets==0')

    formulas['sr_ttbar_emu_os'] = ('is_opp_sign and l0_is_el and l1_is_mu and '+common_req_sr
                                   ).replace('n_jets==0',
                                             'n_b_jets>=1 and n_f_jets==0')
    formulas['sr_ttbar_mue_os'] = formulas['sr_ttbar_emu_os'].replace('l0_is_el and l1_is_mu',
                                                                      'l0_is_mu and l1_is_el')
    common_req_cr2 = ('is_opp_sign'
                      ' and '+pt_req+
                      ' and n_jets==0'+
                      ' and dphi_l1_met<0.7'
                      ' and dphi_l0_l1<2.3'
                      ' and dphi_l0_met>2.5')
    formulas['cr2_emu_os'] = ('l0_is_el and l1_is_mu and '+common_req_cr2)
    formulas['cr2_mue_os'] = ('l0_is_mu and l1_is_el and '+common_req_cr2)

    common_cr34 = ('is_opp_sign and '+pt_req+' and n_jets==0')
    common_cr3_os = (common_cr34+' and     dphi_l1_met<0.4 and (dphi_l0_l1<2.3 or dphi_l0_met<2.5)')
    common_cr4_os = (common_cr34+' and 0.4<dphi_l1_met<0.7 and (dphi_l0_l1<2.3 or dphi_l0_met<2.5)')

    formulas['cr3_emu_os'] = ('l0_is_el and l1_is_mu and '+common_cr3_os)
    formulas['cr3_mue_os'] = ('l0_is_mu and l1_is_el and '+common_cr3_os)
    formulas['cr4_emu_os'] = ('l0_is_el and l1_is_mu and '+common_cr4_os)
    formulas['cr4_mue_os'] = ('l0_is_mu and l1_is_el and '+common_cr4_os)

    formulas['base_emu_os'] = ('l0_is_el and l1_is_mu and is_opp_sign and '+pt_req+' and n_jets==0')
    formulas['base_mue_os'] = ('l0_is_mu and l1_is_el and is_opp_sign and '+pt_req+' and n_jets==0')

    common_req_cr_old = ('is_opp_sign'+
                         ' and '+pt_req+
                         ' and n_jets==0'+
                         ' and (dphi_l1_met>0.7 or dphi_l0_l1<2.3 or dphi_l0_met<2.5)')

    formulas['cr_old_emu_os'] = ('l0_is_el and l1_is_mu and '+common_req_cr_old)
    formulas['cr_old_mue_os'] = ('l0_is_mu and l1_is_el and '+common_req_cr_old)

    formulas['sr_40_emu_os']       = formulas['sr_emu_os'      ].replace('l0_pt>45.0', 'l0_pt>40.0')
    formulas['sr_40_mue_os']       = formulas['sr_mue_os'      ].replace('l0_pt>45.0', 'l0_pt>40.0')
    formulas['sr_40_ttbar_emu_os'] = formulas['sr_ttbar_emu_os'].replace('l0_pt>45.0', 'l0_pt>40.0')
    formulas['sr_40_ttbar_mue_os'] = formulas['sr_ttbar_mue_os'].replace('l0_pt>45.0', 'l0_pt>40.0')
    formulas['sr_40_ttbar_emu_os'] = formulas['sr_ttbar_emu_os'].replace('l0_pt>45.0', 'l0_pt>40.0')
    formulas['sr_40_ttbar_mue_os'] = formulas['sr_ttbar_mue_os'].replace('l0_pt>45.0', 'l0_pt>40.0')
    formulas['cr_40_ttbar_emu_os'] = formulas['cr_ttbar_emu_os'].replace('l0_pt>45.0', 'l0_pt>40.0')
    formulas['cr_40_ttbar_mue_os'] = formulas['cr_ttbar_mue_os'].replace('l0_pt>45.0', 'l0_pt>40.0')

    common_req_cr5_os = ('is_opp_sign and'
                         ' l0_pt>45.0 and l1_pt>12.0 and'
                         ' dphi_l1_met < 0.7 and '
                         ' (dphi_l0_l1 < 2.3 or dphi_l0_met < 2.5 or (l0_pt-l1_pt)<7.0) and'
                         ' n_jets==0')
    formulas['cr5_emu_os'] = ('l0_is_el and l1_is_mu and '+common_req_cr5_os)
    formulas['cr5_mue_os'] = ('l0_is_mu and l1_is_el and '+common_req_cr5_os)

    return formulas
#___________________________________________________________
