import math
from rootUtils import importRoot
r = importRoot()
import fakeUtils as fakeu

def histogram(variable, histoname, title=None, nx=None, xmin=None, xmax=None, ny=None, ymin=None, ymax=None,
              binsx=None, binsy=None) :
    twopi = +2.0*math.pi
    mpi = -math.pi
    ppi = +math.pi
    h = None
    attrs_1d = {
        'onebin'         :{'title':';; entries',                               'nx': 1, 'xmin':  0.5, 'xmax':  1.5},
        'njets'          :{'title':';N_{jets}; entries',                       'nx':10, 'xmin': -0.5, 'xmax':  9.5},
        'pt0'            :{'title':';p_{T,l0} [GeV]; entries/bin',             'nx':48, 'xmin':  0.0, 'xmax':240.0},
        'pt1'            :{'title':';p_{T,l1} [GeV]; entries/bin',             'nx':48, 'xmin':  0.0, 'xmax':240.0},
        'd_pt0_pt1'      :{'title':';p_{T,l0}-p_{T,l1} [GeV]; entries/bin',    'nx':24, 'xmin':  0.0, 'xmax':120.0},
        'eta0'           :{'title':';#eta_{l0}; entries/bin',                  'nx':26, 'xmin': -2.6, 'xmax': +2.6},
        'eta1'           :{'title':';#eta_{l1}; entries/bin',                  'nx':26, 'xmin': -2.6, 'xmax': +2.6},
        'phi0'           :{'title':';#phi_{l0} [rad]; entries/bin',            'nx':20, 'xmin':  mpi, 'xmax':  ppi},
        'phi1'           :{'title':';#phi_{l1} [rad]; entries/bin',            'nx':20, 'xmin':  mpi, 'xmax':  ppi},
        'mcoll'          :{'title':';m_{coll,l0,l1} [GeV]; entries/bin',       'nx':80, 'xmin':  0.0, 'xmax':400.0},
        'mcollcoarse'    :{'title':';m_{coll,l0,l1} [GeV]; entries/bin',       'nx':40, 'xmin':  0.0, 'xmax':400.0},
        'mll'            :{'title':';m_{l0,l1} [GeV]; entries/bin',            'nx':24, 'xmin':  0.0, 'xmax':240.0},
        'ptll'           :{'title':';p_{T,l0+l1} [GeV]; entries/bin',          'nx':24, 'xmin':  0.0, 'xmax':240.0},
        'met'            :{'title':';MET [GeV]; entries/bin',                  'nx':24, 'xmin':  0.0, 'xmax':240.0},
        'dphil0met'      :{'title':';|#Delta#phi(l0, met)| [rad]; entries/bin','nx':10, 'xmin':  0.0, 'xmax':  ppi},
        'dphil1met'      :{'title':';|#Delta#phi(l1, met)| [rad]; entries/bin','nx':10, 'xmin':  0.0, 'xmax':  ppi},
        'nsj'            :{'title':';N_{jets,20<pt<30};entries/bin',           'nx':10, 'xmin': -0.5, 'xmax':  9.5},
        'drl0csj'        :{'title':';#DeltaR(l0, j_{close,soft});entries/bin', 'nx':10, 'xmin':  0.0, 'xmax':  2.0},
        'drl1csj'        :{'title':';#DeltaR(l1, j_{close,soft});entries/bin', 'nx':10, 'xmin':  0.0, 'xmax':  2.0},
        'l0_d0Sig'       :{'title':';d_{0 sig, l0}; entries/bin',              'nx':40, 'xmin':-10.0, 'xmax':+10.0},
        'l0_z0Sin'       :{'title':';z_{0, l0}sin#theta; entries/bin',         'nx':40, 'xmin':-10.0, 'xmax':+10.0},
        'l0_etCone'      :{'title':';E_{T,cone, l0} [GeV]; entries/bin',       'nx':40, 'xmin':-10.0, 'xmax':+30.0},
        'l0_ptCone'      :{'title':';p_{T,cone, l0} [GeV]; entries/bin',       'nx':50, 'xmin':  0.0, 'xmax':+25.0},
        'l0_etConeCorr'  :{'title':';E_{T, cone, corr, l0}; entries/bin',      'nx':60, 'xmin':-10.0, 'xmax':+20.0},
        'l0_ptConeCorr'  :{'title':';p_{T, cone, corr, l0}; entries/bin',      'nx':50, 'xmin':  0.0, 'xmax':+25.0},
        'l1_d0Sig'       :{'title':';d_{0 sig, l1}; entries/bin',              'nx':40, 'xmin':-10.0, 'xmax':+10.0},
        'l1_z0Sin'       :{'title':';z_{0, l1}sin#theta; entries/bin',         'nx':40, 'xmin':-10.0, 'xmax':+10.0},
        'l1_etCone'      :{'title':';E_{T,cone, l1} [GeV]; entries/bin',       'nx':40, 'xmin':-10.0, 'xmax':+30.0},
        'l1_ptCone'      :{'title':';p_{T,cone, l1} [GeV]; entries/bin',       'nx':50, 'xmin':  0.0, 'xmax':+25.0},
        'l1_etConeCorr'  :{'title':';E_{T, cone, corr, l1}; entries/bin',      'nx':60, 'xmin':-10.0, 'xmax':+20.0},
        'l1_ptConeCorr'  :{'title':';p_{T, cone, corr, l1}; entries/bin',      'nx':50, 'xmin':  0.0, 'xmax':+25.0},
        'deta_jj'        :{'title':';|#Delta#eta_{jj}|; entries/bin',          'nx':50, 'xmin':  0.0, 'xmax':+5.2},
        'm_jj'           :{'title':';m_{jj} [GeV]; entries/bin',               'nx':50, 'xmin':  0.0, 'xmax':+1000},
        }
    attrs_2d = {
        'mcoll_vs_pt1'     :{'title':'; p_{T,l1} [GeV]; m_{coll,l0,l1} [GeV]',      'nx':240, 'xmin':0.0, 'xmax':240.0, 'ny':400, 'ymin':0.0, 'ymax':400.0},
        'pt0_vs_pt1'       :{'title':'; p_{T,l1} [GeV]; p_{T,l0} [GeV]',            'nx':240, 'xmin':0.0, 'xmax':240.0, 'ny':240, 'ymin':0.0, 'ymax':240.0},
        'met_vs_pt1'       :{'title':'; p_{T,l1} [GeV]; MET [GeV]',                 'nx':240, 'xmin':0.0, 'xmax':240.0, 'ny':24, 'ymin':0.0, 'ymax':240.0},
        'dphil0l1_vs_pt1'  :{'title':'; p_{T,l1} [GeV]; #Delta#phi(l0, l1) [rad]',  'nx':240, 'xmin':0.0, 'xmax':240.0, 'ny':40, 'ymin':0.0, 'ymax':4.0},
        'dphil0met_vs_pt1' :{'title':'; p_{T,l1} [GeV]; #Delta#phi(l0, met) [rad]', 'nx':240, 'xmin':0.0, 'xmax':240.0, 'ny':40, 'ymin':0.0, 'ymax':4.0},
        'dphil1met_vs_pt1' :{'title':'; p_{T,l1} [GeV]; #Delta#phi(l1, met) [rad]', 'nx':240, 'xmin':0.0, 'xmax':240.0, 'ny':40, 'ymin':0.0, 'ymax':4.0},
        'pt_eta' : {'title':';p_{T} [GeV]; #eta',
                    'nx':len(fakeu.ptBinEdges())-1,  'binsx':fakeu.ptBinEdges(),
                    'ny':len(fakeu.etaBinEdges())-1, 'binsy':fakeu.etaBinEdges(),}
        }


    if variable in attrs_1d:
        attrs = attrs_1d[variable]
        title = title if title else attrs['title']
        nx    = nx    if nx    else attrs['nx']
        xmin  = xmin  if xmin  else attrs['xmin']
        xmax  = xmax  if xmax  else attrs['xmax']
        h = r.TH1F(histoname, title, nx, xmin, xmax)
        h.Sumw2()
        h.SetDirectory(0)
    elif variable in attrs_2d:
        attrs = attrs_2d[variable]
        title = title if title else attrs['title']
        nx    = nx    if nx    else attrs['nx'   ]
        ny    = ny    if ny    else attrs['ny'   ]
        xmin  = xmin  if xmin  else attrs['xmin' ] if 'xmin' in attrs else None
        xmax  = xmax  if xmax  else attrs['xmax' ] if 'xmax' in attrs else None
        ymin  = ymin  if ymin  else attrs['ymin' ] if 'ymin' in attrs else None
        ymax  = ymax  if ymax  else attrs['ymax' ] if 'ymax' in attrs else None
        binsx = binsx if binsx else attrs['binsx'] if 'binsx' in attrs else None
        binsy = binsy if binsy else attrs['binsy'] if 'binsy' in attrs else None
        h = (r.TH2F(histoname, title, nx, binsx, ny, binsy) if binsx!=None and binsy!=None else
             r.TH2F(histoname, title, nx, xmin, xmax, ny, ymin, ymax))
        h.Sumw2()
        h.SetDirectory(0)
    else:
        print "settings.histogram: unknown variable, cannot book histogram"
    return h
