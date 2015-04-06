#!/bin/env python

# Pick objects from several files, and write them to one output file
#
# The objects to be picked+written and their names can be specified in
# the translation dictionary (see --print-default-dict)
#
# davide.gerbaudo@gmail.com
# 2014

import os
import pprint
import sys

import ROOT as r
r.gROOT.SetBatch(1)


def main():
    if len(sys.argv)<2:
        print 'usage:'
        print "{} --print-default-dict".format(sys.argv[0])
        print 'or'
        print "{} matrix_ingredients.py".format(sys.argv[0])
        return
    elif sys.argv[1]=='--print-default-dict':
        pprint.pprint(get_default_input_dict())
        return
    input_parameters = eval(open(sys.argv[1]).read())
    output_filename = input_parameters['output_filename']
    input_dict = input_parameters['input_dict']
    inputFiles = dict((fn, r.TFile.Open(fn, 'read')) for fn in input_dict.keys())
    if not all(v for k,v in inputFiles.iteritems()):
        raise IOError('Invalid input files:\n'+pprint.pformat(inputFiles))
    if os.path.exists(output_filename):
        raise IOError("Existing output filename: {}\nRefusing to overwrite".format(output_filename))
    out = r.TFile.Open(output_filename, 'recreate')
    out.cd()
    for fn, objList in input_dict.iteritems():
        inFile = inputFiles[fn]
        for oldName, newName in objList:
            newName = newName if newName else oldName
            obj = inFile.Get(oldName)
            if not obj:
                print "missing '{0}' from '{1}'".format(oldName, inFile.GetName())
                continue
            if hasattr(obj,'SetName') : obj.SetName(newName)
            obj.Write(newName)
    obj = r.TNamed('input_parameters', str(input_parameters))
    obj.Write()
    out.Close()
    print "objects written to %s"%out.GetName()


def get_default_input_dict():
    output_filename = './out/FakeMatrices/FakeMatrix_Apr_06.root'

    input0 = '/gdata/atlas/gerbaudo/wh/Susy2013_Nt_01_04/SusyMatrixMethod/data/FinalFakeHist_Apr_10.root'
    basedir = '/gdata/atlas/gerbaudo/wh/Susy2013_Nt_01_04_fake_dev/SusyTest0/run'
    input1 = basedir+'/'+'out/fake/weigtedmatrix_Jul_26/fake_matrices_el.root'
    input2 = basedir+'/'+'out/fake/weigtedmatrix_Jul_26/fake_matrices_mu.root'
    #('orig_name', 'new_name')  # none = keep the orig name
    input_dict = { input0 : [('mu_real_eff_CR_SSInc1j',    'mu_real_eff_emuInc')
                         ,('mu_real_eff2d_CR_SSInc1j', 'mu_real_eff2d_emuInc')
                         ,('el_real_eff_CR_SSInc1j',   'el_real_eff_emuInc')
                         ,('el_real_eff2d_CR_SSInc1j', 'el_real_eff2d_emuInc')
                         ,('el_real_up', None)
                         ,('el_real_down', None)
                         ,('mu_real_up', None)
                         ,('mu_real_down', None)
                         ,('el_HFLFerr', None)
                         ,('mu_HFLFerr', None)
                         ,('el_datamc', None)
                         ,('mu_datamc', None)
                         ,('el_region', None)
                         ,('mu_region', None)
                         ,('el_eta_sys', None)
                         ,('mu_eta_sys', None)
                         ],
               input1 : [('el_fake_pt_emu',      'el_fake_rate_emuInc')
                         ,('el_fake_pt_eta_emu', 'el_fake_rate2d_emuInc')
                         ],
               input2 : [('mu_fake_pt_emu',      'mu_fake_rate_emuInc')
                         ,('mu_fake_pt_eta_emu', 'mu_fake_rate2d_emuInc')
                         ],
               }
    return {'output_filename':'Fake_Matrix.root',
            'input_dict' : input_dict
            }

if __name__=='__main__':
    main()
