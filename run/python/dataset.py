# Class and functions to handle the datasets (data and simulated ones)

# davide.gerbuado@gmail.com
# Jun 2014

import glob
import os
import re
import tempfile
import unittest
import utils

class Dataset(object):
    """
    A dataset with the attributes that can be extracted from its name.
    Notes on implementation:
    - instantiate regex once, as static datamembers
    - define parsing methods as @classmethod, so that the class can be
      used as a base class and extended to parse additional attributes
    """
    attributes_data = ['period', 'stream']
    attributes_simu = ['dsid', 'dsname']
    stem_data = 'data12_8TeV'
    stem_simu = 'mc12_8TeV'
    # static regex used to tokenize attributes
    regex_data = re.compile('group\.phys-susy\.data12\_8TeV\.'
                            'period(?P<period>[A-Z])\.'
                            'physics\_(?P<stream>\w+)\.'
                            'PhysCont')
    regex_simu = re.compile('mc12\_8TeV\.'
                            '(?P<dsid>\d+)\.'
                            '(?P<dsname>.*)')
    verbose_parsing = False
    @classmethod
    def parse_datasets_from_file(cls, filename=''):
        datasets = []
        def clean_line(l) : return l.strip()
        def is_empty(l) : return not l
        def is_comment(l) : return l.startswith('#')
        with open(filename) as input_file:
            for line in input_file.readlines():
                line = clean_line(line)
                if is_empty(line) or is_comment(line) : continue
                ds = Dataset(line)
                if ds.is_valid : datasets.append(ds)
        return datasets
    @classmethod
    def parse_files_in_dir(cls, dirname='', ext='.txt'):
        files = glob.glob(os.path.join(dirname, '*'+ext))
        return [d for f in files for d in cls.parse_datasets_from_file(f)]
    @classmethod
    def parse_dataset_name(cls, dsname=''):
        return (cls.regex_data.search(dsname) if cls.is_data(dsname) else
                cls.regex_simu.search(dsname) if cls.is_simulation(dsname) else
                None)
    @classmethod
    def is_data(cls, dsname) : return cls.stem_data in dsname
    @classmethod
    def is_simulation(cls, dsname) : return cls.stem_simu in dsname
    @property
    def attributes(self):
        return (Dataset.attributes_data if Dataset.is_data(self.name) else
                Dataset.attributes_simu if Dataset.is_simulation(self.name) else
                [])
    def __init__(self, name):
        """
        Behaviour: always set the name; set the attributes when we can
        parse. The idea is that one might want to use this as a base
        class, but define a different type of parsing.
        """
        self.name = name
        self.is_valid = False
        self.is_data = Dataset.is_data(name)
        self.is_simulation = Dataset.is_simulation(name)
        atts = self.attributes
        match = Dataset.parse_dataset_name(name)
        if match:
            kwargs = dict([(g, match.group(g)) for g in atts])
            self.__dict__.update(kwargs)
            self.is_valid = True
        elif Dataset.verbose_parsing:
            print "Cannot parse '%s'; attributes not set"%name
    def build_filelist(self, base_input_dir, out_dir, verbose=False):
        success = False
        getCommandOutput = utils.getCommandOutput
        cmd = 'find '+base_input_dir+' -maxdepth 1 -name "*'+self.name+'*"'
        directories = getCommandOutput(cmd)['stdout'].split()
        dest_filename = out_dir+'/'+self.name+'.txt'
        if len(directories)==1:
            directory = directories[0]
            if verbose and os.path.exists(dest_filename) : print "overwriting %s"%dest_filename
            cmd = 'ls '+directory+'/*.root* > '+dest_filename
            out = getCommandOutput(cmd)
            success = out['returncode']==0
        elif len(directories)>1:
            print "build_filelist: ambiguous output for '%s'"%self.name
            print '\n\t : '.join(directories)
            print "you need to generate it by hand: %s"%('ls '+base_input_dir+'/[dir]/*.root* > '+self.name+'.txt')
        else:
            print "build_filelist: missing '%s' from %s"%(self.name, base_input_dir)
        return success

# testing
#_____________________________________
class testDataset(unittest.TestCase):
    def testValidAttributes(self):
        "check that we can access the parsed attributes"
        knownValues = [('group.phys-susy.data12_8TeV.periodL.physics_Egamma.PhysCont',
                        {'period':'L',
                         'stream':'Egamma'}),
                       ('mc12_8TeV.160155.PowhegPythia8_AU2CT10_ggH125_ZZ4lep',
                        {'dsid':'160155',
                         'dsname':'PowhegPythia8_AU2CT10_ggH125_ZZ4lep'}),
                       ]
        for dsname, values in knownValues:
            dataset = Dataset(dsname)
            attributes = dataset.attributes
            for att in attributes:
                self.assertEqual(values[att], getattr(dataset, att))
    def testIsValid(self):
        "check that we can construct Dataset object with and without attributes"
        dsnames = ['group.phys-susy.data12_8TeV.periodL.physics_Egamma.PhysCont',
                   'mc12_8TeV.160155.PowhegPythia8_AU2CT10_ggH125_ZZ4lep',
                   'mc12_8TeV.160155PowhegPythia8_AU2CT10_ggH125_ZZ4lep', # missing . after dsid
                   'group.phys-susy.data12_8TeV.periodLabfkjhs.physics_Egamma.PhysCont_typo', # typo at the end
                   ]
        datasets = [Dataset(n) for n in dsnames]
        num_valid = sum(d.is_valid for d in datasets)
        num_invalid = sum(not d.is_valid for d in datasets)
        self.assertEqual(num_valid, 2)
        self.assertEqual(num_invalid, 2)
    def testParsing(self):
        "check that we can parse a file"
        dummy_file_content = """
# just a comment
group.phys-susy.data12_8TeV.periodL.physics_Egamma.PhysCont
mc12_8TeV.160155.PowhegPythia8_AU2CT10_ggH125_ZZ4lep
# now an invalid dataset
mc12_8TeV.16bazzinga
and then an invalid line
"""
        with tempfile.NamedTemporaryFile() as temp:
            temp.write(dummy_file_content)
            temp.flush()
            datasets = Dataset.parse_datasets_from_file(temp.name)
            num_valid = sum(d.is_valid for d in datasets)
            num_invalid = sum(not d.is_valid for d in datasets)
            self.assertEqual(num_valid, 2)
            self.assertEqual(num_invalid, 0) # they're dropped by parse_datasets_*
            
            
#_____________________________________
if __name__ == "__main__":
    # Dataset.verbose_parsing = True # toggle verbose if you see errors
    unittest.main()
