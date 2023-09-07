import numpy as np 
from configobj import ConfigObj
from distutils.util import strtobool

class ParameterReader():
    def __init__(self,file):
        self.fromfile = ConfigObj(file)
        self.dic = {}
        dic = self._update_dic()
        
    def _update_dic(self):
        self.dic['ipath'] = self.fromfile['Input']['ipath']
        self.dic['ifile'] = self.fromfile['Input']['ifile']
        self.dic['opath'] = self.fromfile['Output']['opath']
        self.dic['ofile'] = self.fromfile['Output']['ofile']
        self.dic['Rs'] = float(self.fromfile['Stellar']['Rs'])
        self.dic['stellar_file'] = self.fromfile['Stellar']['stellar_file']
        self.dic['save'] = bool(strtobool(self.fromfile['Output']['save']))
if __name__=='__main__':
    read = ParameterReader('myparam.par')
    print(read.dic.keys())
    for kk in read.dic.keys():
        print(kk,read.dic[kk])