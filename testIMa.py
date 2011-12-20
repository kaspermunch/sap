
import sys

from SAP.PostAnalysis import IMa
from SAP import Fasta

class Bunch(object):
    """
    Generic class for lumping attributes.
    """
    def __init__(self, **keywords):
        self.__dict__.update(keywords)

    def has(self, keyword):
        return self.__dict__.has_key(keywords)

    def __str__(self):
        s = 'Bunch: '
        for k, v in self.__dict__.items():
            s += "%s: %s  " % (k, v)
        return s

options = Bunch()
pwd = "/Users/kasper/projects/sap/analyses/imatest/prelim/"
options.homologcache = pwd + "homologcache"
options.treestatscache = pwd + "treestatscache"

ass = IMa.Assignment(options)

"prelim_BOTW283_05_USNM_621383_Dendroica_kirtlandii_"


queryFastaRecord = Fasta.Record(title="HCBR116_03_1B_1613_Dendroica_caerulescens_", sequence='GGAATGGTAGGTACCGCCCTAAGCCTCCTCATTCGAGCAGAACTAGGCCAACCTGGAGCCCTTCTAGGAGACGACCAAGTCTACAACGTAGTTGTCACGGCCCATGCTTTCGTAATAATTTTCTTTATAGTTATGCCGATTATAATCGGAGGGTTCGGAAACTGACTAGTCCCCCTAATAATCGGAGCCCCAGACATAGCGTTCCCACGAATAAACAACATAAGCTTCTGACTACTCCCACCATCATTCCTTCTCCTCCTAGCATCCTCCACAGTTGAAGCAGGTGTAGGCACAGGCTGAACAGTATACCCCCCACTAGCTGGCAACTTGGCCCACGCCGGAGCCTCAGTCGACCTCGCAATCTTCTCCCTACACCTAGCCGGTATTTCCTCAATCCTCGGGGCAATCAACTTCATTACAACAGCAATTAATATGAAACCTCCTGCCCTCTCACAATACCAAACCCCACTATTCGTCTGATCAGTCCTAATCACTGCAGTCCTCTTACTCCTTTCCCTTCCAGTTCTAGCTGCAGGAATCACAATGCTCCTCACAGACCGCAACCTTAACACCACATTCTTCGACCCTGCTGGAGGAGGAGATCCCGTCCTATATCAACATCTCTTCTGATTCTTTGGTCACCCAGAAGTTTACATCCTAATCCTC')

ass.run(queryFastaRecord, 'prelim')




# cmd = "im -h"
# outputPrefix = '/tmp/kasper'
# 
# arguments = cmd.split(' ')
# 
# retval = IM.runprogram(arguments, outputPrefix)
