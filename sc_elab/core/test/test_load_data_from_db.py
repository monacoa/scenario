


import sys
#import datetime as dtime

"""
import math
import numpy as np

#from sc_elab import numpy as np
from mdates import holidays
from mdates import daycount
from mdates import busdayrule

from scipy import optimize
from scipy.optimize import minimize
from scipy.optimize import fmin


import matplotlib.pyplot as plt

import datetime
import calendar

from datetime import datetime as dtime
from cookielib import DAYS

from sc_elab.core import funzioni_base as fb
from dateutil.relativedelta import relativedelta
from cmath import isnan
"""
from sc_elab.core import funzioni_bond_fitting as bf
from sc_elab.core import funzioni_base as fb
from sc_elab.core import SwpCurve as sw_curve
#from sc_elab.excel_hook import interfacce_base as ib

import time
import datetime

    


def FQ(label):
    print ('------------- FIN QUI TUTTO OK  %s ----------' %(label))
    sys.exit()



if __name__ == "__main__":
    

    #print 'AAA'
    
    #cc = sw_curve.CdsCurve()
    #ib.bootstrap_cds_from_xls()
    
    #FQ(999)
    
    
    cc = sw_curve.CdsCurve()
    opt_download = {}
    
    ref_date     = '2017-03-01'
    codeSeg      = '%LS'
    tipo_modello = 'LIN'
    valuta       = 'EUR'

    tipo_modello = 'LIN'
    opt_download['refDate'] = ref_date
    opt_download['valuta'] = valuta


    model_list    = []
    code_seg_list = []

    model_list    = ['SVE', 'NS', 'CIR', 'LIN']
    code_seg_list = ['%LS', '%DS', '%LFS', '%DFS']

    #model_list    = ['SVE']
    #code_seg_list = ['%LS']


    opt_download['codeSeg'] = codeSeg
    opt_download['tipo_modello'] = tipo_modello

    flag_break = False
    for mdl in model_list:
        for code_seg in code_seg_list:

            opt_download['codeSeg'] = code_seg
            opt_download['tipo_modello'] = mdl

            flag_download = cc.loadBenchDataFromDB(opt_download)

            print 'tipo_modello: ', mdl
            print 'code_seg: ', code_seg
            print 'flag_download: ', flag_download 
            print '-----------------------------------------------------'
            
            if (flag_download == 1):
                
                flag_break = True
                break

                
        if (flag_break == True):
                break   
            
            
        
    #print 'flag_download: ', flag_download
    #cc.loadBenchDataFromDB

    FQ(999)

