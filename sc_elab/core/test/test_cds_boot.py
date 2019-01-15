


import sys
#import datetime as dtime
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


def FQ(label):
    print ('------------- FIN QUI TUTTO OK  %s ----------' %(label))
    sys.exit()


def load_data_fromFile(inputRetFile):

    fin = open(inputRetFile, 'r')
    listInput = fin.readlines()

    nodoVec = []
    matVec = []
    valVec = []
    usaNodoVec = []
    tipoSegmentoVec = []
    
    data_out = {}
    
    for i in range(1, len(listInput)):
    
        line_splittedTmp = listInput[i].split("\t")        

        nodoTmp         = str(line_splittedTmp[0])
        valTmp          = float(line_splittedTmp[2])
        usaNodoTmp      = str(line_splittedTmp[3])
        tipoSegmentoTmp = str(line_splittedTmp[4].split("\n")[0]) 

        matDateTmp      = line_splittedTmp[1].split('-')
        
        try:
            
            int(matDateTmp[0])
            
        except:
            
            matDateTmp      = line_splittedTmp[1].split('/')


        matDateTmp_dd = int(matDateTmp[0])
        matDateTmp_mm = int(matDateTmp[1])
        matDateTmp_yy = int(matDateTmp[2])
                
        dateTmp = dtime.date(matDateTmp_yy, matDateTmp_mm, matDateTmp_dd)

        nodoVec.append(nodoTmp)
        matVec.append(dateTmp)
        valVec.append(valTmp)
        usaNodoVec.append(usaNodoTmp)
        tipoSegmentoVec.append(tipoSegmentoTmp)

        
    data_out['Nodo'] = nodoVec
    data_out['MatDate'] = matVec
    data_out['ValoreNodo'] = valVec
    data_out['UsaNodo'] = usaNodoVec
    data_out['TipoSegmento'] = tipoSegmentoVec
        
    return data_out


def load_cds_data_fromFile(inputRetFile):

    fin = open(inputRetFile, 'r')
    listInput = fin.readlines()

    nodoVec = []
    matVec = []
    valVec = []
    
    data_out = {}
    
    for i in range(1, len(listInput)):
    
        line_splittedTmp = listInput[i].split("\t") 

        matTimeTmp       = float(line_splittedTmp[1])             
        nodoTmp          = str(line_splittedTmp[0])
        valTmp           = float(line_splittedTmp[2])
        
        matVec.append(matTimeTmp)
        nodoVec.append(nodoTmp)
        valVec.append(valTmp)

        
    data_out['Nodo'] = nodoVec
    data_out['MatTimes'] = matVec
    data_out['ValoreNodo'] = valVec
        
    return data_out

def load_simple_curves_fromFile(inputRetFile):

    basis = '6M'
    fin = open(inputRetFile, 'r')
    listInput = fin.readlines()

    nodoVec = []
    matVec = []
    valVec = []
    
    data_out = {}
    
    for i in range(1, len(listInput)):
    
        line_splittedTmp = listInput[i].split("\t")
        
        

        valTmp          = float(line_splittedTmp[1])

        matDateTmp      = line_splittedTmp[0].split('-')
        
        try:
            
            int(matDateTmp[0])
            
        except:
            
            matDateTmp      = line_splittedTmp[0].split('/')


        matDateTmp_dd = int(matDateTmp[0])
        matDateTmp_mm = int(matDateTmp[1])
        matDateTmp_yy = int(matDateTmp[2])
                
        dateTmp = datetime.date(matDateTmp_yy, matDateTmp_mm, matDateTmp_dd)

        matVec.append(dateTmp)
        valVec.append(valTmp)
        
    data_out['MatDate'] = matVec
    data_out['ValoreNodo'] = valVec
    data_out['Basis'] = basis
    
    return data_out



def load_zc_data_fromFile(inputFile):

    fin = open(inputFile, 'r')
    listInput = fin.readlines()

    matTVec = []
    matDVec  = []
    valVec  = []
    
    data_out = {}
    
    for i in range(1, len(listInput)):
    
        line_splittedTmp = listInput[i].split("\t")        

        valTmp          = float(line_splittedTmp[1])
        matDateTmp      = line_splittedTmp[0].split('-')
        
        try:
            
            int(matDateTmp[0])
            
        except:
            
            matDateTmp      = line_splittedTmp[0].split('/')

        matDateTmp_dd = int(matDateTmp[0])
        matDateTmp_mm = int(matDateTmp[1])
        matDateTmp_yy = int(matDateTmp[2])
        
        dateTmp = dtime.date(matDateTmp_yy, matDateTmp_mm, matDateTmp_dd)
        
        if (i == 1):date0 = dateTmp

        ggTmp = (dateTmp - date0).days
        timeTmp = ggTmp/365.2425
        
        matDVec.append(dateTmp)
        matTVec.append(timeTmp)
        valVec.append(valTmp)

    """ 
    data_out['MatDate']  = matDVec[1:]
    data_out['MatTimes'] = matTVec[1:]
    data_out['ValoreZC'] = valVec[1:]
    """
    
    data_out['MatDate']  = matDVec
    data_out['MatTimes'] = matTVec
    data_out['ValoreZC'] = valVec


    return data_out






def dumpDataOnFile(cds_times_n, surv_prob, fileOutSurv):
    
        
    n_times = len(cds_times_n)
    fileO = open(fileOutSurv, 'w')

    fileO.write("Data\t Valore\n")
    
    for i in range(0, n_times):

        dateTmp  =  cds_times_n[i]
        valueTmp =  surv_prob[i]   

        fileO.write("%s\t %s\n" %(dateTmp, valueTmp))
             

    fileO.close()
 
    
    

def convert2Df(data_raw_ref, base_ref):
    
    
    data_discount = []
    n_dates = len(data_raw_ref['MatDate'])
    
    for i in range(0, n_dates):
        
        
        dateTmp = data_raw_bench['MatDate'][i]
        valTmp = data_raw_bench['ValoreNodo'][i]
        
        if (i == 0):date0 = dateTmp
        
        ggTmp = (dateTmp - date0).days
        timeTmp = ggTmp/365.2425
        
        dfTmp = np.exp(-timeTmp*valTmp/base_ref)
        
        
        data_discount.append(dfTmp)
      

    data_raw_ref['DiscountFactor'] = data_discount
    
    return data_raw_ref


def set_model_prms(data_raw_bench):
    
    data_raw_bench['Model'] = 2
    
    dict_params = {}
    
    dict_params['const1'] = [1.0]
    dict_params['const2'] = [1.0]
    dict_params['beta0']  = [0.03]
    dict_params['beta1']  = [0.0]
    dict_params['beta2']  = [0.0]
    dict_params['beta3']  = [0.0]
    
    
    data_raw_bench['prms'] = dict_params

    
    return data_raw_bench

#def listaDateOUT

def dumpResultsOnFile(res, fileOutRes):
    

    n_times = len(res['zcSpread'])
    fileO = open(fileOutRes, 'w')

    fileO.write("Date\t    Survival\t    MargDef\t    HazardRate\t     zeroSpread\t    pySpread\n")
    
    for i in range(0, n_times):

        zspreadTmp = res['zcSpread'][i]
        pySpreadTmp =  res['pySpread'][i]

        dataOutTmp = res['outputDates'][i]

        survTmp =  res['survProbCum'][i]
        mdTmp = res['marginalDefault'][i]
        hrTmp =  res['hazardRate'][i]

        fileO.write("%s\t %s\t %s\t %s\t %s\t %s\n" %(dataOutTmp, survTmp, mdTmp, hrTmp, zspreadTmp, pySpreadTmp))
             
    fileO.close()
 
def load_zc_data_fromFile(inputFile):

    fin = open(inputFile, 'r')
    listInput = fin.readlines()

    matTVec = []
    matDVec  = []
    valVec  = []
    
    data_out = {}
    
    for i in range(1, len(listInput)):
    
        line_splittedTmp = listInput[i].split("\t")
        

        valTmp          = float(line_splittedTmp[1])
        matDateTmp      = line_splittedTmp[0].split('-')
        
        try:
            
            int(matDateTmp[0])
            
        except:
            
            matDateTmp      = line_splittedTmp[0].split('/')

        matDateTmp_dd = int(matDateTmp[0])
        matDateTmp_mm = int(matDateTmp[1])
        matDateTmp_yy = int(matDateTmp[2])
        
        dateTmp = datetime.date(matDateTmp_yy, matDateTmp_mm, matDateTmp_dd)
        #dateTmp = dtime.date(matDateTmp_yy, matDateTmp_mm, matDateTmp_dd)
        
        if (i == 1):date0 = dateTmp

        ggTmp = (dateTmp - date0).days
        timeTmp = ggTmp/365.2425
        
        matDVec.append(dateTmp)
        matTVec.append(timeTmp)
        valVec.append(valTmp)

    """ 
    data_out['MatDate']  = matDVec[1:]
    data_out['MatTimes'] = matTVec[1:]
    data_out['ValoreZC'] = valVec[1:]
    """
    
    data_out['MatDate']  = matDVec
    data_out['MatTimes'] = matTVec
    data_out['ValoreZC'] = valVec


    return data_out



if __name__ == "__main__":
    

    #inputCDSData = r"input_test\cds_data_v0.txt"
    #inputBenchFile = r"input_test\raw_py_bench_curve.txt"
    #inputSwapFile = r"input_test\raw_py_swap_curve.txt"

    #inputCDSData = r"input_test\cds_data_flat_v0.txt"
    #inputCDSData = r"input_test\cds_data_v1.txt"
    #-----------------------------------------------------------
    

    #inputCDSData = r"input_test\boot_cds\test_2\cds_data_v2.txt"
    #inputCDSData = r"input_test\boot_cds\test_2\cds_data_mediobanca_27_11_17.txt"
    #inputCDSData = r"input_test\boot_cds\test_2\cds_IGT_USD_27_11.txt"
    inputCDSData = r"input_test\boot_cds\test_2\cds_IGT_EUR_27_11.txt"

    #inputSwapFile = r"input_test\boot_cds\test_2\raw_py_swap_curve_v3.txt"
    #inputBenchFile = r"input_test\boot_cds\test_2\raw_py_bench_curve_v3.txt"

    #inputSwapFile = r"input_test\boot_cds\test_2\raw_py_swap_curve_USD_27_11.txt"
    #inputBenchFile = r"input_test\boot_cds\test_2\raw_py_bench_curve_USD_27_11.txt"

    inputSwapFile = r"input_test\boot_cds\test_2\raw_py_swap_curve_EUR_27_11.txt"
    inputBenchFile = r"input_test\boot_cds\test_2\raw_py_bench_curve_EUR_27_11.txt"

    #fileOutRes = 'output_test/cds_data/test_2/cds_data_mediobanca_eur_27_11.txt'
    #fileOutRes = 'output_test/cds_data/test_2/cds_data_IGT_USD_27_11.txt'
    fileOutRes = 'output_test/cds_data/test_2/cds_data_IGT_EUR_27_11.txt'

    # --------------------- LOAD --------------------------------    
    import funzioni_boot_cds as f_cds

    #opt_dict       = f_cds.set_data_opt_for_cds(data_ref)   
    data_raw_cds   = load_cds_data_fromFile(inputCDSData)
    data_raw_swp   = load_simple_curves_fromFile(inputSwapFile)
    data_raw_bench = load_simple_curves_fromFile(inputBenchFile)
    
    dataRef = data_raw_swp['MatDate'][0]
    
    opt_dict       = f_cds.set_data_opt_for_cds(dataRef)   


    data_raw_bench['Type'] = 'swap'
    
    
    data_raw_zc =   inputSwapFile 
    data_raw = load_zc_data_fromFile(data_raw_zc)
    
    
    model_fit = '0' #0: LIN, 1: AVD, 2: SVE, 3: CIR, 4: NS 
    flag_make_graph = False    
    opt_dict_f = fb.set_data_opt_for_fit(model_fit, flag_make_graph, dataRef)


    c_dates = data_raw['MatDate']
    c_rates = data_raw['ValoreZC']
    c_rates = np.array(c_rates)/100.0
    
    dict_model_par = fb.fitting(c_dates, c_rates, opt_dict_f)
    
    #print 'dict_model_par: ', dict_model_par
    #fb.FQ(999)
    
    
    #-----------------------------------------


    data_raw_bench = set_model_prms(data_raw_bench)
    
    data_raw_bench['Model'] = 0    
    data_raw_bench['prms'] = dict_model_par
        
    #print 'len(prms[a]):', len(data_raw_bench['prms']['a'])


    print data_raw_bench.keys()
    
    
    
    data_raw_bench = convert2Df(data_raw_bench, 100.0)
    data_raw_swp   = convert2Df(data_raw_swp, 100.0)


    """
    print 'data_raw_swp.keys(): ', data_raw_swp.keys()
    print 'data_raw_bench.keys(): ', data_raw_bench.keys()
    print 'data_raw_cds.keys(): ', data_raw_cds.keys()
    print 'data_raw_cds[matTimes]: ', data_raw_cds['MatTimes']
    """

    
    print '----------------------------------------------------------------'
    print 'data_raw_swp[Basis]: ', data_raw_swp['Basis']
    print 'data_raw_swp[DiscountFactor]: ', data_raw_swp['DiscountFactor']
    
     
    print 'opt_dict: ', opt_dict
    print '-------------------------------------------'
    
    print 'data_raw_swp_DF: ', data_raw_bench.keys()

    data_raw_bench['ValoreNodo'] = []
    data_raw_swp['ValoreNodo'] = []

    
    print 'bench DF: ', data_raw_bench['DiscountFactor']
    print 'bench node: ', data_raw_bench['ValoreNodo']
    print 'bench model: ', data_raw_bench['Model']
    print 'bench prms: ', data_raw_bench['prms']['a']


    print '===================================================='
    print 'bench DF: ', data_raw_swp['DiscountFactor']
    print 'bench node: ', data_raw_swp['ValoreNodo']
    
    #print 'data_raw_cds: ', data_raw_cds

    
     
    #fb.FQ(999)
    
    

    #res = fb.boot_cds(opt_dict, data_raw_cds, data_raw_bench, data_raw_swp)
    res = f_cds.boot_cds(opt_dict, data_raw_cds, data_raw_bench, data_raw_swp)
    

    dumpResultsOnFile(res, fileOutRes) 

    """    
    print 'res[pySpread]: ', res['pySpread']
    print '----------------------------------------'
    print 'res[hazardRate]: ', res['hazardRate']
    print '----------------------------------------'
    print 'res[marginalDefault]: ', res['marginalDefault']
    print '----------------------------------------'
    print 'res[zcSpread]: ', res['zcSpread']
    print '----------------------------------------'
    print 'res[outputDates]: ', res['outputDates']


    """



