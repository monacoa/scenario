


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
                
        #dateTmp = dtime.date(matDateTmp_yy, matDateTmp_mm, matDateTmp_dd)
        dateTmp = datetime.date(matDateTmp_yy, matDateTmp_mm, matDateTmp_dd)

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


"""
def set_data_opt_for_cds():

    # bench_boot_method 
    # interp_method 
    
    opt_dict = {}

    opt_dict['DataRef']        = datetime.date(2005, 12, 31)

    opt_dict['MKT']            = 'de'
    opt_dict['tenor']          = 3 # 
    opt_dict['Basis']          = 'ACT/365' #ACT/365
    opt_dict['interp']         = '2' # 
    opt_dict['BusConv']        = 'modfollow'
    opt_dict['fixingDays']     = 2
    opt_dict['compounding']    = 0   #0 = semplice, 1 = composto, 2 = continuo
    opt_dict['ReocveryRate']   = 0.4
    opt_dict['hr_bootMethod']  = 1 #0 = LCS, 1 = CHR
    opt_dict['Basis']          = 'ACT/360' 
    opt_dict['BusConv']        = 'follow'


    opt_dict['opt_path_graph']  =  'C:\\'
    

    return opt_dict
"""

def set_opt_for_fit(fit_model, data_ref):
    
    
    opt_dict = {}
    
    opt_dict['interp'] = fit_model
    opt_dict['DataRef'] = data_ref
    
    opt_dict['MakeGraph'] = False
    
    

    opt_dict['opt_fwd_tenor']   = "1M"
    
    
    opt_dict['opt_path_graph']  =  'C:\\'
    
    #opt_dict['fit_type']        = "boot"
    opt_dict['fit_type']        = "py"
    

    
    if (opt_dict['interp'] == '2'): # 
    
        bound_min_sve = [0.001,  0.001, -0.05, -0.05, -0.05, -0.05]
        bound_max_sve = [100.0,  100.0,  10.0,  10.0,  10.0,  10.0]
        bound_x0_sve  = [  1.0,    1.0,  0.03,  0.03,  0.03,  0.03]

        opt_dict['bound_min_sve'] = bound_min_sve
        opt_dict['bound_max_sve'] = bound_max_sve
        opt_dict['x0']        = bound_x0_sve

    elif (opt_dict['interp'] == '3'): 

        bound_min_cir = [0.00001,   0.00001,   0.00001,   0.00001]
        bound_max_cir = [10.0, 10.0,  10.0, 0.05]
        bound_x0_cir  = [0.01, 0.03, 0.03, 0.1]

        opt_dict['bound_min_cir'] = bound_min_cir
        opt_dict['bound_max_cir'] = bound_max_cir
        opt_dict['x0']        = bound_x0_cir

    else:
        
        pass

    return opt_dict

def computeDateOutputCDS(raw_data, data_opt):
    
    ln = len(raw_data['MatTimes'])
    time_w_y   = raw_data['MatTimes'][ln-1]
    n_steps = int(time_w_y)
    ln = 31
    
    date_0 = data_opt['DataRef']
    date_0 = dtime.fromordinal(date_0.toordinal())
    
    dateOut = []
    timeOut = []
    for i in range(0, ln):
        
        timeTmp = i
        #timeTmp = int(raw_data['MatTimes'][i])
        
        dateTmp = date_0 + relativedelta(months=12*i)
        #dateTmp = date_0 + relativedelta(months=12*timeTmp)
        
        dd = int(dateTmp.day)
        mm = int(dateTmp.month)
        yy = int(dateTmp.year)
        
        dateDateTmp = datetime.date(yy, mm, dd)
        
        dateOut.append(dateDateTmp)
        timeOut.append(float(timeTmp))
        
    return  dateOut, timeOut
"""
"""
def computeTimesFromDates(dateList):

    matTimes = []
    dateRef  = dateList[0]
    
    ln = len(dateList)
    
    matTimes.append(0.0)
    
    for i in range(1, ln):
        
        timeTmp = (dateList[i] - dateRef).days/365.2425

        matTimes.append(timeTmp)

    return matTimes
"""
"""
def computePYRates(time_ref, df_ref, pyFreq, time_0, matTime):

    

    if (type(matTime) == list):

        py_out = []
        ln = len(matTime)
        for i in range(0, ln):
            
            
            matTimeTmp = matTime[i]
            if (matTimeTmp < pyFreq): matTimeTmp = pyFreq
            
            
            py_outTmp  = computePYRate(time_ref, df_ref, pyFreq, time_0, matTimeTmp)
            py_out.append(py_outTmp)


    else:    
        py_outTmp  = computePYRate(time_ref, df_ref, pyFreq, time_0, matTime)
        py_out = py_outTmp

    py_out = np.asanyarray(py_out)
    
    return py_out
"""
"""
def computePYRate(refTimes, refDiscounts, pyFreq, time_0, matTime):
    
    
    n_payments = int(matTime/pyFreq)
    
    
    refRates = np.zeros(len(refTimes))
    refRates[1:] = -np.log(refDiscounts[1:])/refTimes[1:]
    refRates[0] = refRates[1]
    
    
    r_0   = np.interp(time_0, refTimes, refRates)
    r_end = np.interp(matTime, refTimes, refRates)

    #z_0 = np.exp(-r_0*time_0)
    #z_end = np.exp(-r_end*matTime)
    
    n_0 = int(time_0/pyFreq)
    z_0   = 1.0/(1.0 + pyFreq*r_0)**n_0
    z_end = 1.0/(1.0 + pyFreq*r_end)**n_payments
    
    Ann = 0
    for i in range(1, n_payments+1):
        
        dtTmp = pyFreq
        t_i = dtTmp*i
        r_i = np.interp(t_i, refTimes, refRates)
        z_i = 1.0/(1.0 + dtTmp*r_i)**i
        #z_i = np.exp(-r_i*t_i)
        
        annTmp = float(z_i*dtTmp)
        Ann =  Ann + annTmp

    if (n_payments == 1):
        py_out = r_0

    else:
        py_out = (z_0 - z_end)/Ann
    
    return py_out


def bootstrap_lcs(cds_times, cds_values, cds_pay_times, RR, df_risk_free):
    # la funzione riceve in input i fattori di sconto risk-free gia interpolati#

    n_cf      = len(cds_pay_times)

    
    cds_val_n = np.interp(cds_pay_times, cds_times, cds_values)
    
    dt_pay = cds_pay_times[1:] - cds_pay_times[:n_cf-1]
    dt_pay  = np.insert(dt_pay, 0, dt_pay[0])
    
    surv = np.zeros(n_cf)
    hr   = np.zeros(n_cf)
    
    hr[0]   = 0.0
    surv[0] = 1.0

    hr[1]   = (1.0/dt_pay[0])*np.log(1.0 + (cds_val_n[1]*dt_pay[1])/(1.0 - RR))
    hr[0]   = hr[1]
    
    surv[1] = surv[0]*(np.exp(-hr[1]*(dt_pay[1])))


    #print ''
    for i in range(2, n_cf):
    #for i in range(2, 20):
        
        den1 = 0.0
        #print 'RR:', RR

        for k in range(1, i-1):
            #den1 = den1 + (1.0 - RR)*(1.0 - np.exp(-hr[k]*dt_pay[k]) - cds_val_n[i]*dt_pay[k]*np.exp(-hr[k]*dt_pay[k-1]))*df_risk_free[k]*surv[k - 1]
            den1 = den1 + (1.0 - RR)*(1.0 - np.exp(-hr[k]*dt_pay[k]) - cds_val_n[i]*dt_pay[k]*np.exp(-hr[k]*dt_pay[k-1]))*df_risk_free[k]*surv[k-1]

            #if (i == 3):
    
            #    print '**********************'
            #    print 'hr[k]: ', hr[k]
            #    print 'dt_pay[k]: ', dt_pay[k]
            #    print 'df_risk_free[k]: ', df_risk_free[k]
            #    print 'surv[k - 1]: ', surv[k - 1]
    
            #    print 'cds_val_n: ', cds_val_n[i]
                
        
        #for k=2:i-1
        #    den1 = den1 + (1.0 - RR)*(1.0 - np.exp(-h[k]*dt_pay[k]) - cds_interp(i)*delta_t(k)*np.exp(-h(k)*dt_pay[k-1]))*df(k)*s(k-1)

            #den1 = den1 + (1.0 - RR)*(1.0 - np.exp(-hr[k]*(td_pay[k])) - cds_interp(i)*delta_t(k)*exp(-h(k)*(times(k)-times(k-1)))) * discount_factors(k) * s(k-1);
                        
        
        #print 'df_risk_free[i]: ', df_risk_free[i]
        
        num  = (cds_val_n[i]*dt_pay[i] + (1 - RR)*surv[i - 1])*df_risk_free[i]
        den2 = (1.0 - RR)*surv[i - 1]*df_risk_free[i]

        
        den  = den1 + den2
        u_dt = 1.0/dt_pay[i]
        
        hr[i] = u_dt*np.log(num/den)

        surv[i] = surv[i - 1]*np.exp(-hr[i]*dt_pay[i])

        
        #print 'cds_val_n[i]: ', cds_val_n[i]
        #print 'dt_pay[i]: ', dt_pay[i]
        #print 'num: ', num
        #print 'den: ', den
        #print 'den1: ', den1
        #print 'den2: ', den2
        


    return surv, hr
    





def fromTime2Node(val_in):

    val_out = val_in

    return val_out

def buildDataRawForBoot(date, tempi, valori_tassi):


    data_raw_n = {}
    
    data_raw_n['TipoSegmento'] = {}
    data_raw_n['ValoreNodo'] = {}
    data_raw_n['UsaNodo'] = {}
    data_raw_n['MatDate'] = {}
    data_raw_n['Nodo'] = {}
    
    ln_tassi = len(date)
    for i in range(0, ln_tassi):
        
        dateTmp = date[i]
        timeTmp = tempi[i]  
        tassiTmp = valori_tassi[i]
        nodoTmp = fromTime2Node(timeTmp)

        data_raw_n['UsaNodo'][i] = 'Y'
        
        if (timeTmp <= 1):
            
            segTmp = 'D'
        else:
            segTmp = 'S'
                
        data_raw_n['TipoSegmento'][i] = segTmp
        data_raw_n['ValoreNodo'][i] = tassiTmp
        data_raw_n['MatDate'][i] = dateTmp
        data_raw_n['Nodo'][i] = nodoTmp

    return data_raw_n 


def fromDates2Times(date_list):

    target_times = []
    
    n_dates = len(date_list)
    
    date_0 = date_list[0]
    
    for i in range(0, n_dates):

    
        dateTmp      = date_list[i]
        diffdaysTmp  = (dateTmp - date_0).days
        diffTimesTmp = diffdaysTmp/365.2425
        
        target_times.append(diffTimesTmp)
        
    return target_times


def fromTimes2Dates(refDates, times_list):

    target_dates = []
    
    n_times = len(times_list)
    
    
    for i in range(0, n_times):
    
        timeTmp = times_list[i]
        daysTmp = int(365.2425*timeTmp)
        
        target_datesTmp = refDates  + datetime.timedelta(days=daysTmp)
        target_dates.append(target_datesTmp)
        
    return target_dates








def dumpDataOnFile(cds_times_n, surv_prob, fileOutSurv):
    
        
    n_times = len(cds_times_n)
    fileO = open(fileOutSurv, 'w')

    fileO.write("Data\t Valore\n")
    
    for i in range(0, n_times):

        dateTmp  =  cds_times_n[i]
        valueTmp =  surv_prob[i]   

        fileO.write("%s\t %s\n" %(dateTmp, valueTmp))
             

    fileO.close()
 
    
    
 



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

def showPrms(dict_model_par, dump_out):
    
    modelPrmsList = dict_model_par.keys()
    n_prms = len(modelPrmsList)
    
    dateList = dict_model_par['Dates']
    ln = len(dateList)
    
    if (n_prms == 6):
        fit_model = 'AVD'
    elif(n_prms == 7):
        fit_model = 'SVE'
    elif(n_prms == 5):
        
        try:
            dict_model_par['const1']
            fit_model = 'NS'
        except:
            fit_model = 'CIR'
            
    elif(n_prms == 3):
        fit_model = 'LIN'
    else:
        fit_model = '??'

        
    print 'fit_model: ', fit_model
    

    fout = open(dump_out, "w")
    
    fout.write('Dates \t')
    
    for i in range(0, len(modelPrmsList)):
        
        prmsNameTmp = modelPrmsList[i]

        if (prmsNameTmp == 'Dates'):
            pass
        else:
            fout.write(str(prmsNameTmp))
            fout.write("\t")
        
    fout.write("\n")

    
    for i in range(0, ln):
        
        
        dateTmp = dateList[i]
        #print 'n_prms: ', n_prms

        fout.write(str(dateTmp))
        fout.write("\t")
        
        for j in range(0, n_prms):
            
            prmsNameTmp = modelPrmsList[j]
    
            if (prmsNameTmp == 'Dates'):
                pass
            else:
                prmsValTmp  = dict_model_par[prmsNameTmp][i]
                fout.write(str(prmsValTmp[0]))
                fout.write("\t")

                print '%s,  %s: %s' %(dateTmp, prmsNameTmp, prmsValTmp)

        fout.write("\n")
            
        print '-------------------------------------------'
    
    fout.close()
    print 'fit_model: ', fit_model
            
            
            
        
        
        
    
    

def set_prms_manually(opt_dict, x0_prms, min_prms, max_prms):

    """    
    map_model_fit = {}
    map_model_fit['NS'] = '4'
    
    modelRef = map_model_fit[model_fit]
    """
    
    if (opt_dict['interp'] == '4'):
        
        opt_dict['bound_max_ns'] = max_prms
        opt_dict['bound_min_ns'] = min_prms
        opt_dict['x0'] = x0_prms

    elif (opt_dict['interp'] == '3'):
        
        opt_dict['bound_max_cir'] = max_prms
        opt_dict['bound_min_cir'] = min_prms
        opt_dict['x0'] = x0_prms


    return  opt_dict 


if __name__ == "__main__":
    


    """
    # ---------------------- NS PRMS ------------------------------
    x0_prms = [5.0003037,    0.0207240,    -0.0250637,    -0.0030637]
    min_prms = [5.0003037,    0.0207240,    -0.0250637,    -0.0030637]
    max_prms = [5.0003037,   0.0207240,    -0.0250637,    -0.0030637]
    """
    
    
    # ------------ CHK FITTING -------------------------------------

    
    #inputCurveFile = r"input_test\raw_short_swap_curves.txt"
    #inputCurveFile = r"input_test\zc_curve_to_fit_cs_v2.txt"
    inputCurveFile = r"input_test\fitting_data\test_swap_curve_12_07_17_D.txt"    
    inputCurveFile = r"input_test\fitting_data\test_swap_curve_12_07_17_x100_D.txt"    

    data_raw = load_zc_data_fromFile(inputCurveFile)
    
    dataRef = data_raw['MatDate'][0]
    
    model_fit = '0' #0: LIN, 1: AVD, 2: SVE, 3: CIR, 4: NS 
    flag_make_graph = True
    
    opt_dict = fb.set_data_opt_for_fit(model_fit, flag_make_graph, dataRef)

    """
    x0_prms = [5.0003037,    0.0207240,    -0.0250637,    -0.0030637]
    min_prms = [5.0003037,    0.0207240,    -0.0250637,    -0.0030637]
    max_prms = [5.0003037,   0.0207240,    -0.0250637,    -0.0030637]
    """

    #x0_prms = [0.000000,    0.130126,    0.027817,    0.241468]
    #min_prms = [0.000000,    0.130126,    0.027817,    0.241468]
    #max_prms = [0.000000,    0.130126,    0.027817,    0.241468]

    
    #opt_dict = set_prms_manually(opt_dict, x0_prms, min_prms, max_prms)


    c_dates = data_raw['MatDate']
    c_rates = data_raw['ValoreZC']
    c_rates = np.array(c_rates)/100.0
    
    
    
    
    dict_model_par = fb.fitting(c_dates, c_rates, opt_dict)
    dump_out = 'test_out_prms.txt'
    
    showPrms(dict_model_par, dump_out)
    print dict_model_par
    
    
    


