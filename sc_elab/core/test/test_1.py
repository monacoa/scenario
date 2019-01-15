


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

def set_data_opt_for_fit():
    
    """
    {'fit_type': 'boot', 
     'bound_max_cir': [0.02, 10.0, 1.0, 0.5], 
     'bound_max_sve': [100.0, 100.0, 10.0, 10.0, 10.0, 10.0], 
     'bound_min_cir': [0.0, 0.001, 0.005, 0.001], 
     'interp': '0', 
     'cir_params': ['r0', 'kappa', 'theta', 'sigma'], 
     'opt_path_graph': 'C://', 
     'sve_params': ['const1', 'const2', 'beta0', 'beta1', 'beta2', 'beta3'], 
     'bound_min_sve': [0.001, 0.001, -0.05, -0.05, -0.05, -0.05], 
     'opt_fwd_tenor': '1'}    

    """
    
    opt_dict = {}
    
    opt_dict['MakeGraph'] = False
    
    opt_dict['interp'] = '2' # 
    
    """
    opt_dict['interp'] == '0' # 'LIN'
    opt_dict['interp'] == '1' # 'AVD'
    opt_dict['interp'] == '2' # 'SVE'
    opt_dict['interp'] == '3' # 'CIR'
    """

    opt_dict['opt_fwd_tenor']   = "1M"
    
    """ 
    opt_dict['opt_fwd_tenor']   = "3M"
    opt_dict['opt_fwd_tenor']   = "6M"
    opt_dict['opt_fwd_tenor']   = "12M"
    """
    
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

"""


def set_data_opt():
    
    
    
    data_opt = {}
    data_opt['Basis'] = {}
    data_opt['Basis']['D']  = 'ACT/360'
    data_opt['Basis']['L']  = 'ACT/360'
    data_opt['Basis']['S']  = '30/360'
    data_opt['Basis']['F']  = 'ACT/360'

    data_opt['Basis']['TN'] = 'ACT/360'
    data_opt['Basis']['O/N'] = 'ACT/360'
    
    data_opt['BusConv'] = {}
    data_opt['ParConvexity'] = {}
    data_opt['ParConvexity']['A'] = 0.000005
    data_opt['ParConvexity']['B'] = 0.000005

    data_opt['FutureTenor'] = 90

    data_opt['TenorSwap'] = '6M'
    data_opt['Convexity'] = str(1) # 1 = attivo, 0 = inattivo
    data_opt['GapFutures'] = str(1) # 1 = previous spot rate, 0 = next forward rate
    data_opt['SwapGapMethod'] = str(0) # 0 = Constant Forward Rate, 1 = Linear Swap Rate
    data_opt['InterpLinFutures'] = str(1) # 1 = interp. esponenziale sui fattori di sconto, 0 = interp. line. sui tassi forward

    data_opt['RefDate'] = datetime.date(2017, 02, 15)
    #data_opt['RefDate'] = datetime.date(2016, 10, 02)

    data_opt['BusConv']['D'] = 'follow'
    data_opt['BusConv']['L'] = 'follow'
    data_opt['BusConv']['S'] = 'follow'
    data_opt['BusConv']['F'] = 'follow'
    data_opt['BusConv']['TN'] = 'follow'
    data_opt['BusConv']['O/N'] = 'follow'
    
    data_opt['MKT']          = 'us'

    data_opt['MakeGraph'] = True
    data_opt['SaveGraph'] = True
    data_opt['RegimeOutput'] = str(2) # 0 = Semplice, 1 = Composto, 2 = Continuo


    return data_opt

"""
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

#def bootstrap_lcs(cds_val_n, py_bench_rates):
"""
 
"""    
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
    
"""




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



"""
def fromCurveToSpread(df_bench_values, zc_bench_dates, prms_bench, bench_model, py_risky_val, py_risky_dates, risky_model, targetDates, targetTimes):
    
    
    dataRef = zc_bench_dates[0]
    fitting_opt_dict = set_opt_for_fit(risky_model, dataRef)    
    zc_bench_times   = fromDates2Times(zc_bench_dates)
    py_risky_times   = fromDates2Times(py_risky_dates)

    refDates = zc_bench_dates[0]
    
    prms_risky_for_py      = fb.fitting(py_risky_dates, py_risky_val, fitting_opt_dict)
    py_risky_val_n         = fb.makeRatesFromModel(py_risky_times, py_risky_val, prms_risky_for_py, targetTimes, risky_model)
    
    
    zc_risk_free           = fb.makeRatesFromModel(zc_bench_times, df_bench_values, prms_bench, targetTimes, bench_model)
    
    pyFreq = 0.25
    py_risk_free = computePYRates(zc_bench_times, df_bench_values, pyFreq, 0.0, targetTimes)
    
    py_spread = py_risky_val_n - py_risk_free
    
    
    ln_spread = len(py_spread)
    
    py_spread_n = []
    
    for i in range(0, ln_spread):
        
        py_spreadTmp = py_spread[i]
        
        if (py_spreadTmp <= 0.0):
            py_spreadTmp = 0.0
        else:
            py_spreadTmp = py_spreadTmp
            
        py_spread_n.append(py_spreadTmp)

    py_spread_n = np.asarray(py_spread_n)
        
    py_risky_val_n = py_risk_free + py_spread_n
    
    data_opt_for_boot = fb.set_data_opt()
    data_opt_for_boot['RefDate'] = dataRef
    data_opt_for_boot['SwapGapMethod'] = 1
    
    #n_rates = len(py_risky_val_n)
    
    target_dates = fromTimes2Dates(refDates, targetTimes)

    #data_raw_risky = buildDataRawForBoot(target_dates, targetTimes, py_risky_val_n)
    
    
    #--------------- IMPLEMENTA LA PARTE DI BOOT ----------------------
    
    zc_risky_val_n = py_risky_val_n
        
    zc_spread_n = zc_risky_val_n - zc_risk_free
        
        
    return py_spread_n, zc_spread_n
"""

"""
def func_chr(hr,b,cds_rr,s_n,t,z,cds_m,delta_t):
    
    dummy1 = 0
    dummy2 = 0
    
    for i in range(0,len(t)-1):
        dummy1 = dummy1 + ( np.exp(-hr*(t[i]-t[0])) - np.exp(-hr*(t[i+1]-t[0])) )*z[i]  
        dummy2 = dummy2 + delta_t[i]*z[i]*np.exp(-hr*(t[i+1]-t[0]))
    
    f = b - (1-cds_rr)*s_n*dummy1 + s_n*cds_m*dummy2
    
    return f
"""

"""
def func_h2(hr,cds_rr,t,z,cds_spread,delta_t):

    dummy1 = z[1:]*np.exp(-hr*(t[0:-1]-t[0]))*(1-np.exp(-hr*(t[1:]-t[0:-1])))
    dummy2 = z[1:]*delta_t[1:]*np.exp(-hr*(t[1:]-t[0]))
    
    f = cds_spread - (1-cds_rr) * sum(dummy1)/sum(dummy2)
    
    return f
"""

"""
def bootstrap_chr(maturities_cds, adjusted_cds, recovery_rate, times, discount_factors):
                                            
    # numero di flussi di cassa
    n_cf = len(times) 
                                                    
    # Inizializzo variabili
    h = np.zeros(n_cf)
    #exit_flag = 5.0*np.ones(n_cf)
    s = np.zeros(n_cf)
    
    # Imposto i dati iniziali
    s[0] = 1.0
    delta_t = times[1:] - times[:n_cf-1]
    delta_t = np.insert(delta_t, 0, 0)
    
    #k = find(maturities_cds(1)>times,1,'last');
    k =  fb.find_indx_next(maturities_cds[0], times)

    # trovo l'hazard rate costante fino alla prima maturity
    #options = optimset('fzero');
    h_star = 0.007
        
    h[1] = optimize.newton(func_h2, h_star, args=(recovery_rate,times[:k+1],discount_factors[:k+1],adjusted_cds[0],delta_t[:k+1]))
    #[h(2),exit_flag(2)] = fzero(@func_h2,h_star,options,recovery_rate,times(1:k),discount_factors(1:k),adjusted_cds(1),delta_t(1:k))
    
    #h(2) = 1/times(2)*log(1+adjusted_cds(1)*delta_t(2)/(1-recovery_rate));
    h[0]= h[1]
    s[1] = s[0]*np.exp(-h[1]*(times[1]-times[0]))

    for i in range(2,k+1):
            h[i]=h[1]
            s[i] = s[i-1]*np.exp( -h[1]*(times[i]-times[i-1]))


    
    j = k
    for i in range(k+1, n_cf):
    #for i in range(k+1, n_cf+1):

        if (i!=j+1):
            continue
        n = i-1
        s_n = s[i-1]

        #m = find(maturities_cds>times(i),1,'first')
        m =  fb.find_indx_next(times[i], maturities_cds)

        if (m == None):
        #if (m.size == 0):
            h[i] = h[i-1]
            #exit_flag[i]=exit_flag[i-1]
            s[i] = s_n*(np.exp(-h[i]*(times[i]-times[n])))
            continue
        
        cds_m = adjusted_cds[m]
        
        #j = find(maturities_cds(m)>times,1,'last')
        
        if (maturities_cds[m] > times[len(times)-1]):
            j = len(times)-1
        else:
            j =  fb.find_indx_next(maturities_cds[m], times)
            #j = j - 1

        b1 = 0
        b2 = 0
        
        for k in range(1,n+1):
            b1 = b1 + delta_t[k]*discount_factors[k]*s[k]
            b2 = b2 + discount_factors[k]*(s[k-1]- s[k])

        b = cds_m*b1 - (1.0 - recovery_rate)*b2
        #options = optimset('fzero')
        
       
        
        #h[i] = optimize.newton(func_chr, h[i-1], args=(b,recovery_rate,s_n,times[n:j+1],discount_factors[n+1:j+1],cds_m,delta_t[n:j]))
        h[i] = optimize.newton(func_chr, h[i-1], args=(b,recovery_rate,s_n,times[n:j+1],discount_factors[n+1:j+1],cds_m,delta_t[n:j]))

        #[h(i),exit_flag(i)] = fzero(@func_chr,h(i-1),options,b,recovery_rate,s_n,times(n:j),discount_factors(n+1:j),cds_m,delta_t(n:j-1));
        
        s[i] = s_n*np.exp(-h[i]*(times[i]-times[i-1]))

        for k in range(n+1,j+1):
            h[k] = h[i]
            #exit_flag[k] = exit_flag[i]
            s[k] = s_n*np.exp( -h[i]*(times[k]-times[n]))
          
    
    
    
    t_last = times[len(times)-1]
    t_end  = 30.0
    dt     = 0.25

    n_to_end = int((t_end - t_last)/dt) 

    h_ref = h[len(times)-1]
    s_n   = s[len(times)-1]

    for k in range(1, n_to_end + 1):
        h = np.append(h, h_ref)
        tTmp = t_last + k*dt
        times = np.append(times, tTmp)
        
        sTmp = s_n*np.exp( -h_ref*(tTmp - t_last))
        s    = np.append(s, sTmp)



    surv_p = s
    hazard_rates = h

    
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    # function: func_h2
    # funzione di cui si vuole trovare una radice: cds_spread1 - a = 0
    
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    # function: func_chr
    # funzione di cui si vuole trovare una radice: b - a = 0
    

    return surv_p, hazard_rates, times
"""


def dumpDataOnFile(cds_times_n, surv_prob, fileOutSurv):
    
        
    n_times = len(cds_times_n)
    fileO = open(fileOutSurv, 'w')

    fileO.write("Data\t Valore\n")
    
    for i in range(0, n_times):

        dateTmp  =  cds_times_n[i]
        valueTmp =  surv_prob[i]   

        fileO.write("%s\t %s\n" %(dateTmp, valueTmp))
             

    fileO.close()
 
    
    

"""
def boot_cds(opt_dict, raw_data, bench_data, swap_data):
    

    output_data = {}
    
    #--> definisco struttura output: lista date
    #outputDates
    
    dateOut, timeOut = computeDateOutputCDS(raw_data, opt_dict)
    benchTimes       = computeTimesFromDates(bench_data['MatDate'])
    
    #--> interpolo la curva dei cds sulle scadenze di OUTPUT
    # CDS_NEW
    
    cds_times_n = timeOut
    cds_values  = raw_data['ValoreNodo']
    cds_times   = raw_data['MatTimes']
    
    cds_values  = np.asarray(cds_values)/10000.0
    
    cds_val_n = np.interp(cds_times_n, cds_times, cds_values)

    #--> costruisco i tassi PY della curva benchmark alle scadenze delle date di OUTPUT
    #PY_NEW
    
    
    benchDf      = bench_data['DiscountFactor']
    prmsBench    = bench_data['prms']
    zcBenchDates = bench_data['MatDate']
    
    model_risky  = bench_data['Model']
    model_rf     = bench_data['Model']
    
    RR = opt_dict['ReocveryRate']

    model_risky = str(model_risky)
    model_rf    = str(model_rf)
    
    refDate = zcBenchDates[0]

    
    time_0 = 0.0
    
    cds_pay_freq = opt_dict['tenor']/12.0
    pyFreq = cds_pay_freq
    
    
    #rr_tmp = -np.log(benchDf)/benchTimes
    
    py_bench_rates = computePYRates(benchTimes, benchDf, pyFreq, time_0, timeOut)
    
    benchDf = np.asarray(benchDf)
    benchTimes = np.asarray(benchTimes)

    benchZcRates = -1.0/benchTimes*np.log(benchDf)
    benchZcRates[0] = benchZcRates[1]
    
    tipo_curva_bench = bench_data['Type']
    
    if(tipo_curva_bench != 'swap'):

        swapTimes_risk_free = computeTimesFromDates(swap_data['MatDate'])
        swapDf_risk_free    = swap_data['DiscountFactor']

        py_swap_rates   = computePYRates(swapTimes_risk_free, swapDf_risk_free, pyFreq, time_0, timeOut)

        py_spread_adj   = py_swap_rates - py_bench_rates
        cds_adj         = cds_val_n + py_spread_adj
        cds_val_n       = cds_adj
        

    py_risky = cds_val_n + py_bench_rates
    
    py_risky_times = timeOut
    py_risky_dates = dateOut
    
    model_risky = '0'
    
    pySpread, zcSpread = fromCurveToSpread(benchDf, zcBenchDates, prmsBench, model_rf, py_risky, py_risky_dates, model_risky, dateOut, timeOut)
    
    

    
    # calcolo numero di flussi di cassa
    
    ln = len(cds_times)
    n_cf = int(cds_times[ln-1]*(1.0/cds_pay_freq))
    

    cds_pay_dates = []
    
    mkt_code = opt_dict['MKT']
    basis_cds = opt_dict['Basis']
    day_conv_cds = opt_dict['BusConv'] 

    #from datetime import datetime



    mkt_ref  = holidays.get_calendar(mkt_code)

    for i in range(0, n_cf):

        cds_pay_dateTmp = refDate + relativedelta(months=int(cds_pay_freq*12*i))      
        #cds_pay_dateTmp = refDate + datetime.timedelta(weeks=int(cds_pay_freq*12*4*i))
        cds_pay_dateTmp = busdayrule.rolldate(cds_pay_dateTmp, mkt_ref, day_conv_cds)

        cds_pay_dates.append(cds_pay_dateTmp)
        
    cds_pay_times = daycount.yearfractions(cds_pay_dates, basis_cds)
    cds_pay_times = np.asarray(cds_pay_times)

    boot_method   = opt_dict['hr_bootMethod']
    
    
    benchZcRates_n = np.interp(cds_pay_times, benchTimes, benchZcRates)
    benchDf_n = np.exp(-benchZcRates_n*cds_pay_times)
    df_risk_free = benchDf_n
    
    dt = 0.25
    times_n = []
    t_last = np.around(cds_pay_times[len(cds_pay_times)-1], decimals = 0)
    #n_t = int(t_last/dt)

    n_t = int(t_last/dt)    
    for i in range(0, n_t+1):
        
        tTmp = dt*i
        times_n.append(tTmp)
        
    
    log_z = -np.log(df_risk_free[1:])/cds_pay_times[1:]
    log_z_int    = np.interp(times_n, cds_pay_times[1:], log_z)
    df_risk_free = np.exp(-log_z_int*times_n)
    
    cds_pay_times = times_n
    cds_pay_times = np.array(times_n)


    boot_method == 1
    
    if (boot_method == 0):
                        
        surv_prob_all, hr_values_all = bootstrap_lcs(cds_times_n, cds_val_n, cds_pay_times, RR, df_risk_free)

        log_surv = np.log(surv_prob_all)

        hr_values   = np.interp(cds_times_n, cds_pay_times, hr_values_all)
        log_surv    = np.interp(cds_times_n, cds_pay_times, log_surv)
        
        surv_prob = np.exp(log_surv) 

        
        #fileOutSurv = 'dataSurv_lcs_v0.txt'
        #fileOutHR   = 'dataaH_lcs_v0.txt'
        
        #dumpDataOnFile(cds_times_n, surv_prob, fileOutSurv)
        #dumpDataOnFile(cds_times_n, hr_values, fileOutHR)
        
            
    elif (boot_method == 1):
        
        cds_pay_times = np.array(cds_pay_times)

        df_risk_free = np.exp(-cds_pay_times*0.01)

        ln = len(cds_pay_times)
        surv_prob_all, hr_values_all, times_n = bootstrap_chr(cds_times, cds_values, RR, cds_pay_times, df_risk_free)
        
        hr_values   = np.interp(cds_times_n, times_n, hr_values_all)
        surv_prob   =  np.interp(cds_times_n, times_n, surv_prob_all)

        
        #fileOutHR         = 'output_test/dataH_chr_v2.txt'
        #fileOutSurv       = 'output_test/dataSurv_chr_v2.txt'
        #fileOutZCSpread   = 'output_test/dataSpread_zc_v2.txt'
        #fileOutPYSpread   = 'output_test/dataSpread_py_v2.txt'
        #fileOutMDefault   = 'output_test/dataMDefault_py_v2.txt'


        #fileOutHR         = 'output_test/dataH_chr_flat_v0.txt'
        #fileOutSurv       = 'output_test/dataSurv_chr_flat_v0.txt'
        #fileOutZCSpread   = 'output_test/dataSpread_zc_flat_v0.txt'
        #fileOutPYSpread   = 'output_test/dataSpread_py_flat_v0.txt'
        #fileOutMDefault   = 'output_test/dataMDefault_py_flat_v0.txt'
        
        #dumpDataOnFile(cds_times_n, surv_prob, fileOutSurv)
        #dumpDataOnFile(cds_times_n, hr_values, fileOutHR)
        #dumpDataOnFile(cds_times_n, zcSpread, fileOutZCSpread)
        #dumpDataOnFile(cds_times_n, pySpread, fileOutPYSpread)
        
    else:
        
        return 'metodo non gestito!!!'


    # calcolo le probabilita' marginali di default
    
    
    ln = len(surv_prob)
    marg_d   = 1.0 - surv_prob[1:ln]/surv_prob[0:ln-1]
    marg_d_n = np.insert(marg_d, 0, 0)

    #dumpDataOnFile(cds_times_n, marg_d_n, fileOutMDefault)


    output_data['outputDates']     = dateOut
    output_data['zcSpread']        = zcSpread
    output_data['pySpread']        = pySpread
    output_data['hazardRate']      = hr_values
    output_data['survProbCum']     = surv_prob
    output_data['marginalDefault'] = marg_d_n
    
    return output_data     
    

    #------ set up Output---------------------

"""


def convert2Df(data_raw_ref):
    
    
    data_discount = []
    n_dates = len(data_raw_ref['MatDate'])
    
    for i in range(0, n_dates):
        
        
        dateTmp = data_raw_bench['MatDate'][i]
        valTmp = data_raw_bench['ValoreNodo'][i]
        
        if (i == 0):date0 = dateTmp
        
        ggTmp = (dateTmp - date0).days
        timeTmp = ggTmp/365.2425
        
        dfTmp = np.exp(-timeTmp*valTmp)
        
        
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

if __name__ == "__main__":
    """
    cc = CdsCurve()

    cc.ref_date       = datetime.date(day=int(curve_date[-2:]), month=int(curve_date[5:7]), year=int(curve_date[:4]))
    cc.description    = curve_des

    cc.ratingProvider = app.new_window.new_window.rating.get()
    cc.sectorProvider = app.new_window.new_window.sector.get()
    cc.loadDataFromDB()
    writeCurveOnXls(cc, nameSheet, xla, curve_type)
    """

    #inputCDSData = r"input_test\cds_data_v0.txt"
    #inputBenchFile = r"input_test\raw_py_bench_curve.txt"
    #inputSwapFile = r"input_test\raw_py_swap_curve.txt"

    #inputCDSData = r"input_test\cds_data_flat_v0.txt"
    #inputCDSData = r"input_test\cds_data_v1.txt"
    #-----------------------------------------------------------
    
    """
    inputCDSData = r"input_test\cds_data_flat_v0.txt"
    inputSwapFile = r"input_test\raw_py_swap_curve_v1.txt"
    inputBenchFile = r"input_test\raw_py_bench_curve_v1.txt"
    """

    
    #inputCDSData = r"input_test\cds_data_v2.txt"
    inputCDSData = r"input_test\cds_data_flat_v0.txt"
    inputSwapFile = r"input_test\raw_py_swap_curve_v1.txt"
    inputBenchFile = r"input_test\raw_py_bench_curve_v1.txt"
    
    
    opt_dict       = fb.set_data_opt_for_cds()   
    data_raw_cds   = load_cds_data_fromFile(inputCDSData)
    data_raw_swp   = load_simple_curves_fromFile(inputSwapFile)
    data_raw_bench = load_simple_curves_fromFile(inputBenchFile)
    
    
    data_raw_bench['Type'] = 'GOV'
    
    
    data_raw_bench = set_model_prms(data_raw_bench)
    
    """
    print 'data_raw_bench: ', data_raw_bench
    print 'data_raw_swp: ', data_raw_swp
    FQ(888)
    """
    
    data_raw_bench = convert2Df(data_raw_bench)
    data_raw_swp   = convert2Df(data_raw_swp)
    
    #print 'data_raw_bench: ', data_raw_bench
    #print 'data_raw_swp: ', data_raw_swp
    

    print 'data_raw_bench: ', data_raw_bench.keys()
    print 'data_raw_cds: ', data_raw_cds
    print 'data_raw_bench: ', data_raw_bench.keys()
    print 'data_raw_swp: ', data_raw_swp.keys()
    
    #print 'data_raw_swp[Basis]: ', data_raw_swp.keys()
    print 'data_raw_swp[ValoreNodo]: ', data_raw_swp['ValoreNodo']
    
    #'Basis', 'ValoreNodo', 'DiscountFactor', 'MatDate', 'prms', 'Model', 'Type'
    #'DiscountFactor', 'ValoreNodo', 'MatDate', 'Basis'

    
    FQ(111)
    
    
    res = fb.boot_cds(opt_dict, data_raw_cds, data_raw_bench, data_raw_swp)
    
    
    #print 'res: ', res.keys()

    """    
    print 'res[pySpread]: ', res['pySpread']
    print '----------------------------------------'
    """

    print 'res[survProbCum]: ', res['survProbCum']
    print '----------------------------------------'

    """
    print 'res[hazardRate]: ', res['hazardRate']
    print '----------------------------------------'

    print 'res[marginalDefault]: ', res['marginalDefault']
    print '----------------------------------------'

    print 'res[zcSpread]: ', res['zcSpread']
    print '----------------------------------------'

    print 'res[outputDates]: ', res['outputDates']
    """
    FQ(229)
    
