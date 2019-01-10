


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

from datetime import datetime as ddtime
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
        
                
        #dateTmp = ddtime.date(matDateTmp_yy, matDateTmp_mm, matDateTmp_dd)
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
        
        #dateTmp = ddtime.date(matDateTmp_yy, matDateTmp_mm, matDateTmp_dd)
        dateTmp = datetime.date(matDateTmp_yy, matDateTmp_mm, matDateTmp_dd)
        
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


def boot3s_elab_v2(data_opt, data_raw):

    #%--------------------------------------------------------------------------
    #%-------------- RETRIEVE DATA  -----------------------------------------
    #%--------------------------------------------------------------------------


    ref_field       = 'UsaNodo'
    val_not_allowed = ['N', 'n']

    data_raw_p = purge_data(data_raw, ref_field, val_not_allowed)
    
    dict1s, dict_f, dict_s = select_segments(data_raw_p)
    
    flag_f = (len(dict_f['Nodo'])>0)
    flag_s = (len(dict_s['Nodo'])>0)
    
    chkDataCoherence(dict1s, dict_f, dict_s)
    
    #FQ(2233)
    s1 = list(set(dict1s['TipoSegmento']))
    s2 = list(set(dict_f['TipoSegmento']))
    s3 = list(set(dict_s['TipoSegmento']))
    

    codeCurve = s1[0]
    if (len(s2)>0): codeCurve = codeCurve + s2[0] 
    if (len(s3)>0): codeCurve = codeCurve + s3[0] 

    seg1_dates  = dict1s['MatDate']
    seg1_values = dict1s['ValoreNodo']
    seg1_values = seg1_values/100.0

    n1s  = len(seg1_dates)

    if (flag_f == True):

        futures_start_dates = dict_f['MatDate']
        futures_values      = dict_f['ValoreNodo']
        nf   = len(futures_start_dates)

        basis_f = data_opt['Basis']['F']
        basis_f = convert_basis(basis_f)
        day_conv_f = data_opt['BusConv']['F']
    
    
    if (flag_s == True):
        
        swap_dates          = dict_s['MatDate']
        swap_val            = dict_s['ValoreNodo']

        swap_val = swap_val/100.0
        nsw = len(swap_dates)
        
        basis_s = data_opt['Basis']['S']
        basis_s = convert_basis(basis_s)
        day_conv_s = data_opt['BusConv']['S']

        par_a = data_opt['ParConvexity']['A']
        par_b = data_opt['ParConvexity']['B']
    
        par_convexity    = [par_a, par_b]
    
        tenor_swap       = data_opt['TenorSwap']
        tenor_swap       = convertNodeToMnth(tenor_swap)

    #%--------------------------------------------------------------------------
    #%-------------- RETRIEVE OPTIONS SETUP ------------------------------------
    #%--------------------------------------------------------------------------

    setting_default = set_data_default()

    ref_date = data_opt['RefDate'] 
    
    flag_convexity   = data_opt['Convexity']
    flag_futures_gap = data_opt['GapFutures']
    flag_method_swap = data_opt['SwapGapMethod']
    flag_interp1     = setting_default['InterpLinFutures']

    flag_make_graph = data_opt['MakeGraph']
    flag_save_graph = setting_default['SaveGraph']
    regime_output   = data_opt['RegimeOutput']   
    future_tenor    = data_opt['FutureTenor'] 
    
    '------- retrive conventions ------------------------'

    mkt_code = data_opt['MKT']
    mkt_ref  = holidays.get_calendar(mkt_code)

    day_conv_tn = setting_default['BusConv']['TN']
    day_conv_on = setting_default['BusConv']['O/N']

    tipo_s1        = dict1s['TipoSegmento'][0] 
    
    basis_s1       = data_opt['Basis'][tipo_s1]
    flag_regime_1s = setting_default['RegimeRate'][tipo_s1]
    
    basis_fix = setting_default['BasisFix']
    basis_ref = setting_default['BasisRef']
    
    '------------------ inizializzazione data output -----------------------'
    
    dtype  = [('Times', float), ('Dates', datetime.date), ('Df', float), ('Rates', float)]
    data_0 = np.array([], dtype=dtype)

    #%--------------------------------------------------------------------------
    #%-------------- START ELABORATION -----------------------------------------
    #%--------------------------------------------------------------------------

    '-------------- generazione scadenza futures ------------------------------'
    
    if (flag_f == True):


        futures_end_dates = increase_datetime_list(futures_start_dates, future_tenor)
        futures_end_dates = busdayrule.rolldates(futures_end_dates, mkt_ref, day_conv_f)
        futures_end_dates = conv_to_dates(futures_end_dates)
        
    '------ generazione dei tempi a partire dalle date  ----------------------------'

    if (len(seg1_values) > 0):
        
        seg1_dates_n = seg1_dates
        seg1_dates_n.insert(0, ref_date)
        
        seg1_times = daycount.yearfractions(seg1_dates_n, basis_s1)
        
        seg1_dates = seg1_dates[1:]
        seg1_times = seg1_times[1:]
        
        seg1_times = np.asarray(seg1_times)
        
        
    
    if (flag_f == True):

        futures_end_dates_n = futures_end_dates
        futures_end_dates_n.insert(0, ref_date)

        futures_start_dates_n = futures_start_dates
        futures_start_dates_n.insert(0, ref_date)

        futures_start_dates = futures_start_dates[1:]
        futures_end_dates   = futures_end_dates[1:]

        futures_end_times    = daycount.yearfractions(futures_end_dates_n, basis_f)
        futures_start_times  = daycount.yearfractions(futures_start_dates_n, basis_f)

        futures_end_times   = futures_end_times[1:]
        futures_start_times = futures_start_times[1:]
        
        futures_end_tx      = daycount.yearfractions(futures_end_dates_n, basis_ref)
        futures_start_tx    = daycount.yearfractions(futures_start_dates_n, basis_ref)
        
        futures_end_tx   = futures_end_tx[1:]
        futures_start_tx = futures_start_tx[1:]



    if (flag_s == True):

        
        swap_dates_n = swap_dates
        swap_dates_n.insert(0, ref_date)
        
        swap_times           = daycount.yearfractions(swap_dates_n, basis_s)
        swap_times_x         = daycount.yearfractions(swap_dates_n, basis_ref)

        swap_times           = np.asanyarray(swap_times)
        swap_times_x         = np.asanyarray(swap_times_x)

        swap_dates           = swap_dates[1:]

    #%-----------------------------------------------------------------------------------------'
    #%------ generazione fattori di sconto per dati OverNight, Tomorrow Next ------------------'
    #%-----------------------------------------------------------------------------------------'


    ref_date_next = ref_date + datetime.timedelta(days = 1)    
    on_date    = busdayrule.rolldate(ref_date_next, mkt_ref, day_conv_on)

    on_date_next = on_date + datetime.timedelta(days = 1)
    tn_date = busdayrule.rolldate(on_date_next, mkt_ref, day_conv_tn) 
    
    fix_date1      = tn_date
    fix_time1      = daycount.yearfrac(ref_date, fix_date1, basis_fix)

    on_date_n2     = on_date + datetime.timedelta(days = 2)    
    fix_swap_date  = busdayrule.rolldate(on_date_n2,  mkt_ref, day_conv_s) 
    fix_swap_time  = daycount.yearfrac(ref_date, fix_swap_date, basis_s)
    
    yy_fix = fix_swap_date.year
    mm_fix = fix_swap_date.month
    dd_fix = fix_swap_date.day
    
    fix_swap_date = datetime.date(yy_fix, mm_fix, dd_fix)

    
    #%--------------------------------------------------------------------------
    #%-------------- ELABORAZIONE 1mo SEGMENTO: DEPOSITI, LIBOR ----------------
    #%--------------------------------------------------------------------------


    seg1_df = np.zeros(n1s + 1)
    

    if (n1s > 0):


        #------ generazione fattori di sconto alle date ON/TN7FIX --------------------------------
            
        seg1_df = np.ndarray(n1s)
        
        seg1_df, df_fix_seg1, id_start = compute_first_df(seg1_df, seg1_dates, seg1_times, seg1_values, fix_time1, on_date, tn_date)
      
    
        seg1_times_m_fix = seg1_times[id_start:] - fix_time1
        seg1_times_m_fix = seg1_times[id_start:] - fix_time1
      
        seg1_val_m_fix = seg1_values[id_start:]
        
      
        if (flag_regime_1s == 0): # regime semplice
    
            seg1_df_ = df_fix_seg1*df_simple(seg1_val_m_fix, seg1_times_m_fix)
    
        elif (flag_regime_1s == 1): # regime composto
    
            seg1_df_ = df_fix_seg1*df_cmp(seg1_val_m_fix, seg1_times_m_fix)
    
        elif (flag_regime_1s == 2): # regime composto
    
            seg1_df_ = df_fix_seg1*df_cont(seg1_val_m_fix, seg1_times_m_fix)
        else:
            
            print 'Regime non gestito'
        
        
        seg1_df[id_start:] = seg1_df_

        '----------- set output -------------------------------------'

        seg1_r = compute_rates(seg1_df, seg1_times, regime_output)
        
        seg1_r     = np.insert(seg1_r, 0, seg1_r[0])
        seg1_df    = np.insert(seg1_df, 0, 1.0)
        seg1_times = np.insert(seg1_times, 0, 0.0)
        seg1_dates.insert(0, ref_date)
        
    
    data_merge = merge_data(data_0, seg1_dates, seg1_times, seg1_r, seg1_df)
    
    
    
    #%--------------------------------------------------------------------------
    #%-------------- ELABORAZIONE SEGMENTO FUTURES -----------------------------
    #%--------------------------------------------------------------------------

    if (flag_f == True):
  
        futures_rates = (100.0 - futures_values)/100.0  #-------- tasso future
    
        #-------- Correzione: "Convexity"------------------------------------------------'
        
        if (int(flag_convexity) == 1):
            conv_correction = compute_convexity(par_convexity, futures_start_tx, futures_end_tx)
            
            futures_rates = futures_rates - conv_correction

        #----------------------------------------------------------------------------------------
        #------------------ generazione di Z alla prima data "start" del futures -----------------
        #-----------------------------------------------------------------------------------------
    
        futures_end_df   = np.ndarray(nf)

        futures_df       = np.ndarray(2*nf)
        futures_times    = np.ndarray(2*nf)

        futures_dates    = []
    
        for i in range(0, nf):
            
            futures_start_df_tmp, futures_end_df_tmp, futures_end_df = compute_df_future(seg1_times, 
                                                                         seg1_values, 
                                                                         seg1_df,  
                                                                         futures_rates, 
                                                                         flag_futures_gap, 
                                                                         futures_start_times, 
                                                                         futures_end_times,  
                                                                         futures_end_df, 
                                                                         flag_interp1, 
                                                                         i)
            
            
            futures_df[2*i]     = futures_start_df_tmp
            futures_df[2*i + 1] = futures_end_df_tmp

            futures_times[2*i] = futures_start_times[i]
            futures_times[2*i + 1] = futures_end_times[i]
            
            futures_dates.append(futures_start_dates[i])
            futures_dates.append(futures_end_dates[i])
            
        
        futures_rates_boot = compute_rates(futures_df, futures_times, regime_output)
        

    
        df_array = futures_df
        r_array  = futures_rates_boot
        t_array  = futures_times
        d_array  = futures_dates

        '------------------- merge ------------------'

        data_merge = merge_data(data_merge, d_array, t_array, r_array, df_array)
    
    else:

        futures_end_times = [999]
        futures_start_times = [999]
        futures_times       = [999]
        futures_rates       = [999]
        futures_df          = [999]
        


    merge_dates = data_merge['Dates']
    merge_times = data_merge['Times']
    merge_df    = data_merge['Df']

    #print'merge_times: ', len(merge_times) 
    
    #swap_discount_factor = np.zeros(len(swap_dates))

    #%----------------------------------------------------------------------------------------------
    #%-------------- ELABORAZIONE SEGMENTO SWAP: FATTORI DI SCONTO A PARTIRE DAI TASSI SWAP ---------
    #%----------------------------------------------------------------------------------------------
    
    #'merge_df: ', merge_df
    #FQ(2233)
    
    if (flag_s == True):

        swap_discount_factor = np.zeros(len(swap_dates))

        index_last = find_indx(fix_swap_date, merge_dates)

        t_o =  merge_times[index_last]
        t_n =  merge_times[index_last+1]

        df_o = merge_df[index_last]
        df_n = merge_df[index_last + 1]
        
        t_target = fix_swap_time
        
        df_fix_swap = df_from_interp_df_exp(t_target, df_n, df_o, t_n, t_o)

    
        #%--------------------------- Def. variabili ----------------------------------------------------------
    
        dt_time_swap     = swap_times[1:]- swap_times[0:len(swap_times)-1]  # differenze tempi
        
        
        index_condition_irr  = dt_time_swap.round(0) > 1

        idx_start_irregular  = np.where(index_condition_irr)  # vettore indici date "irregolari" (array, dtype=int64)        
        idx_start_irregular  = idx_start_irregular[0]         # vettore indici date "irregolari"

        idx_end_irregular=[]
        #idx_start_irregular = index_irregular_s
        
        for i in range(0, len(idx_start_irregular)-1):
            
                irr_s = idx_start_irregular[i]
                irr_e = idx_start_irregular[i+1]
            
                if (irr_e != irr_s + 1):
        
                    idx_end_irregular.append(irr_s + 1)
                else: 
                    idx_end_irregular.append(irr_e)
                    
        if (len(idx_start_irregular) > 0):
            idx_end_irregular.append(len(dt_time_swap))
        
            idx_end_irregular = np.asarray(idx_end_irregular)
        
        start_irregular_time = []
        end_irregular_time = []
        end_irregular_time_int = []
        
        for i in range(0, len(idx_end_irregular)):
        
            n_i = idx_start_irregular[i]  #% indice inizio data irregolare
            m_i = idx_end_irregular[i]    #% indice fine data irregolare
            
            
            swap_times_n  = swap_times_x[n_i]  #% tempo inizio data irregolare
            swap_times_m  = swap_times_x[m_i]  #% tempo fine data irregolare
    
            t_swt_n = np.round(swap_times_n)   #% n. anni dell'n-ma data "irregolare"
            t_swt_m = np.round(swap_times_m)   #% m. anni dell'm-ma data "irregolare"

            start_irregular_time.append(t_swt_n)
            end_irregular_time.append(t_swt_m)
            end_irregular_time_int.append(int(t_swt_m))
        
        start_irregular_time = np.asarray(start_irregular_time)
        end_irregular_time = np.asarray(end_irregular_time)
        end_irregular_time_int = np.asarray(end_irregular_time_int)
        
    
        
        index_irregular  = np.where(index_condition_irr)  # vettore indici date "irregolari" (array, dtype=int64)        
        index_irregular  = index_irregular[0]         # vettore indici date "irregolari"

        
        ln_irr = len(index_irregular)

        if (ln_irr == 0):
            index_irregular = np.insert(index_irregular, 0, len(dt_time_swap))
        else:
            index_irregular = np.append(index_irregular, index_irregular[ln_irr-1] + 1) # viene aggiunto un indice in piu'

        n_new_swap       = int(round(swap_times_x[len(swap_times_x)-1]/(tenor_swap/12.0)))

        
        swap_dates_new   = []
        times_swap_new   = np.zeros(n_new_swap + 1)
        times_swap_x_new = np.zeros(n_new_swap + 1)
    
        for j in range(0, n_new_swap + 1):
        
            swap_dates_tmp = add_months(fix_swap_date, tenor_swap*j)            

            swap_dates_tmp = busdayrule.rolldate(swap_dates_tmp, mkt_ref, day_conv_s) #% genero vettore delle date corrette DAY CONVENTION
            swap_dates_new.append(swap_dates_tmp.date())
            
            times_swap_tmp = daycount.yearfrac(ref_date, swap_dates_tmp, basis_s)
            times_swap_new[j] = times_swap_tmp
            
            times_swap_x_tmp = daycount.yearfrac(fix_swap_date, swap_dates_tmp, basis_ref); # basis_ref al posto di basis_s
            times_swap_x_new[j] = times_swap_x_tmp
    
        # ---------------------------------------------------------------------

        tenors_0  = times_swap_new[0]
        tenors_1  = times_swap_new[1:] - times_swap_new[0:len(times_swap_new)-1]
    
        tenors    = np.insert(tenors_1, 0, tenors_0)
        
        swap_times  = swap_times[1:]
        #index_start = find_indx(swap_times[0], times_swap_new) # indice di inizio del segmento degli swap
        


        times_swap_       = np.round(times_swap_x_new/(tenor_swap/12.0))*(tenor_swap/12.0)
        
        sw_0              = np.round(swap_times[0]/(tenor_swap/12.0))*(tenor_swap/12.0)
        index_start       = find_indx_toll(sw_0, times_swap_, 0.001)
        #index_start       = find_indx(swap_times[0], times_swap_)

        times_swap_tmp1s  = times_swap_new[:index_start]   # tempi fino al primo tasso swap
        times_swap_tmp2s  = times_swap_new[index_start:]  # tempi oltre il primo tasso swap

        dates_swap_tmp1s  = swap_dates_new[:index_start]   # date fino al primo tasso swap        
        dates_swap_tmp2s  = swap_dates_new[index_start:]  # date oltre il primo tasso swap
        
        nsw1s             = len(times_swap_tmp1s)
        nsw2s             = len(times_swap_tmp2s)
    
        df_tmp1s          = np.zeros(nsw1s) # fattori di sconto fino al primo tasso swap

        dates_swap_2s_ordinal = from_date_to_ordinal(dates_swap_tmp2s)
        swap_dates_ordinal    = from_date_to_ordinal(swap_dates)
        

        #%-------------------------------------------------------------------------------------
        #%----------------------------- metodo lsr --------------------------------------------
        #%-------------------------------------------------------------------------------------
    
        #%----------------------- Interoplazione lineare dei tassi swap alle date di interesse ----------------------
    
        swap_val_new  = np.interp(dates_swap_2s_ordinal, swap_dates_ordinal, swap_val)#; % interpolazione


        #%------------------------------------------------------------------------------------------
        #%----------------------- Generazione fattori di sconto alle date dei pagamenti "Swap"
        #%--------------------------------------------------------------------------------------------
        
        #print'merge_times: ', merge_times[len(merge_times)-1]
        #print'merge_times: ', merge_times
       

        for i in range(0, nsw1s):
            
            if (merge_dates[len(merge_dates)-1] < dates_swap_tmp1s[i]):  #%--- NO INTERPOLAZIONE ---------------
    
                #%------------ gap colmato mantenendo costante il tasso fwd -----------
    
                if all(merge_df == 1):
                    fwd = swap_val[0]
                    df_tmp1s[i] = np.exp(-fwd*(times_swap_tmp1s[i]- merge_times[len(merge_times)-1]))
            
                else:

                    fwd         = -np.log(merge_df[len(merge_df)-1]/merge_df[len(merge_df)-2])/(merge_times[len(merge_times)-1] - merge_times[len(merge_times)-2])
                    df_tmp1s[i] = np.exp(-fwd*(times_swap_tmp1s[i]- merge_times[len(merge_times)-1]))
    
            else: #%--- INTERPOLAZIONE --------------------------
    
                overlap_flag = (dates_swap_tmp1s[i] in merge_dates)
    
                if (overlap_flag == False): # se il fattore di sconto non e' stato calcolato
    
                    #----------intepolazione esponenziale sui fattori di sconto ------------
    
                    index_1 = find_indx(dates_swap_tmp1s[i], merge_dates)

                    t_o     = merge_times[index_1]
                    df_o    = merge_df[index_1]

                    t_n     = merge_times[index_1 + 1]
                    df_n    = merge_df[index_1 + 1]
                    
                    t_target = times_swap_tmp1s[i]

                    df_tmp1s[i] = df_from_interp_df_exp(t_target, df_n, df_o, t_n, t_o)

                    

                else:
    
                    #%---------- viene selezionato il fattore di sconto gia' calcolato -------
    
                    index_0 = find_indx_equal(dates_swap_tmp1s[i], merge_dates)                    
                    df_tmp1s[i] = merge_df[index_0]
    
        
        
        t_z_all      = np.concatenate((times_swap_tmp1s, times_swap_tmp2s))
        z_all        = np.concatenate((df_tmp1s, np.zeros(nsw2s)))
    
        #%-------------------------------------------------------------------------------------
        #%--------------- Generazione dei fattori di sconto alle date dei tassi swap ----------
        #%-------------------------------------------------------------------------------------
        
        for i in range(0,nsw2s):
            
            z_all[nsw1s + i] = compute_z_from_swap_rate(z_all, df_fix_swap, tenors, swap_val_new[i], tenors[nsw1s + i])
            
        #%--------------------------------  set LSR method OUTPUT --------------------------------

        if (int(flag_method_swap) == 1): #% considero il metodo: Linear Swap Rate: 'LSR'

            swap_times_tmp = np.insert(times_swap_new, 0, 0)
            swap_times_tmp = np.round(swap_times_tmp, 2)
            z_out = z_all

            #%------------------------------------------------------------------------------------
            #%--------------------------------  CFR METHOD ----------------------------------------%
            #%------------------------------------------------------------------------------------
        
        
        elif (int(flag_method_swap) == 0) and (ln_irr > 1): # Considero il metodo "Constant fwd rate" 'CFR'
            
            index_irregular_first = int(np.round(swap_times_x[index_irregular[0]]/(tenor_swap/12.0))) # % indice prima data irregolare
            #index_irregular_last  = int(np.round(swap_times_x[index_irregular[len(index_irregular)-1]]/(tenor_swap/12.0))) #; % indice ultima data irregolare


            index_irregular_0     =  max(len(df_tmp1s)+1, index_irregular_first + 1)

            #z_out = np.zeros(max(1, index_irregular_last + 1))
            z_out = np.zeros(len(z_all))
                
            z_out[:index_irregular_0] = z_all[: index_irregular_0]

            
            times_swap_new  = times_swap_new + fix_swap_time

            k = 0
            
            #print 'end_irregular_time: ', end_irregular_time_int
            #print '----------------------------------------------------'
            #print 'nsw1s: ', nsw1s 
            #print 'nsw1s + nsw2s: ', nsw1s + nsw2s 
            
            for i in range(nsw1s, nsw2s+nsw1s):     #% ciclo su tutte le date irregolari
            
                t_ref = times_swap_[i]
                
                if ((t_ref)  in  end_irregular_time_int):
                
                    m_i = idx_end_irregular[k]    #% indice fine data irregolare
                    
                    t_swt_n = start_irregular_time[k]
                    t_swt_m = end_irregular_time[k]

                    swp_tmp =  swap_val[m_i - 1]
                    
                    z_out = compute_z_from_cfr(z_out, tenors, swp_tmp, times_swap_new, t_swt_n, t_swt_m, times_swap_, df_fix_swap)
    
                    k = k + 1
                else:

                    
                    t_swt_n = t_ref
                    
                    nz  = find_indx_n(t_swt_n, times_swap_) #% indice del vettore times_out corrispondente a n anni
                    
                    swp_tmp =  swap_val_new[i-nsw1s]
                    
                    z_out_n  = z_out[:nz]
                    tenors_n = tenors[:nz]

                    z_out[nz] = compute_z_from_swap_rate(z_out_n, df_fix_swap, tenors_n, swp_tmp, tenors[i])

                    
                    """
                    dt_ref = []
                    
                    t_n = times_swap_[i]
                    t_o = times_swap_[i-1]

                    dt_ref.append(t_o)  
                    dt_ref.append(t_n)  

                    dt_ref = np.asarray(dt_ref)
                    dz_ref = z_out[nz-1:nz+1]

                    r_out = compute_rates(dz_ref, dt_ref, regime_output)
                    
                    
                    print 'r_out: ', r_out[0] 
                    print 't_n: ', t_n
                    print '-----------------------------------'
                    """
        else:
            
            z_out = z_all

            #times_swap_       = np.round(times_swap_x_new/(tenor_swap/12.0))*(tenor_swap/12.0)            


            
        #%---------------------------------------------------------------------
        #%------------------ set swap output ----------------------------------
        #%---------------------------------------------------------------------
        
        times_swap_         = np.round(times_swap_x_new/(tenor_swap/12.0))*(tenor_swap/12.0)            
        times_swap_n        = times_swap_[1:]
        swap_times_adj_n    = np.round(swap_times_x[1:]/(tenor_swap/12.0))*(tenor_swap/12.0)
        
        
        #swp__ = np.round(swap_times_x/(tenor_swap/12.0), 0)
        #swap_times_adj = np.round(swp__*(tenor_swap/12.0),1)
        
        #if (swap_times_adj[0] == 0): swap_times_adj = swap_times_adj[1:]

        #print '--------------------------------------------'
        #print '--------------------------------------------'
        
        #print 'len(z_out): ', len(z_out)
        #print 'len(swap_discount_factor): ', len(swap_discount_factor)
        #print 'index_irregular_last N: ', index_irregular_last

        #print 'times_swap_: ',times_swap_n
        #print 'z_out: ', z_out
                    
        
        #print 'nsw: ', nsw
        
        for i in range(0, nsw):
            
            #indxTmp = find_indx_n(swap_times_adj[i], times_swap_[1:])
            #indxTmp = find_indx_n(swap_times_adj[i], times_swap_[1:])
            indxTmp = find_indx_n(swap_times_adj_n[i], times_swap_n)
            
            #print 'swap_times_adj: ',swap_times_adj_n[i]            
            #print 'indxTmp: ',indxTmp
            #print '--------------------------------------------'
                        
            #swap_discount_factor[i] = z_out[indxTmp+1]
            swap_discount_factor[i] = z_out[indxTmp+1]
        
        
            #print 'swap_discount_factor: ', swap_discount_factor
            #print 'swap_times_adj: ', swap_times_adj
            #print 'regime_output: ', regime_output

        #swap_rates_out = compute_rates(swap_discount_factor, swap_times_adj, regime_output)
        
        #print 'swap_discount_factor : ', swap_discount_factor
        #print 'swap_times_adj_n: ', swap_times_adj_n
        #print 'regime_output: ', regime_output
        
        swap_rates_out = compute_rates(swap_discount_factor, swap_times_adj_n, regime_output)
        
        #print 'swap_rates_out: ', swap_rates_out


        df_array = swap_discount_factor
        r_array  = swap_rates_out
        t_array  = swap_times
        d_array  = swap_dates

        data_merge = merge_data(data_merge, d_array, t_array, r_array, df_array)

    else:
        
        swap_times = [999]
        swap_val   = [999]
        swap_discount_factor =[999]
            
        
    #%-------------------------------------------------------------------------
    #%-------------- GRAFICI DI CONTROLLO  ------------------------------------
    #%-------------------------------------------------------------------------
    
    merge_dates = data_merge['Dates'] 
    merge_times = data_merge['Times']
    merge_rates = data_merge['Rates']
    merge_df    = data_merge['Df']


    
    flag_make_graph = 1
    if (flag_make_graph == 1):

        graphrates(seg1_times[1:], seg1_values, futures_end_times, futures_rates, swap_times, swap_val, merge_times, merge_rates)
        graphdf(seg1_times, seg1_df, futures_times, futures_df, swap_times, swap_discount_factor)


 
    #%-------------------------------------------------------------------------
    #%-------------- SET OUTPUT  ROUTINES -------------------------------------
    #%-------------------------------------------------------------------------

    nodi_out = convertTimeToNode(merge_times)

    data_elab_out = {}
    
    data_elab_out['DateScadenza']    = merge_dates        
    data_elab_out['DiscountFactors'] = merge_df
    data_elab_out['TassiZC']         = merge_rates * 100.
    data_elab_out['Tempi']           = merge_times
    data_elab_out['Nodi']            = nodi_out
    data_elab_out['CodiceCurva']     = codeCurve
    
    return data_elab_out



if __name__ == "__main__":
    
    
    
    # ------------ CHK BOOTSTRAP -------------------------------------
        
    #inputCurveFile = r"input_test\raw_swap_curves_n2.txt"
    inputCurveFile = r"input_test\raw_short_swap_curves.txt"
    data_opt = set_data_opt()
    data_raw = load_data_fromFile(inputCurveFile)

    
    data_elab_out = fb.boot3s_elab_v2(data_opt, data_raw)

    print 'data_elab_out[DiscountFactors]: ', data_elab_out['DiscountFactors']
    print '-----------------------------------------------------------------------'
    print 'data_elab_out[TassiZC]: ', data_elab_out['TassiZC']        

    
    

    
    


