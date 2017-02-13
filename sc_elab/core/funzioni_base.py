

import sys
from mdates import holidays
from mdates import daycount
from mdates import busdayrule

import datetime as dtime


#import datetime


import numpy as np
import math



def FQ(label):
    print ('------------- FIN QUI TUTTO OK  %s ----------' %(label))
    sys.exit()




def compute_rates(df_in, time_df_in, regime_output):
    
    
    n_dmerge = len(df_in)
    out_rates = np.zeros(n_dmerge)
    
    if (regime_output == 1):     #% regime composto
        #regime = 'composto';
        for i in range(0, n_dmerge):
            out_r_tmp    = ( (1.0/df_in[i])**(1.0/time_df_in[i]) - 1.0) #% regime composto
            out_rates[i]= out_r_tmp

    elif (regime_output == 0): #% regime semplice
        #regime = 'semplice';
        for i in range(0, n_dmerge):
            out_r_tmp    = ( (1.0/df_in[i]) - 1.0 )/time_df_in[i] #% semplice
            out_rates[i]= out_r_tmp

    elif (regime_output == 2): #% regime continuo
        #regime = 'continuo';
        for i in range(0, n_dmerge):
            out_r_tmp    = - np.log(df_in[i])/time_df_in[i] #% semplice
            out_rates[i]= out_r_tmp
    
    return out_rates 




def fun_test_0():

    aa = 123
    bb = 126

    bb = 12312313126
    bb = 126
    
    
    AA = 'Modifica non esistente'
    
    'Speriamo si generi un conflitto!!'
    

    return 100

<<<<<<< HEAD

def purge_data(data_out, ref_field, val_not_allowed):
    
    data_out_n = {}
    key_list = data_out.keys()

    for kk in key_list:    
        data_out_n[kk] = []
    
    ln = len(data_out[ref_field])

    for kk in key_list:    
        new_vec = []
        
        for i in range(ln):
        
            valRefTmp = data_out[ref_field][i]
        
            if (valRefTmp == val_not_allowed):
                continue
            else:
                new_vec.append(data_out[kk][i])
         
        data_out_n[kk] = new_vec
    
    return data_out_n
        

def select_segments(data_raw_p):
    
    
    n_mat = len(data_raw_p['MatDate']) 
    
    valore_1s = []
    mat_1s    = []
    nodo_1s   = []

    valore_f = []
    mat_f    = []
    nodo_f   = []

    valore_s = []
    mat_s    = []
    nodo_s   = []
    
    data_1s ={}
    data_futures ={}
    data_swap ={}
    
    for i in range(0, n_mat):
        
        ts_tmp = data_raw_p['TipoSegmento'][i]
        vn_tmp = data_raw_p['ValoreNodo'][i]
        mat_tmp = data_raw_p['MatDate'][i]
        nd_tmp = data_raw_p['Nodo'][i]

        if (ts_tmp == 'D') or (ts_tmp == 'L'):

            valore_1s.append(vn_tmp)
            mat_1s.append(mat_tmp)
            nodo_1s.append(nd_tmp)
                    
        elif (ts_tmp == 'F'):
        
            valore_f.append(vn_tmp)
            mat_f.append(mat_tmp)
            nodo_f.append(nd_tmp)

        elif (ts_tmp == 'S'):      

            valore_s.append(vn_tmp)
            mat_s.append(mat_tmp)
            nodo_s.append(nd_tmp)

        else:


            print 'Caso non gestito!!!'   
            
        data_1s['ValoreNodo'] = valore_1s
        data_1s['MatDate'] = mat_1s
        data_1s['Nodo'] = nodo_1s
    
        data_futures['ValoreNodo'] = valore_f
        data_futures['MatDate'] = mat_f
        data_futures['Nodo'] = nodo_f
        
        data_swap['ValoreNodo'] = valore_s
        data_swap['MatDate'] = mat_s
        data_swap['Nodo'] = nodo_s
        
    return data_1s, data_futures, data_swap
    
def df_cont(spot_rate, time_ref):

    try:
        len(spot_rate)
        df_out = math.exp(-time_ref*spot_rate)

    except:
        
        df_out = np.ndarray(len(time_ref))
        for i in range(0, len(time_ref)):

            df_tmp = math.exp(-float(time_ref[i]*spot_rate[i]))
            df_out[i] = df_tmp

    return df_out



def boot3s_elab_n(data_opt, data_raw):
    
    
    ref_field       = 'UsaNodo'
    val_not_allowed = 'N'

    data_raw_p = purge_data(data_raw, ref_field, val_not_allowed)
    
    dict1s, dict_f, dict_s = select_segments(data_raw_p)
    
    
    refDate = data_opt['RefDate'] 

    
    rate1s = dict1s['ValoreNodo']
    rate2s = dict_f['ValoreNodo']
    rate3s = dict_s['ValoreNodo']

    dates1s = dict1s['MatDate']
    dates2s = dict_f['MatDate']
    dates3s = dict_s['MatDate']

    nodo1s = dict1s['Nodo']
    nodo2s = dict_f['Nodo']
    nodo3s = dict_s['Nodo']
    

    basis1s = data_opt['Basis']['D']    
    basis2s = data_opt['Basis']['F']
    basis3s = data_opt['Basis']['S']


    regime_output   = data_opt['RegimeOutput']
    
    if (len(dates1s)>0): dates1s.insert(0, refDate)
    if (len(rate1s)>0): rate1s.insert(0, rate1s[0])
    if (len(nodo1s)>0): nodo1s.insert(0, '0')

    if (len(dates2s)>0): dates2s.insert(0, refDate)
    if (len(dates3s)>0): dates3s.insert(0, refDate)
    
    times1s = daycount.yearfractions(dates1s, basis1s);
    times2s = daycount.yearfractions(dates2s, basis2s);
    times3s = daycount.yearfractions(dates3s, basis3s);
    
    if (len(times2s)>0): times2s = times2s[1:]
    if (len(times3s)>0): times3s = times3s[1:]

    
    df1s = df_cont(rate1s, times1s)
    df2s = df_cont(rate2s, times2s)
    df3s = df_cont(rate3s, times3s)
        
    if (len(dates2s)>0): dates2s = dates2s[1:]
    if (len(dates3s)>0): dates3s = dates3s[1:]

    
    rates1s = compute_rates(df1s[1:], times1s[1:], regime_output)
    rates2s = compute_rates(df2s, times2s, regime_output)
    rates3s = compute_rates(df3s, times3s, regime_output)

    #FQ(222)

    if (len(rates1s)>0): rates1s = np.insert(rates1s, 0, rates1s[0])
    if (len(rates2s)>0): rates2s = np.insert(rates2s, 0, rates2s[0])
    if (len(rates3s)>0): rates3s = np.insert(rates3s, 0, rates3s[0])

    merge_rates = []
    merge_dates = []
    merge_times = []
    merge_nodi  = []
    merge_df    = []


    if (len(rates1s)>0):
        for i in range(0, len(dates1s)):
            
            
            dateTmp  = dates1s[i]
            rateTmp  = rates1s[i]
            timesTmp = times1s[i]
            nodoTmp  = nodo1s[i]
            dfTmp    = df1s[i]
            
            merge_rates.append(rateTmp)
            merge_dates.append(dateTmp)
            merge_times.append(timesTmp)
            merge_df.append(dfTmp)
            merge_nodi.append(nodoTmp)
            



    if (len(rates2s)>0):

        for i in range(0, len(dates2s)):
            
            timesTmp = times2s[i]
            dateTmp = dates2s[i]
            rateTmp = rates2s[i]
            nodoTmp  = nodo2s[i]
            dfTmp   = df2s[i]
    
            merge_rates.append(rateTmp)
            merge_dates.append(dateTmp)
            merge_times.append(timesTmp)
            merge_df.append(dfTmp)
            merge_nodi.append(nodoTmp)
    
    if (len(rates3s)>0):

        for i in range(0, len(dates3s)):
            
            timesTmp = times3s[i]
            dateTmp = dates3s[i]
            rateTmp = rates3s[i]
            nodoTmp  = nodo3s[i]
            dfTmp   = df3s[i]
            
            merge_rates.append(rateTmp)
            merge_dates.append(dateTmp)
            merge_times.append(timesTmp)
            merge_df.append(dfTmp)
            merge_nodi.append(nodoTmp)



    data_elab_out = {}
    
    data_elab_out['DateScadenza']    = merge_dates        
    data_elab_out['DiscountFactors'] = merge_df
    data_elab_out['TassiZC']         = merge_rates
    data_elab_out['Tempi']           = merge_times
    data_elab_out['Nodi']            = merge_nodi

    
    return data_elab_out

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

        matDateTmp_dd = int(matDateTmp[0])
        matDateTmp_mm = int(matDateTmp[1])
        matDateTmp_yy = int(matDateTmp[2])
        
        
        dateTmp = dtime.date(matDateTmp_yy, matDateTmp_mm, matDateTmp_dd)

        """
        print 'dateTmp: ', dateTmp
        print 'nodoTmp: ', nodoTmp
        print 'matDateTmp: ', matDateTmp
        print 'valTmp: ', valTmp
        print 'usaNodoTmp: ', usaNodoTmp
        print 'tipoSegmentoTmp: ', tipoSegmentoTmp
        """

        
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



def set_data_opt():
    
    data_opt = {}
    
    data_opt['Basis'] = {}
    data_opt['Basis']['D']  = 'ACT/360'
    data_opt['Basis']['L']  = 'ACT/360'
    data_opt['Basis']['S']  = 'ACT/360'
    data_opt['Basis']['F']  = 'ACT/360'
    data_opt['Basis']['TN'] = 'ACT/360'
    data_opt['Basis']['O/N'] = 'ACT/360'
    
    data_opt['BasisFix'] = 'ACT/360'
    data_opt['BasisRef'] = 'ACT/ACT'
    
    data_opt['BusConv'] = {}
    data_opt['ParConvexity'] = {}
    data_opt['ParConvexity']['A'] = 1.0
    data_opt['ParConvexity']['B'] = 0.25

    data_opt['FutureTenor'] = 90

    data_opt['TenorSwap'] = 6
    data_opt['Convexity'] = True
    data_opt['GapFutures'] = True
    data_opt['SwapGapMethod'] = True
    data_opt['InterpLinFutures'] = True

    data_opt['RefDate'] = dtime.date(2016, 10, 02)

    data_opt['BusConv']['D'] = 'follow'
    data_opt['BusConv']['L'] = 'follow'
    data_opt['BusConv']['S'] = 'follow'
    data_opt['BusConv']['F'] = 'follow'
    data_opt['MKT']          = 'de'


    data_opt['MakeGraph'] = True
    data_opt['SaveGraph'] = True
    data_opt['RegimeOutput'] = 1

    data_opt['RegimeRate'] ={}
    data_opt['RegimeRate']['D'] = 0
    data_opt['RegimeRate']['L'] = 1


    return data_opt



if __name__ == "__main__":

    
    data_opt = {}
    data_raw = {}
    
    data_opt = set_data_opt()    
       
    #------------------------------------------------------------
    
    inputCurveFile = r"input_test\raw_swap_curves.txt"
    data_raw = load_data_fromFile(inputCurveFile)
    

    # -------- STRUTTURA DATA INPUT --------------------

    #print 'data_raw[Nodo]: ', data_raw['Nodo']
    

    
    """
    data_raw['UsaNodo'] = ['N', 'Y', ....]
    data_raw['Nodo']     = ['O/N', 'T/N', '1W', '2W', '1M', '2M', '3M', '4M', '5M'...]
    data_raw['ValoreNodo']  = [float, float, ..]
    data_raw['TipoSegmento']  = ['D', 'L', 'F', 'S']
    data_raw['MatDate']  = [datetime.date(2016, 10, 3), datetime.date(2016, 10, 4), .....]
    """
    
    
    
    data_elab_out = boot3s_elab_n(data_opt, data_raw)


    # -------- STRUTTURA DATA INPUT --------------------
    

    #print data_elab_out['DiscountFactors']
    
    print 'data_elab_out[Nodi]: ', data_elab_out['Nodi']
    #print 'data_elab_out: ', data_elab_out.keys()
    
    """
    data_elab_out['DiscountFactors'] = [float, float, float,...]
    data_elab_out['Tempi']     = [float, float, float,...]
    data_elab_out['DateScadenza']  =  [datetime.date(2016, 12, 30), datetime.date(2016, 10, 4), .....], prima data data di rieferimento
    data_elab_out['Nodi']  = ['0', '1M', '2M', '3M', '6M', '1Y', '10Y', '12Y', '15Y', '20Y', '25Y', '30Y', '40Y', '50Y']
    data_elab_out['TassiZC']  = [[float, float, float,...]
    """


    
    #FQ(332)

