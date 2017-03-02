import sys
import datetime as dtime
import math
import numpy as np

#from sc_elab import numpy as np
from sc_elab.core.mdates import holidays
from sc_elab.core.mdates import daycount
from sc_elab.core.mdates import busdayrule

from scipy import optimize
from scipy.optimize import minimize
from scipy.optimize import fmin


import matplotlib.pyplot as plt

import datetime
import calendar



def FQ(label):
    print ('------------- FIN QUI TUTTO OK  %s ----------' %(label))
    sys.exit()


def fitting():
    #---
    #esemplificativo dell'output richiesto per il fitting, da riempire solo i parametri del modello selezionato
    #---
    res = {'Dates'  : [datetime.date(2017,12,29), datetime.date(2017,12,30),datetime.date(2017,12,31)],
           'a'      : [0.0,0.0],
           'b'      : [0.9, 0.9],
           'c'      : [0.111, 0.111],
           'd'      : [0.222, 0.222],
           'e'      : [0.333, 0.333],
           'const1' : [0.44, 0.44],
           'const2' : [0.55, 0.55],
           'beta0'  : [0.1, 0.1],
           'beta1'  : [0.2, 0.2],
           'beta2'  : [0.3, 0.3],
           'beta3'  : [0.4, 0.4]
           }
    return res


def convertNodeToMnth(nodeRef):

    dict_N2T = {} 
    dict_N2T['1M'] = 1 
    dict_N2T['3M'] = 3
    dict_N2T['6M'] = 6
    dict_N2T['1Y'] = 1
    
    timeNode = dict_N2T[nodeRef]

    return timeNode


def convertTimeToNode(time_ref):
    

    node_out = []
    for i in range(len(time_ref)):
    
        time_tmp = time_ref[i]
        
        n_y = (time_tmp*1.0)
        n_m = float(time_tmp*12.0)
        n_w = float(time_tmp*52.0)
        #n_d = float(time_tmp*365.0)

        n_y_int = int(time_tmp)*1.0
        n_m_int = int(time_tmp)*12.0
        n_w_int = int(time_tmp)*52.0
        #n_d_int = int(time_tmp)*365.0

        #print 'n_m: ', n_m
        #print 'n_m_int: ', n_m_int
        
        if (time_tmp < 1.0/12.0):

            if abs(n_w - n_w_int)> 1:
    
                    unit = 'W'
                    n_unit = int(n_w)
            else:
                    unit = 'M'
                    n_unit = int(n_m)
        else:
        
            if abs(n_m - n_m_int)> 1:
    
                    unit = 'M'
                    n_unit = int(n_m)
            else:
                    unit = 'Y'
                    n_unit = int(n_y_int)

        nodeTmp = str(n_unit) + unit
        node_out.append(nodeTmp)
        
    return node_out 
        

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

        #print 'line_splittedTmp: ', line_splittedTmp
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

    
def select_segments(data_raw_p):
    
    
    n_mat = len(data_raw_p['MatDate']) 
    
    valore_1s = []
    mat_1s    = []
    nodo_1s   = []
    ts_1s     = []

    valore_f = []
    mat_f    = []
    nodo_f   = []
    ts_f   = []

    valore_s = []
    mat_s    = []
    nodo_s   = []
    ts_s   = []
    
    data_1s_dict ={}
    data_futures_dict ={}
    data_swap_dict ={}
    
    for i in range(0, n_mat):
        
        ts_tmp = data_raw_p['TipoSegmento'][i]
        vn_tmp = data_raw_p['ValoreNodo'][i]
        mat_tmp = data_raw_p['MatDate'][i]
        nd_tmp = data_raw_p['Nodo'][i]

        if (ts_tmp == 'D') or (ts_tmp == 'L'):

            ts_1s.append(ts_tmp)
            valore_1s.append(vn_tmp)
            mat_1s.append(mat_tmp)
            nodo_1s.append(nd_tmp)
                    
        elif (ts_tmp == 'F'):
        
            ts_f.append(ts_tmp)
            valore_f.append(vn_tmp)
            mat_f.append(mat_tmp)
            nodo_f.append(nd_tmp)

        elif (ts_tmp == 'S'):      

            ts_s.append(ts_tmp)
            valore_s.append(vn_tmp)
            mat_s.append(mat_tmp)
            nodo_s.append(nd_tmp)

        else:


            print 'Caso non gestito!!!'
            
            
            
        data_1s_dict['TipoSegmento'] = ts_1s
        data_1s_dict['ValoreNodo'] = np.asarray(valore_1s)
        data_1s_dict['MatDate'] = mat_1s
        data_1s_dict['Nodo'] = nodo_1s
    
        data_futures_dict['TipoSegmento'] = ts_f
        data_futures_dict['ValoreNodo'] = np.asarray(valore_f)
        data_futures_dict['MatDate'] = mat_f
        data_futures_dict['Nodo'] = nodo_f
        
        data_swap_dict['TipoSegmento'] = ts_s
        data_swap_dict['ValoreNodo'] = np.asarray(valore_s)
        data_swap_dict['MatDate'] = mat_s
        data_swap_dict['Nodo'] = nodo_s
        
     
    return data_1s_dict, data_futures_dict, data_swap_dict
    
        



def retrive_segment(data_raw, tipo_dato):
    
    date_out = []
    val_out = []
    
    for i in range(0, len(data_raw['MatDate'])):
        
        dateTmp = data_raw['TipoSegmento']
        tipoSTmp = data_raw['TipoSegmento']
        valTmp   = data_raw['TipoSegmento']
        
        if (tipoSTmp == tipo_dato):
        
            date_out.append(dateTmp)
            val_out.append(valTmp)

    return date_out, val_out



def increase_datetime_list(date_list_start, n_days):
    
    
    date_list_end = []
    
    for i in range(0, len(date_list_start)):
        
        dateTmp = date_list_start[i] + datetime.timedelta(days = n_days)
        date_list_end.append(dateTmp)
        
        
        
    return date_list_end
        
        
        
    

def map_segment(data_raw):    
                
    list_s = set(data_raw['TipoSegmento'])
    list_s = list(list_s)
    
    n_s    = len(list_s)

    seg_info = {}
    for i in range(n_s):
        
        TmpS = list_s[i]
        if (TmpS == 'D'):
            seg_info['1'] = 'D'
            
        elif (TmpS == 'L'):
            seg_info['1'] = 'L'

        elif (TmpS == 'F'):
            seg_info['2'] = 'F'

        else: 
        
            if (n_s == 3):
                seg_info['3'] = 'S'
            else:
                seg_info['2'] = 'S'

 
    return seg_info, n_s



def df_simple(spot_rate, time_ref):
    

    try:
        df_out = 1.0/(1.0 + spot_rate*time_ref)

    except:
        
        df_out = np.ndarray(len(time_ref))
        for i in range(0, len(time_ref)):

            df_tmp = 1.0/(1.0 + float(spot_rate[i])*float(time_ref[i]))
            df_out[i] = df_tmp

    return df_out
    
def df_cont(spot_rate, time_ref):

    try:
        df_out = math.exp(-time_ref*spot_rate)

    except:
        
        df_out = np.ndarray(len(time_ref))
        for i in range(0, len(time_ref)):

            df_tmp = math.exp(-float(time_ref[i]*spot_rate[i]))
            df_out[i] = df_tmp

    return df_out

def df_cmp(spot_rate, time_ref):

    try:
        df_out = 1.0/(1 + spot_rate)**(time_ref)

    except:
        
        df_out = np.ndarray(len(time_ref))
        for i in range(0, len(time_ref)):

            df_tmp = 1.0/(1.0 + spot_rate[i])**time_ref[i]
            df_out[i] = df_tmp

    return df_out


def find_indx(x_target, x_list_ref):

    ln = len(x_list_ref) - 1

    for i in range(0, len(x_list_ref)):

        if (x_list_ref[0] > x_target): 
            i_ref = None
            break
         
        elif (x_list_ref[i] == x_target) and (i == 0): 
            i_ref = i
            break

        elif (x_list_ref[ln] < x_target): 
            i_ref = None
            break
        
        elif (x_list_ref[i] >= x_target):
            i_ref = i - 1
            break
            
        else:

            continue

    return i_ref


def find_indx_toll(x_target, x_list_ref, toll):

    ln = len(x_list_ref) - 1

    for i in range(0, len(x_list_ref)):

        if (abs(x_list_ref[i] - x_target) < toll ): 
            i_ref = i
            break

        elif (x_list_ref[0] > x_target): 
            i_ref = 0
            break
         
        elif (x_list_ref[ln] < x_target): 
            i_ref = ln
            break

        else:

            continue

    return i_ref



def find_indx_n(x_target, x_list_ref):

    ln = len(x_list_ref) - 1

    for i in range(0, len(x_list_ref)):

        if (x_list_ref[0] > x_target): 
            i_ref = None
            break

        elif (x_list_ref[ln] < x_target): 
            i_ref = None
            break
 
        elif (x_list_ref[i] == x_target): 
            i_ref = i
            break

        elif (x_list_ref[i] > x_target):
            i_ref = i - 1
            break
            
        else:

            continue

    return i_ref




def find_indx_equal(x_target, x_list_ref):

    
    i_ref = None

    for i in range(0, len(x_list_ref)):


        if (x_list_ref[i] == x_target):
            i_ref = i
            break
            
        else:

            continue

    return i_ref






def df_from_interp_df_exp(t_target, df_n, df_o, t_n, t_o):

    A =  df_o           
    B = (1.0/(t_o - t_n))*np.log(df_n/df_o)
    
    df_target = A*np.exp(-B*(t_target - t_o))

    return   df_target 





def df_from_interp_R_lin(t_target, df_n, df_o, t_n, t_o):
              
    r_n = np.log(df_n)/t_n
    r_o = np.log(df_o)/t_o
    
    r_target = r_o + (r_n - r_o)*(t_target - t_o)/(t_n - t_o)
    df_target = np.exp(-r_target*t_target)

    return df_target



def add_months(sourcedate, months):
    month = sourcedate.month - 1 + months
    year = int(sourcedate.year + month / 12 )
    month = month%12 + 1
    day = min(sourcedate.day,calendar.monthrange(year,month)[1])
    return datetime.date(year,month,day)


def compute_df_future(seg1_times, seg1_val, seg1_df,  futures_rates, flag_futures_gap, futures_start_time, futures_end_time,  futures_end_df_x, flag_interp1, indx_ref):

    
    t_s_target   = futures_start_time[indx_ref]
    t_e_target   = futures_end_time[indx_ref]
    
    
    index_seg1_last_s = find_indx_toll(t_s_target, seg1_times, 0.1)
    index_seg1_last_e = find_indx_toll(t_e_target, seg1_times, 0.1)

    #index_seg1_last_s = find_indx_n(t_s_target, seg1_times)
    #index_seg1_last_e = find_indx_n(t_e_target, seg1_times)
    
    
    index_cross     = find_indx_equal(futures_start_time[indx_ref], seg1_times)
    


    if (indx_ref == 0):

        print 'AAAA'
        """
        t_o  = seg1_times[index_seg1_last_s]
        t_oo  = seg1_times[index_seg1_last_s-1]

        df_o = seg1_df[index_seg1_last_s]
        df_oo = seg1_df[index_seg1_last_s-1]
        """
    else:

        df_o = futures_end_df_x[indx_ref-1]
        t_o  = futures_end_time[indx_ref-1]
        
        futures_rates_old = futures_rates[indx_ref - 1]

    futures_rates_target = futures_rates[indx_ref]
    
    
    if (index_cross != None):
        futures_start_1 = seg1_df[index_cross]
        
        
        futures_start_df_0 = futures_start_1
         
        futures_end_df_0 = futures_start_1/(1.0 + futures_rates_target*(futures_end_time[indx_ref] - t_s_target) )
        
        

    elif (seg1_times[len(seg1_times)-1] > t_s_target): #% controllo della sovrapposizione: seg1 e futures

        t_o_s  = futures_end_time[indx_ref-1]
        t_o_e  = futures_end_time[indx_ref-1]

        #t_o_s  = seg1_times[index_seg1_last_s]
        t_n_s  = seg1_times[index_seg1_last_s + 1]
        df_n_s = seg1_df[index_seg1_last_s + 1] 
        df_o_s = seg1_df[index_seg1_last_s]        
        
        if (index_seg1_last_e != len(seg1_times) - 1):    

            t_o_e  = seg1_times[index_seg1_last_e] 
            t_n_e  = seg1_times[index_seg1_last_e + 1]
            df_n_e = seg1_df[index_seg1_last_e + 1] 
            df_o_e = seg1_df[index_seg1_last_e]        


        
        
        t_o_s1_s  = seg1_times[index_seg1_last_s]
        df_o_s1_s = seg1_df[index_seg1_last_s]

        t_o_s1_e  = seg1_times[index_seg1_last_e]
        df_o_s1_e = seg1_df[index_seg1_last_e]
        
        if (t_o_s1_s > t_o_s):
            
            t_o_s  = t_o_s1_s
            df_o_s = df_o_s1_s
        

        if (index_seg1_last_e != None) and (t_o_s1_e > t_o_e):
            
            t_o_e  = t_o_s1_e
            df_o_e = df_o_s1_e
        


        if (flag_interp1 == 1): #%----- SI SOVRAPPOSIZIONE

            #%----------intepolazione esponenziale sui fattori di scont ---------------
            futures_start_df_0 =  df_from_interp_df_exp(t_s_target, df_n_s, df_o_s, t_n_s, t_o_s)

            if (seg1_times[len(seg1_times)-1] > t_e_target):
                
                futures_end_df_0 = df_from_interp_df_exp(t_e_target, df_n_e, df_o_e, t_n_e, t_o_e)
                
            else:
                futures_rates_target = futures_rates[indx_ref]                    
                futures_end_df_0 = futures_start_df_0/(1.0 + futures_rates_target*(futures_end_time[indx_ref] - t_s_target))


        else: #%-------------- interpolazione lineare sui tassi -----------------------
            
            futures_start_df_0 = df_from_interp_R_lin(t_s_target, df_n_s, df_o_s, t_n_s, t_o_s)

        
            if (seg1_times[len(seg1_times)-1] > t_e_target):
                
                futures_end_df_0 = df_from_interp_R_lin(t_e_target, df_n_e, df_o_e, t_n_e, t_o_e)

            else:        
                futures_rates_target = futures_rates[indx_ref]                    
                futures_end_df_0 = futures_start_df_0/(1.0 + futures_rates_target*(futures_end_time[indx_ref] - t_s_target))
        
        
        

    else:  # ---- NO SOVRAPPOSIZIONE

        #t_n = futures_start_time[indx_ref + 1]
        t_o = futures_start_time[indx_ref - 1]

        df_o = futures_end_df_x[indx_ref-1]
        t_o  = futures_end_time[indx_ref-1]

        
        #--------------------------

        #df_n = seg1_df[index_seg1_last + 1] 
        futures_rates_target = futures_rates[indx_ref]
        
        #if (indx_ref == 0):
        
        if (indx_ref == 0) and (len(seg1_val) == 0):
                
            R_start1 = futures_rates_target
            futures_start_df_0 = np.exp(-R_start1*futures_start_time[indx_ref])
        
        elif (flag_futures_gap == 1) and (len(seg1_val) > 0): #------------ gap colmato mantenedo costante il tasso spot

            R_start1            = -np.log(df_o)/t_o #CHK indici 
            
            futures_start_df_0 =  np.exp(-R_start1*t_s_target)
            
        else:  #------------ gap colmato mantenedo costante il tasso fwd ---------------------------------------------
            
            #fwd_start1         = -np.log(df_n/df_o )/(t_n - t_o)
            
            if (indx_ref == 0):
                
                t_o  = seg1_times[index_seg1_last_s]
                t_oo  = seg1_times[index_seg1_last_s-1]
        
                df_o = seg1_df[index_seg1_last_s]
                df_oo = seg1_df[index_seg1_last_s-1]
                
                fwd_start1         = -np.log(df_o/df_oo)/(t_o - t_oo)
                futures_start_df_0 = df_o*np.exp(-fwd_start1*(t_s_target - t_o ))

            else:
                fwd_start1         = futures_rates_old
                futures_start_df_0 = df_o/(1.0 + fwd_start1*(t_s_target - t_o ))
            
        


        futures_end_df_0 = futures_start_df_0/(1 + futures_rates_target*(t_e_target - t_s_target) )

    futures_end_df_x[indx_ref] = futures_end_df_0


    #dr = 1.0/(t_e_target - t_s_target)*(futures_start_df_0/futures_end_df_0 - 1.0)
    
    
    return futures_start_df_0, futures_end_df_0, futures_end_df_x 



def compute_convexity(par_convexity, futures_start_tx, futures_end_tx):
    
    nf_n = len(futures_start_tx) 
    conv_correction  = np.ndarray(nf_n)
    mean_convexity   = par_convexity[0]
    vol_convexity    = par_convexity[1]
        
    for i in range(0, nf_n):

        t1 = futures_start_tx[i]
        t2 = futures_end_tx[i]

        if (mean_convexity == 0): #---------limite a-->0 (Hull& White) --> (Ho Lee)
            conv_corrTmp = 0.5*t1*t2*vol_convexity**2
        else:
            #%--------------- Calcolo della correzione mediante il modello di Hull&White
            
            B1  = ( 1 - np.exp(- mean_convexity*(t2 - t1)) )/mean_convexity
            B2  = ( 1 - np.exp(- mean_convexity*(t1) ) )/mean_convexity
            conv_corrTmp  = B1/(t2 - t1)*( B1*(1 - np.exp(-2*mean_convexity*t1)) + 2*mean_convexity*B2**2)*vol_convexity**2/(4*mean_convexity)

        conv_correction[i] = conv_corrTmp
            
            
    return conv_correction

def compute_rates(df_in, time_df_in, regime_output):
    
    
    regime_output = int(regime_output)
    
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


def from_date_to_ordinal(date_dates):


    serial_dates = []

    for i in range(0, len(date_dates)):
        
        
        dateTmp = date_dates[i].toordinal()
        
        serial_dates.append(dateTmp)


    return  serial_dates  

def compute_first_df(seg1_df, seg1_dates, seg1_times, seg1_val, fix_time1, on_date, tn_date): 


    seg1_df_0 = df_simple(seg1_val[0], seg1_times[0])
    
    if ( (seg1_dates[0] != on_date) and (seg1_dates[0] != tn_date) ):

        seg1_df[0] = 1/(1 + seg1_val[0]*seg1_times[0]);
        spot_seg1_rate_1 = -math.log(seg1_df_0)/seg1_times[0];
        df_fix_seg1 = math.exp(-spot_seg1_rate_1*fix_time1);
        id_start = 1

    #---------------- caso 2) si ON, si TN ----------------------------------------------------
    elif (seg1_dates[0] == on_date) and (seg1_dates[1] == tn_date):

        
        seg1_df[0] = 1/(1 + seg1_val[0]*seg1_times[0])
        seg1_df[1] = seg1_df_0*df_simple(seg1_val[1], seg1_times[1] - seg1_times[0])
        df_fix_seg1 = seg1_df[1]
        id_start = 2

    #---------------- caso 3) no ON, si TN ---------------------------------------------------
    elif (seg1_dates[0] != on_date) and (seg1_dates[1] == tn_date):

        seg1_df[0] = seg1_df_0
        df_fix_seg1 = seg1_df[0]
        id_start = 1;

    #---------------- caso 4) si ON, no TN ---------------------------------------------------
    elif (seg1_dates[0] == on_date) and (seg1_dates[1]!= tn_date):

        #seg1_df(1) = 1./( 1 + seg1_val(1).*seg1_times(1) );
        seg1_df[0] = seg1_df_0
        seg1_df[1] = seg1_df_0*df_simple(seg1_val[1], seg1_times[1] - seg1_times[0])
        
        A = seg1_df[0]
        B = 1.0/(seg1_times[0] - seg1_times[1])*math.log(seg1_df[1]/seg1_df[0])
        df_fix_seg1 = A*math.exp(-B*(fix_time1 - seg1_times[0]))
        id_start = 2



        
    else:
        
        df_fix_seg1 = seg1_df_0
        id_start = 1

    return seg1_df, df_fix_seg1, id_start 



def fwd_fun(x,n,m,f1, py_rate_m, time_new):

    #%--------------------------------------------------------------------------
    f0 = 0
    for k in range(n,m):    
        
        f0 = (time_new[k] - time_new[k-1])*np.exp(-x*(time_new[k]-time_new[n])) + f0
    
    #f = py_rate_m*f0 + np.exp(-x*(time_new[m]-time_new[n])) - f1


    f = f1 - x
    
    #%-------------------------------------------------------------------------
    return f



def fwd_fun_0(x, n, m, f1, py_rate_m, time_new):

    #%--------------------------------------------------------------------------
    f0 = 0
    
    for k in range(n+1,m+1):    
        
        #print 'k: ', k
        dt_tmp = (time_new[k] - time_new[k-1])
        t_tmp = time_new[k]-time_new[n]
        n_tmp = k-n
        z_Tmp = (1.0 + x*dt_tmp)**-(n_tmp)
        f0 = dt_tmp*z_Tmp + f0

        """        
        print 't_tmp: ', t_tmp
        print 'dt_tmp: ', dt_tmp
        print 'z_Tmp: ', z_Tmp
        print 'f0: ', f0
        
        print '-----------------------------'

        #f = py_rate_m*f0 + np.exp(-x*(time_new[m]-time_new[n])) - f1
        """

    X = py_rate_m*f0
    Y = (1.0 + x*dt_tmp)**-(m-n)
    

    #Y = np.exp(-x*(time_new[m]-time_new[n]))
    
    f_out = X  + Y - f1

    #print 'X: ', X
    #print 'Y: ', Y
    #print 'f_out: ', f_out
    
    
    #%-------------------------------------------------------------------------
    return f_out





def conv_to_dates(datetimes_to_convert):
    
    converted_dates = []
    
    ln = len(datetimes_to_convert)

    for i in range(0, ln):
        
        convertedTmp = datetimes_to_convert[i].date()
        
        converted_dates.append(convertedTmp)

    return converted_dates


def merge_data(data_to_merge, d_array, t_array, r_array, df_array):


    times_list =  data_to_merge['Times']
    dates_list =  data_to_merge['Dates']
    rates_list =  data_to_merge['Rates']
    df_list    =  data_to_merge['Df']

    for i in range(0, len(d_array)):
        
        d_tmp = d_array[i]
        t_tmp = t_array[i]
        r_tmp = r_array[i]
        df_tmp = df_array[i]
        
        times_list = np.append(times_list, t_tmp)
        dates_list = np.append(dates_list, d_tmp)
        rates_list = np.append(rates_list, r_tmp)
        df_list    = np.append(df_list, df_tmp)

    values_to_fill = []

    for i in range(0, len(times_list)):

        t_tmp = times_list[i]
        d_tmp = dates_list[i]
        r_tmp = rates_list[i]
        df_tmp = df_list[i]
        

        listTmp = (t_tmp, d_tmp, df_tmp, r_tmp)   
        values_to_fill.append(listTmp)

    dtype = [('Times', float), ('Dates', datetime.date), ('Df', float), ('Rates', float)]

    data_out = np.array(values_to_fill, dtype=dtype)
    data_out = np.sort(data_out, order='Times') 
    
    
    return  data_out

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
        



def graphrates(dep_times, dep_rates, fu_times, fu_rates, sw_times, sw_rates, boot_tims, boot_rates):


    g1 = plt.figure()
    plt.plot(dep_times, dep_rates, '-o', label='depositi')

    if (fu_times[0] != 999): plt.plot(fu_times, fu_rates, '-o', label='futures')   
    if (sw_times[0] != 999): plt.plot(sw_times, sw_rates,'-o', label='swap')
    
    plt.plot(boot_tims, boot_rates,'--o', label='bootstrap')
    
    
    plt.xlabel('Tempi [anni]')
    plt.ylabel('Livello tassi')
    plt.grid()
    plt.legend()
    plt.show()
    
    return g1



def graphdf(dep_times, dep_df, fu_times, fu_df, sw_times, sw_df):

    g1 = plt.figure()

    plt.plot(dep_times, dep_df, 'o-', label='depositi')

    if (fu_times[0] != 999): plt.plot(fu_times, fu_df, 'o-', label='futures')
    if (sw_times[0] != 999): plt.plot(sw_times, sw_df,'o-', label='swap')
    
    
    plt.xlabel('Tempi [anni]')
    plt.ylabel('Livello fattore di sconto')
    plt.grid()
    plt.legend()
    plt.show()
    return g1







def compute_z_from_swap_rate(z_all, df_fix_swap, tenors, swap_val_new, tenor_ref):
            
    sum_z      = np.sum(z_all*tenors)
    z_new      = (df_fix_swap - swap_val_new*sum_z)/(1.0 + swap_val_new*tenor_ref)

    return z_new 



def compute_z_from_cfr(z_out, tenors, swap_val, times_swap_new, t_swt_n, t_swt_m, times_swap_, df_fix_swap):
    
    
    nz  = find_indx_n(t_swt_n, times_swap_) #% indice del vettore times_out corrispondente a n anni
    mz  = find_indx_n(t_swt_m, times_swap_) #% indice del vettore times_out corrispondente a m anni

    #fid_cond = (z_out != 0)
    #index_stop = find_indx_p(False, fid_cond) #% selezione dell'(n-1)-mo indice
    
    index_stop = nz + 1
    
    z_n     = z_out*tenors
    z_n     = z_n[:index_stop]
    sum_z   = np.sum(z_n)          #% sommatoria di tutti gli z fino all'(n-1)-mo
    
    z_nz = z_out[nz]    
    z_fix = z_out[0]

    f1    = 1.0/z_nz*(z_fix - swap_val*sum_z) #% Def. f1 per il calcolo di F

    
    #%f1 rappresenta il termine noto della equazione (2.10) pag.38 del
    #%documento modelli_tassi_v02.pdf
    
    x_0   = 0.001
    
    
    
    fwd_target = optimize.newton(fwd_fun_0, x_0, args=(nz, mz, f1, swap_val, times_swap_new))

    
    for k in range(nz+1,mz+1):
        z_out[k] = z_out[nz]*np.exp(-fwd_target*(times_swap_new[k] - times_swap_new[nz])) #; % calcolo il fattore di sconto z_out
    
    
    return z_out

"""
def boot3s_elab(data_opt, data_raw):

    #%--------------------------------------------------------------------------
    #%-------------- RETRIEVE DATA  -----------------------------------------
    #%--------------------------------------------------------------------------


    ref_field = 'UsaNodo'
    val_not_allowed = 'N'

    data_raw_p = purge_data(data_raw, ref_field, val_not_allowed)
    
    dict1s, dict_f, dict_s = select_segments(data_raw_p)

    seg1_dates  = dict1s['MatDate']
    seg1_values = dict1s['ValoreNodo']

    futures_start_dates = dict_f['MatDate']
    futures_values      = dict_f['ValoreNodo']
    
    swap_dates          = dict_s['MatDate']
    swap_val            = dict_s['ValoreNodo']

    n1s  = len(seg1_dates)
    nf   = len(futures_start_dates)
    nsw = len(swap_dates)
    

    #%--------------------------------------------------------------------------
    #%-------------- RETRIEVE OPTIONS SETUP -----------------------------------------
    #%--------------------------------------------------------------------------

    setting_default = set_data_default()

    ref_date = data_opt['RefDate'] 
    
    par_a = data_opt['ParConvexity']['A']
    par_b = data_opt['ParConvexity']['B']

    par_convexity    = [par_a, par_b]

    tenor_swap       = data_opt['TenorSwap']
    flag_convexity   = data_opt['Convexity']
    flag_futures_gap = data_opt['GapFutures']
    flag_method_swap = data_opt['SwapGapMethod']
    flag_interp1     = data_opt['InterpLinFutures']

    flag_make_graph = data_opt['MakeGraph']
    flag_save_graph = setting_default['SaveGraph']
    regime_output   = data_opt['RegimeOutput']   
    future_tenor    = data_opt['FutureTenor'] 
    
    '------- retrive conventions ------------------------'

    mkt_code = data_opt['MKT']
    mkt_ref  = holidays.get_calendar(mkt_code)
    
    basis_s = data_opt['Basis']['S']
    basis_f = data_opt['Basis']['F']

    basis_f = convert_basis(basis_f)
    basis_s = convert_basis(basis_s)

    day_conv_f = data_opt['BusConv']['F']
    day_conv_s = data_opt['BusConv']['S']
    day_conv_tn = data_opt['BusConv']['TN']
    day_conv_on = data_opt['BusConv']['O/N']

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
    
    if (len(futures_start_dates) > 0):


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
    
        
    
    if (len(futures_values) > 0):

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



    if (len(swap_dates)> 0):

        
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

    if (nf > 0):
  
        futures_rates = (100.0 - futures_values)/100.0  #-------- tasso future
    
        #-------- Correzione: "Convexity"------------------------------------------------'
        
        if (flag_convexity == 1):
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
            
        
        future_rates = compute_rates(futures_df, futures_times, regime_output)
    
        df_array = futures_df
        r_array  = future_rates
        t_array  = futures_times
        d_array  = futures_dates

        '------------------- merge ------------------'

        data_merge = merge_data(data_merge, d_array, t_array, r_array, df_array)
    
    else:

        futures_start_times = [999]
        futures_times       = [999]
        futures_rates       = [999]
        futures_df          = [999]
        
        #%----------------------------------------------------------------------------------------------
        #%-------------- ELABORAZIONE SEGMENTO SWAP: FATTORI DI SCONTO A PARTIRE DAI TASSI SWAP ---------
        #%----------------------------------------------------------------------------------------------


    merge_dates = data_merge['Dates']
    merge_times = data_merge['Times']
    merge_df    = data_merge['Df']

    #print'merge_times: ', len(merge_times) 
    
    swap_discount_factor = np.zeros(len(swap_dates))
    
    
    if (nsw > 0):


        index_last = find_indx(fix_swap_date, merge_dates)

        t_o =  merge_times[index_last]
        t_n =  merge_times[index_last+1]

        df_o = merge_df[index_last]
        df_n = merge_df[index_last + 1]
        
        t_target = fix_swap_time
        
        df_fix_swap = df_from_interp_df_exp(t_target, df_n, df_o, t_n, t_o)

    
        #%--------------------------- Def. variabili ----------------------------------------------------------
    
        dt_time_swap     = swap_times[1:]- swap_times[0:len(swap_times)-1]  # differenze tempi
        
        
        index_condition  = dt_time_swap.round(0) > 1
        index_irregular  = np.where(index_condition)  # vettore indici date "irregolari" (array, dtype=int64)        
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
        index_start = find_indx(swap_times[0], times_swap_new) # indice di inizio del segmento degli swap

        times_swap_tmp1s  = times_swap_new[:index_start + 1]   # tempi fino al primo tasso swap
        times_swap_tmp2s  = times_swap_new[index_start  + 1:]  # tempi oltre il primo tasso swap

        dates_swap_tmp1s  = swap_dates_new[:index_start + 1]   # date fino al primo tasso swap        
        dates_swap_tmp2s  = swap_dates_new[index_start  + 1:]  # date oltre il primo tasso swap

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
    
        
        
        z_all        = np.concatenate((df_tmp1s, np.zeros(nsw2s)))
        z_swap2s     = np.zeros(nsw2s) #% fattori di sconto a partire dal primo tasso swap
    
        #%-------------------------------------------------------------------------------------
        #%--------------- Generazione dei fattori di sconto alle date dei tassi swap ----------
        #%-------------------------------------------------------------------------------------
    
    
        for i in range(0,nsw2s):
            sum_z            = np.sum(z_all*tenors);
            z_swap2s[i]      = (df_fix_swap - swap_val_new[i]*sum_z)/(1.0 + swap_val_new[i]*tenors[nsw1s + i])
            z_all[nsw1s + i] = z_swap2s[i]
    
        #%--------------------------------  set LSR method OUTPUT --------------------------------
    
        if (flag_method_swap == 1): #% considero il metodo: Linear Swap Rate: 'LSR'
            #swap_df_tmp    = np.insert(z_all, 0, 1)
            swap_times_tmp = np.insert(times_swap_new, 0, 0)
            swap_times_tmp = np.round(swap_times_tmp, 2)
            z_out = z_all[1:]

            #%------------------------------------------------------------------------------------
            #%--------------------------------  CFR METHOD ----------------------------------------%
            #%------------------------------------------------------------------------------------
        
        elif (flag_method_swap == 0) and (ln_irr > 1): # Considero il metodo "Constant fwd rate" 'CFR'
            
            n_irr = len(index_irregular) # %n. date irregolari
            index_irregular_first = int(np.round(swap_times_x[index_irregular[0]]/(tenor_swap/12.0))) # % indice prima data irregolare
            index_irregular_last  = int(np.round(swap_times_x[index_irregular[len(index_irregular)-1]]/(tenor_swap/12.0))) #; % indice ultima data irregolare


            index_irregular_0     = max(len(df_tmp1s), index_irregular_first)
            index_irregular_last  =  index_irregular_last

            z_out = np.zeros(max(1, index_irregular_last))
                
            z_out[:index_irregular_0] = z_all[1: index_irregular_0 + 1]

            times_swap_       = np.round(times_swap_x_new/(tenor_swap/12.0))*(tenor_swap/12.0)
             
                        
            times_swap_new  = times_swap_new + fix_swap_time


            tenors = tenors[1:]
    
            for i in range(0, n_irr - 1):     #% ciclo su tutte le date irregolari
                
    
                n_i = index_irregular[i]     #% indice dell'n-ma data "irregolare"
                m_i = index_irregular[i+1]   #% indice dell'n-ma data "irregolare"
                
                
                swap_times_n  = swap_times_x[n_i]  #% tempo dell'n-ma data "irregolare"
                swap_times_m  = swap_times_x[m_i]  #% tempo dell'n-ma data "irregolare"
    
                t_swt_n = np.round(swap_times_n)   #% n. anni dell'n-ma data "irregolare"
                t_swt_m = np.round(swap_times_m)   #% m. anni dell'm-ma data "irregolare"

                
                # ------VERSIONE 1

                z_out = compute_z_from_cfr(z_out, tenors, swap_val[m_i-1], times_swap_new, n_i, m_i, t_swt_n, t_swt_m, times_swap_, df_fix_swap)
    

        
        else:
            
            z_out = z_all[1:]
            times_swap_       = np.round(times_swap_x_new/(tenor_swap/12.0))*(tenor_swap/12.0)            


            
        #%---------------------------------------------------------------------
        #%------------------ set swap output ----------------------------------
        #%---------------------------------------------------------------------

        #times_swap_       =  times_swap_[1:]           
        times_swap_       = np.round(times_swap_x_new/(tenor_swap/12.0))*(tenor_swap/12.0)            

        swp__ = np.round(swap_times_x/(tenor_swap/12.0), 0)
        swap_times_adj = np.round(swp__*(tenor_swap/12.0),1)
        
        if (swap_times_adj[0] == 0): swap_times_adj = swap_times_adj[1:]
        
    
        
        for i in range(0, nsw):
            
            indxTmp = find_indx_n(swap_times_adj[i], times_swap_[1:])
            
            
            swap_discount_factor[i] = z_out[indxTmp]
            
        swap_rates_out = compute_rates(swap_discount_factor, swap_times_adj, regime_output)


        df_array = swap_discount_factor
        r_array  = swap_rates_out
        t_array  = swap_times
        d_array  = swap_dates

        data_merge = merge_data(data_merge, d_array, t_array, r_array, df_array)

    else:
        
        swap_times = [999]
        swap_val   = [999]
            
        
    #%-------------------------------------------------------------------------
    #%-------------- GRAFICI DI CONTROLLO  ------------------------------------
    #%-------------------------------------------------------------------------
    
    merge_dates = data_merge['Dates'] 
    merge_times = data_merge['Times']
    merge_rates = data_merge['Rates']
    merge_df    = data_merge['Df']

    
    if (flag_make_graph == 1):

    
        g1 = graphrates(seg1_times[1:], seg1_values, futures_start_times, futures_rates, swap_times, swap_val, merge_times, merge_rates)
        g2 = graphdf(seg1_times, seg1_df, futures_times, futures_df, swap_times, swap_discount_factor)

        '--------------- save graph ------------------------------------------------------------'

    
    if (flag_save_graph == 1):
    
        
        g1.savefig('fig/plot_tassi_%s.png' %(ref_date))
        g2.savefig('fig/plot_discount_%s.png' %(ref_date))
 
    #%-------------------------------------------------------------------------
    #%-------------- SET OUTPUT  ROUTINES -------------------------------------
    #%-------------------------------------------------------------------------

    nodi_out = convertTimeToNode(merge_times)

    data_elab_out = {}
    
    data_elab_out['DateScadenza']    = merge_dates        
    data_elab_out['DiscountFactors'] = merge_df
    data_elab_out['TassiZC']         = merge_rates
    data_elab_out['Tempi']           = merge_times
    data_elab_out['Nodi']            = nodi_out
    
    return data_elab_out

"""

def chkDataCoherence(dict1s, dict_f, dict_s):
    
    s1 = list(set(dict1s['TipoSegmento']))
    #s2 = list(set(dict_f['TipoSegmento']))
    #s3 = list(set(dict_s['TipoSegmento']))

    if (len(s1)< 1):  # caso-1mo segmento non valorizzato       
        msg = 'Missing first Segm...  Check your input!'
        print msg
        raise ValueError(msg)
    
    if (len(s1)> 1):
        msg = 'Found BOTH Dept and Libor Rates... check your input!'
        print msg
        raise ValueError(msg)

    bool_stessi_nodi = False
    for nTmp in dict1s['Nodo']:
        
        if (nTmp in dict_s['Nodo']):
            
            bool_stessi_nodi = True
            nodo_sovrapposto = nTmp 
            
    if (bool_stessi_nodi == True):  # caso-1mo segmento non valorizzato       
        msg = 'Segms overlapping (node: %s)... check your input!'%(nodo_sovrapposto)
        print msg
        raise ValueError(msg)
        
def boot3s_elab_v2(data_opt, data_raw):

    #%--------------------------------------------------------------------------
    #%-------------- RETRIEVE DATA  -----------------------------------------
    #%--------------------------------------------------------------------------


    ref_field = 'UsaNodo'
    val_not_allowed = 'N'

    data_raw_p = purge_data(data_raw, ref_field, val_not_allowed)
    
    dict1s, dict_f, dict_s = select_segments(data_raw_p)
    
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

    futures_start_dates = dict_f['MatDate']
    futures_values      = dict_f['ValoreNodo']
    
    swap_dates          = dict_s['MatDate']
    swap_val            = dict_s['ValoreNodo']

    seg1_values = seg1_values/100.0
    swap_val = swap_val/100.0

    
    n1s  = len(seg1_dates)
    nf   = len(futures_start_dates)
    nsw = len(swap_dates)
    
    

    #%--------------------------------------------------------------------------
    #%-------------- RETRIEVE OPTIONS SETUP -----------------------------------------
    #%--------------------------------------------------------------------------

    setting_default = set_data_default()

    ref_date = data_opt['RefDate'] 
    
    par_a = data_opt['ParConvexity']['A']
    par_b = data_opt['ParConvexity']['B']

    par_convexity    = [par_a, par_b]

    tenor_swap       = data_opt['TenorSwap']
    tenor_swap       = convertNodeToMnth(tenor_swap)
    
    
    
    #print 'tenor_swap: ', tenor_swap


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
    
    basis_s = data_opt['Basis']['S']
    basis_f = data_opt['Basis']['F']

    basis_f = convert_basis(basis_f)
    basis_s = convert_basis(basis_s)

    day_conv_f = data_opt['BusConv']['F']
    day_conv_s = data_opt['BusConv']['S']
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
    
    if (len(futures_start_dates) > 0):


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
        
        
    
    if (len(futures_values) > 0):

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



    if (len(swap_dates)> 0):

        
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

    if (nf > 0):
  
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

        futures_start_times = [999]
        futures_times       = [999]
        futures_rates       = [999]
        futures_df          = [999]
        


    merge_dates = data_merge['Dates']
    merge_times = data_merge['Times']
    merge_df    = data_merge['Df']

    #print'merge_times: ', len(merge_times) 
    
    swap_discount_factor = np.zeros(len(swap_dates))

    #%----------------------------------------------------------------------------------------------
    #%-------------- ELABORAZIONE SEGMENTO SWAP: FATTORI DI SCONTO A PARTIRE DAI TASSI SWAP ---------
    #%----------------------------------------------------------------------------------------------
    
    #'merge_df: ', merge_df
    #FQ(2233)
    
    if (nsw > 0):


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
                    print 'fwd: ',fwd
            
                else:

                    fwd         = -np.log(merge_df[len(merge_df)-1]/merge_df[len(merge_df)-2])/(merge_times[len(merge_times)-1] - merge_times[len(merge_times)-2])
                    df_tmp1s[i] = np.exp(-fwd*(times_swap_tmp1s[i]- merge_times[len(merge_times)-1]))
                    print 'fwd: ',fwd
    
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
    
        #print 'df_fix_swap: ', df_fix_swap
        #FQ(2233)
        
        for i in range(0,nsw2s):
            
            z_all[nsw1s + i] = compute_z_from_swap_rate(z_all, df_fix_swap, tenors, swap_val_new[i], tenors[nsw1s + i])
            
        #%--------------------------------  set LSR method OUTPUT --------------------------------

        #print 'z_all[nsw1s + i]: ', z_all[nsw1s + i]
        #FQ(2223)
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

            #print 'index_irregular_0: ', index_irregular_0
            
            #FQ(2223) 
            times_swap_new  = times_swap_new + fix_swap_time

            k = 0
            
            #print 'end_irregular_time: ', end_irregular_time_int
            #print '----------------------------------------------------'
            #print 'nsw1s: ', nsw1s 
            #print 'nsw1s + nsw2s: ', nsw1s + nsw2s 
            
            for i in range(nsw1s, nsw2s+nsw1s):     #% ciclo su tutte le date irregolari
            
                t_ref = times_swap_[i]
                
                if ((t_ref)  in  end_irregular_time_int):
                
                    #n_i = idx_start_irregular[k]  #% indice inizio data irregolare
                    m_i = idx_end_irregular[k]    #% indice fine data irregolare
                    
                    t_swt_n = start_irregular_time[k]
                    t_swt_m = end_irregular_time[k]

                    swp_tmp =  swap_val[m_i - 1]
                    
                    """
                    print 't_swt_n: ', t_swt_n
                    print 't_swt_m: ', t_swt_m
                    print 'swp_tmp: ', swp_tmp
                    print '------------------------------------'
                    """
                    
                    z_out = compute_z_from_cfr(z_out, tenors, swp_tmp, times_swap_new, t_swt_n, t_swt_m, times_swap_, df_fix_swap)
    
                    k = k + 1
                else:

                    
                    t_swt_n = t_ref
                    
                    nz  = find_indx_n(t_swt_n, times_swap_) #% indice del vettore times_out corrispondente a n anni

                    """
                    print 'i: ', i
                    print 'nz: ', nz
                    print 't_swt_n: ', t_swt_n
                    print '-----------------------------------'
                    """
                    
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
        
        
        print ''
        
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
            
        
    #%-------------------------------------------------------------------------
    #%-------------- GRAFICI DI CONTROLLO  ------------------------------------
    #%-------------------------------------------------------------------------
    
    merge_dates = data_merge['Dates'] 
    merge_times = data_merge['Times']
    merge_rates = data_merge['Rates']
    merge_df    = data_merge['Df']


    
    flag_make_graph = 0
    if (flag_make_graph == 1):


       
        g1 = graphrates(seg1_times[1:], seg1_values, futures_start_times, futures_rates, swap_times, swap_val, merge_times, merge_rates)
        g2 = graphdf(seg1_times, seg1_df, futures_times, futures_df, swap_times, swap_discount_factor)

        '--------------- save graph ------------------------------------------------------------'

    flag_save_graph = 0
    
    if (flag_save_graph == 1):
    
        
        g1.savefig('fig/plot_tassi_%s.png' %(ref_date))
        g2.savefig('fig/plot_discount_%s.png' %(ref_date))
 
    #%-------------------------------------------------------------------------
    #%-------------- SET OUTPUT  ROUTINES -------------------------------------
    #%-------------------------------------------------------------------------

    nodi_out = convertTimeToNode(merge_times)

    data_elab_out = {}
    
    data_elab_out['DateScadenza']    = merge_dates        
    data_elab_out['DiscountFactors'] = merge_df
    data_elab_out['TassiZC']         = merge_rates
    data_elab_out['Tempi']           = merge_times
    data_elab_out['Nodi']            = nodi_out
    data_elab_out['CodiceCurva']     = codeCurve
    
    return data_elab_out








def boot3s_elab_v2_old(data_opt, data_raw):

    #%--------------------------------------------------------------------------
    #%-------------- RETRIEVE DATA  -----------------------------------------
    #%--------------------------------------------------------------------------


    ref_field = 'UsaNodo'
    val_not_allowed = 'N'

    data_raw_p = purge_data(data_raw, ref_field, val_not_allowed)
    
    dict1s, dict_f, dict_s = select_segments(data_raw_p)
    
    #print 'dict1s: ', dict1s
    

    seg1_dates  = dict1s['MatDate']
    seg1_values = dict1s['ValoreNodo']
    
    seg1_values = seg1_values/100.0

    futures_start_dates = dict_f['MatDate']
    futures_values      = dict_f['ValoreNodo']
    
    swap_dates          = dict_s['MatDate']
    swap_val            = dict_s['ValoreNodo']

    swap_val = swap_val/100.0

    n1s  = len(seg1_dates)
    nf   = len(futures_start_dates)
    nsw  = len(swap_dates)
    

    #%--------------------------------------------------------------------------
    #%-------------- RETRIEVE OPTIONS SETUP -----------------------------------------
    #%--------------------------------------------------------------------------

    setting_default = set_data_default()

    ref_date = data_opt['RefDate'] 
    
    par_a = data_opt['ParConvexity']['A']
    par_b = data_opt['ParConvexity']['B']

    par_convexity    = [par_a, par_b]

    tenor_swap       = data_opt['TenorSwap']    
    tenor_swap       = convertNodeToMnth(tenor_swap)

    print '---------------------------------'
    print 'tenor_swap: ', tenor_swap

    
    flag_convexity   = int(data_opt['Convexity'])
    flag_futures_gap = int(data_opt['GapFutures'])
    flag_method_swap = int(data_opt['SwapGapMethod'])
    flag_interp1     = int(setting_default['InterpLinFutures'])

    flag_make_graph = data_opt['MakeGraph']
    flag_save_graph = setting_default['SaveGraph']
    regime_output   = int(data_opt['RegimeOutput'])   
    future_tenor    = data_opt['FutureTenor'] 
    
    '------- retrive conventions ------------------------'

    mkt_code = data_opt['MKT']
    mkt_ref  = holidays.get_calendar(mkt_code)
    
    basis_s = data_opt['Basis']['S']
    basis_f = data_opt['Basis']['F']

    basis_f = convert_basis(basis_f)
    basis_s = convert_basis(basis_s)

    day_conv_f = busdayrule.translateCodeFromDB(data_opt['BusConv']['F'])
    day_conv_s = busdayrule.translateCodeFromDB(data_opt['BusConv']['S'])
    day_conv_tn = busdayrule.translateCodeFromDB(setting_default['BusConv']['TN'])
    day_conv_on = busdayrule.translateCodeFromDB(setting_default['BusConv']['O/N'])
    
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
    
    if (len(futures_start_dates) > 0):


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
    
        
    
    if (len(futures_values) > 0):

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



    if (len(swap_dates)> 0):

        
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

    if (nf > 0):
  
        futures_rates = (100.0 - futures_values)/100.0  #-------- tasso future
    
        #-------- Correzione: "Convexity"------------------------------------------------'
        
        if (flag_convexity == 1):
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
            
        
        future_rates = compute_rates(futures_df, futures_times, regime_output)
    
        df_array = futures_df
        r_array  = future_rates
        t_array  = futures_times
        d_array  = futures_dates

        '------------------- merge ------------------'

        data_merge = merge_data(data_merge, d_array, t_array, r_array, df_array)
    
    else:

        futures_start_times = [999]
        futures_times       = [999]
        futures_rates       = [999]
        futures_df          = [999]
        


    merge_dates = data_merge['Dates']
    merge_times = data_merge['Times']
    merge_df    = data_merge['Df']

    print'len(merge_times): ', len(merge_times) 
    
    swap_discount_factor = np.zeros(len(swap_dates))

    #%----------------------------------------------------------------------------------------------
    #%-------------- ELABORAZIONE SEGMENTO SWAP: FATTORI DI SCONTO A PARTIRE DAI TASSI SWAP ---------
    #%----------------------------------------------------------------------------------------------
    
    
    if (nsw > 0):


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
        
        for i in range(0, len(idx_end_irregular)):
        
            n_i = idx_start_irregular[i]  #% indice inizio data irregolare
            m_i = idx_end_irregular[i]    #% indice fine data irregolare
            
            
            swap_times_n  = swap_times_x[n_i]  #% tempo inizio data irregolare
            swap_times_m  = swap_times_x[m_i]  #% tempo fine data irregolare
    
            t_swt_n = np.round(swap_times_n)   #% n. anni dell'n-ma data "irregolare"
            t_swt_m = np.round(swap_times_m)   #% m. anni dell'm-ma data "irregolare"

            start_irregular_time.append(t_swt_n)
            end_irregular_time.append(t_swt_m)
        
        start_irregular_time = np.asarray(start_irregular_time)
        end_irregular_time = np.asarray(end_irregular_time)
        
        index_irregular  = np.where(index_condition_irr)  # vettore indici date "irregolari" (array, dtype=int64)        
        index_irregular  = index_irregular[0]         # vettore indici date "irregolari"

        print 'index_irregular: ', index_irregular
        print 'idx_start_irregular: ', idx_start_irregular

        FQ(999)
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
        index_start = find_indx(swap_times[0], times_swap_new) # indice di inizio del segmento degli swap

        #times_swap_tmp1s  = times_swap_new[:index_start + 1]   # tempi fino al primo tasso swap
        #times_swap_tmp2s  = times_swap_new[index_start  + 1:]  # tempi oltre il primo tasso swap

        times_swap_tmp1s  = times_swap_new[:index_start ]   # tempi fino al primo tasso swap
        times_swap_tmp2s  = times_swap_new[index_start  :]  # tempi oltre il primo tasso swap

        #dates_swap_tmp1s  = swap_dates_new[:index_start + 1]   # date fino al primo tasso swap        
        #dates_swap_tmp2s  = swap_dates_new[index_start  + 1:]  # date oltre il primo tasso swap

        dates_swap_tmp1s  = swap_dates_new[:index_start ]   # date fino al primo tasso swap        
        dates_swap_tmp2s  = swap_dates_new[index_start  :]  # date oltre il primo tasso swap


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
    
        
        
        z_all        = np.concatenate((df_tmp1s, np.zeros(nsw2s)))
    
        #%-------------------------------------------------------------------------------------
        #%--------------- Generazione dei fattori di sconto alle date dei tassi swap ----------
        #%-------------------------------------------------------------------------------------
    
        for i in range(0,nsw2s):
            
            z_all[nsw1s + i] = compute_z_from_swap_rate(z_all, df_fix_swap, tenors, swap_val_new[i], tenors[nsw1s + i])
            
        #%--------------------------------  set LSR method OUTPUT --------------------------------
    
        if (flag_method_swap == 1): #% considero il metodo: Linear Swap Rate: 'LSR'

            swap_times_tmp = np.insert(times_swap_new, 0, 0)
            swap_times_tmp = np.round(swap_times_tmp, 2)
            z_out = z_all

            #%------------------------------------------------------------------------------------
            #%--------------------------------  CFR METHOD ----------------------------------------%
            #%------------------------------------------------------------------------------------
        
        elif (flag_method_swap == 0) and (ln_irr > 1): # Considero il metodo "Constant fwd rate" 'CFR'
            
            index_irregular_first = int(np.round(swap_times_x[index_irregular[0]]/(tenor_swap/12.0))) # % indice prima data irregolare
            index_irregular_last  = int(np.round(swap_times_x[index_irregular[len(index_irregular)-1]]/(tenor_swap/12.0))) #; % indice ultima data irregolare


            index_irregular_0     =  max(len(df_tmp1s)+1, index_irregular_first + 1)
            index_irregular_last  =  index_irregular_last

            z_out = np.zeros(max(1, index_irregular_last + 1))
                
            z_out[:index_irregular_0] = z_all[: index_irregular_0]

            times_swap_       = np.round(times_swap_x_new/(tenor_swap/12.0))*(tenor_swap/12.0)
             
                        
            times_swap_new  = times_swap_new + fix_swap_time

            k = 0
            for i in range(nsw1s, nsw2s+nsw1s):     #% ciclo su tutte le date irregolari
            
                t_ref = times_swap_[i]

                if (t_ref  in  end_irregular_time):
                
                    #n_i = idx_start_irregular[k]  #% indice inizio data irregolare
                    m_i = idx_end_irregular[k]    #% indice fine data irregolare
                    
                    
                    t_swt_n = start_irregular_time[k]
                    t_swt_m = end_irregular_time[k]

                    swp_tmp =  swap_val[m_i-1]

                    """
                    print 't_swt_n: ', t_swt_n
                    print 't_swt_m: ', t_swt_m
                    print 'swp_tmp: ', swp_tmp
                    """
                    
                    z_out = compute_z_from_cfr(z_out, tenors, swp_tmp, times_swap_new, t_swt_n, t_swt_m, times_swap_, df_fix_swap)

                    
                    """
                    dt_ref.append(t_swt_n)  
                    dt_ref.append(t_swt_m)  
                    dt_ref = np.asarray(dt_ref)
                    dz_ref = z_out[m_i-1:m_i+1]
                    r_out = compute_rates(dz_ref, dt_ref, regime_output)
                    print 'r_out: ', r_out[0] 
                    print 't_n: ', t_n
                    print '-----------------------------------'
                    """
                    
                    k = k + 1
                else:

                    
                    t_swt_n = t_ref
                    
                    nz  = find_indx_n(t_swt_n, times_swap_) #% indice del vettore times_out corrispondente a n anni

                    
                    print 'i: ', i
                    print 'nz: ', nz
                    print 't_swt_n: ', t_swt_n
                    print '-----------------------------------'
                
                    
                    swp_tmp =  swap_val_new[i-nsw1s]
                    
                    z_out_n  = z_out[:nz]
                    tenors_n = tenors[:nz]

                    z_out[nz] = compute_z_from_swap_rate(z_out_n, df_fix_swap, tenors_n, swp_tmp, tenors[i])

                    
                    
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
                    
        else:
            
            
            z_out = z_all[1:]
            times_swap_       = np.round(times_swap_x_new/(tenor_swap/12.0))*(tenor_swap/12.0)            


            
        #%---------------------------------------------------------------------
        #%------------------ set swap output ----------------------------------
        #%---------------------------------------------------------------------

        times_swap_       = np.round(times_swap_x_new/(tenor_swap/12.0))*(tenor_swap/12.0)            

        swp__ = np.round(swap_times_x/(tenor_swap/12.0), 0)
        swap_times_adj = np.round(swp__*(tenor_swap/12.0),1)
        
        if (swap_times_adj[0] == 0): swap_times_adj = swap_times_adj[1:]
        
        for i in range(0, nsw):
            
            #indxTmp = find_indx_n(swap_times_adj[i], times_swap_[1:])
            indxTmp = find_indx_n(swap_times_adj[i], times_swap_[1:])            
            swap_discount_factor[i] = z_out[indxTmp+1]
            
            """
            try:
                swap_discount_factor[i] = z_out[indxTmp + 1]
            except:
                swap_discount_factor[i] = z_out[indxTmp]
            """
        swap_rates_out = compute_rates(swap_discount_factor, swap_times_adj, regime_output)


        df_array = swap_discount_factor
        r_array  = swap_rates_out
        t_array  = swap_times
        d_array  = swap_dates

        data_merge = merge_data(data_merge, d_array, t_array, r_array, df_array)

    else:
        
        swap_times = [999]
        swap_val   = [999]
            
        
    #%-------------------------------------------------------------------------
    #%-------------- GRAFICI DI CONTROLLO  ------------------------------------
    #%-------------------------------------------------------------------------
    
    merge_dates = data_merge['Dates'] 
    merge_times = data_merge['Times']
    merge_rates = data_merge['Rates']
    merge_df    = data_merge['Df']

    
    flag_make_graph = 0
    
    if (flag_make_graph == 1):

    
        #seg1_values = seg1_values/100.0
        #futures_rates = np.asarray(futures_rates)/100.0
        #swap_val = swap_val/100.0
        merge_rates = merge_rates
       
        seg1_values_to_graph = seg1_values*100.0
        futures_rates_to_graph = futures_rates*100.0
        swap_val_to_graph = swap_val*100.0
       
       
        g1 = graphrates(seg1_times[1:], seg1_values_to_graph, futures_start_times, futures_rates_to_graph, swap_times, swap_val_to_graph, merge_times, merge_rates)
        g2 = graphdf(seg1_times, seg1_df, futures_times, futures_df, swap_times, swap_discount_factor)

        '--------------- save graph ------------------------------------------------------------'

    
    flag_save_graph = 0
    
    if (flag_save_graph == 1):
    
        
        g1.savefig('fig/plot_tassi_%s.png' %(ref_date))
        g2.savefig('fig/plot_discount_%s.png' %(ref_date))
 
    #%-------------------------------------------------------------------------
    #%-------------- SET OUTPUT  ROUTINES -------------------------------------
    #%-------------------------------------------------------------------------

    nodi_out = convertTimeToNode(merge_times)

    data_elab_out = {}
    
    data_elab_out['DateScadenza']    = merge_dates        
    data_elab_out['DiscountFactors'] = merge_df
    data_elab_out['TassiZC']         = merge_rates
    data_elab_out['Tempi']           = merge_times
    data_elab_out['Nodi']            = nodi_out
    
    return data_elab_out







  

def boot3s_elab_n(data_opt, data_raw):
    
    
    ref_field = 'UsaNodo'
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

    nodo1s = dict1s['ValoreNodo']
    nodo2s = dict_f['ValoreNodo']
    nodo3s = dict_s['ValoreNodo']
    

    basis1s = data_opt['Basis']['D']    
    basis2s = data_opt['Basis']['F']
    basis3s = data_opt['Basis']['S']


    regime_output   = data_opt['RegimeOutput']
    
    if (len(dates1s)>0): dates1s.insert(0, refDate)
    if (len(rate1s)>0): np.insert(rate1s, 0, rate1s[0])
    if (len(dates2s)>0): dates2s.insert(0, refDate)
    if (len(dates3s)>0): dates3s.insert(0, refDate)
    
    times1s = daycount.yearfractions(dates1s, basis1s);
    times2s = daycount.yearfractions(dates2s, basis2s);
    times3s = daycount.yearfractions(dates3s, basis3s);
    
    if (len(times2s)>0): times2s = times2s[1:]
    if (len(times3s)>0): times3s = times3s[1:]

    
    
    times1s = times1s[1:]
    
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

    rates1s = np.insert(rates1s, 0, rates1s[0])
    times1s = np.insert(times1s, 0, times1s[0])
    df1s = np.insert(df1s, 0, df1s[0])
    nodo1s = np.insert(nodo1s, 0, nodo1s[0])
    
    print 'len(nodo1s): ', len(nodo1s)
    print 'len(dates1s): ', len(dates1s)
    print 'len(df1s): ', len(df1s)
    print 'len(times1s): ', len(times1s)
    print 'len(dates1s): ', len(dates1s)
    
    

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

    
    #print 'data_elab_out: ', data_elab_out
    
    #FQ(223)
    
    return data_elab_out

def convert_basis(basis_to_convert):
    
    validated_basis = {}

    validated_basis['UNDEF']             = 'ACT/365'    
    validated_basis['30/365']            = 'ACT/365'
    validated_basis['ACT/nACT']          = 'ACT/365'
    validated_basis['ACT/365.25']        = 'ACT/365'
    validated_basis['(ACT+1)/360']       = 'ACT/365'
    validated_basis['(ACT+1)/365']       = 'ACT/365'
    validated_basis['ACT/365 (366)']     = 'ACT/365'
    validated_basis['30E+1/360 (ITL)']   = 'ACT/365'
    
    key_list = validated_basis.keys()

    if (basis_to_convert in key_list):
        
        converted_basis = validated_basis[basis_to_convert]
        
        
    else:
        converted_basis = basis_to_convert

    return converted_basis



def set_data_default():
    
    data_default = {}

    data_default['SaveGraph'] = 0
    
    data_default['BasisFix'] = 'ACT/ACT'
    data_default['BasisRef'] = 'ACT/ACT'

    data_default['RegimeRate'] ={}
    data_default['RegimeRate']['D'] = 0 # 0 = Semplice, 1 = Composto, 2 = Continuo
    data_default['RegimeRate']['L'] = 0
    
    data_default['InterpLinFutures'] = 1
    
    data_default['BusConv'] = {}
    data_default['BusConv']['TN'] = 'follow'
    data_default['BusConv']['O/N'] = 'follow'

    
    return data_default


def set_data_opt():
    
    
    
    data_opt = {}
    data_opt['Basis'] = {}
    data_opt['Basis']['D']  = 'ACT/360'
    data_opt['Basis']['L']  = 'ACT/360'
    data_opt['Basis']['S']  = 'ACT/360'
    data_opt['Basis']['F']  = 'ACT/365'
    data_opt['Basis']['TN'] = 'ACT/360'
    data_opt['Basis']['O/N'] = 'ACT/360'
    
    data_opt['BusConv'] = {}
    data_opt['ParConvexity'] = {}
    data_opt['ParConvexity']['A'] = 0.0
    data_opt['ParConvexity']['B'] = 0.0

    data_opt['FutureTenor'] = 90

    data_opt['TenorSwap'] = '6M'
    data_opt['Convexity'] = 1 # 1 = attivo, 0 = inattivo
    data_opt['GapFutures'] = 1 # 1 = previous spot rate, 0 = next forward rate
    data_opt['SwapGapMethod'] = 0 # 0 = Constant Forward Rate, 1 = Linear Swap Rate
    data_opt['InterpLinFutures'] = 1 # 1 = interp. esponenziale sui fattori di sconto, 0 = interp. line. sui tassi forward

    data_opt['RefDate'] = datetime.date(2017, 02, 15)

    data_opt['BusConv']['D'] = 'follow'
    data_opt['BusConv']['L'] = 'follow'
    data_opt['BusConv']['S'] = 'follow'
    data_opt['BusConv']['F'] = 'follow'
    data_opt['BusConv']['TN'] = 'follow'
    data_opt['BusConv']['O/N'] = 'follow'
    
    data_opt['MKT']          = 'de'


    data_opt['MakeGraph'] = True
    data_opt['SaveGraph'] = True
    data_opt['RegimeOutput'] = '2' # 0 = Semplice, 1 = Composto, 2 = Continuo


    return data_opt





