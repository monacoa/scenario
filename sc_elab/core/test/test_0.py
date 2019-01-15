


import sys
import datetime as dtime
import math
import numpy as np

#from sc_elab import numpy as np
from mdates import holidays
from mdates import daycount
from mdates import busdayrule
from scipy import optimize
import matplotlib.pyplot as plt

import datetime
import calendar

def FQ(label):
    print ('------------- FIN QUI TUTTO OK  %s ----------' %(label))
    sys.exit()



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


def map_segment(data_raw):    
                
    list_s = set(data_raw['TipoSegmento'])
    list_s = list(list_s)
    
    print 'list_s: ', list(list_s)  
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
        len(spot_rate)
        df_out = 1.0/(1.0 + spot_rate*time_ref);

    except:
        
        df_out = np.ndarray(len(time_ref))
        for i in range(0, len(time_ref)):

            df_tmp = 1.0/(1.0 + float(spot_rate[i])*float(time_ref[i]));
            df_out.append(df_tmp)

    return df_out
    
def df_cont(spot_rate, time_ref):

    try:
        len(spot_rate)
        df_out = math.exp(-time_ref*spot_rate)

    except:
        
        df_out = np.ndarray(len(time_ref))
        for i in range(0, len(time_ref)):

            df_tmp = math.exp(-float(time_ref[i]*spot_rate[i]))
            df_out.append(df_tmp)

    return df_out

def df_cmp(spot_rate, time_ref):

    try:
        len(spot_rate)
        df_out = 1.0/(1 + spot_rate)**(time_ref)

    except:
        
        df_out = np.ndarray(len(time_ref))
        for i in range(0, len(time_ref)):

            df_tmp = 1.0/(1.0 + spot_rate[i])**time_ref[i]
            df_out.append(df_tmp)

    return df_out


def find_indx(x_target, x_list_ref):

    ln = len(x_list_ref) - 1

    for i in range(0, len(x_list_ref)):

        if (x_list_ref[0] > x_target): 
            i_ref = None
            break 

        elif (x_list_ref[ln] <= x_target): 
            i_ref = None
            break
        
        elif (x_list_ref[i] >= x_target):
            i_ref = i - 1
            break
            
        else:

            continue

    return i_ref


def find_indx_p(ref_cond, list_condition):

    i_out = 0

    ln = len(list_condition)
    for i in range(0, ln):
        
        valueTmp = list_condition[i]
        
        if (valueTmp == ref_cond):
        
            i_out = i
            break
        else:
            continue

    return i_out



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


#def compute_df_at_first_future_date(seg1_times, seg1_val, seg1_df,  futures_rates[0], flag_futures_gap, futures_start_time[0], futures_end_time, flag_interp1):
def compute_df_future(seg1_times, seg1_val, seg1_df,  futures_rates_target, flag_futures_gap, futures_start_time_n, futures_end_time_n, flag_interp1):
    
    
    index_seg1_last = find_indx(futures_start_time_n, seg1_times)
    index_cross     = find_indx(futures_start_time_n, seg1_times)
    

    
    #index_seg1_last = find(seg1_times <= futures_start_time[0], 1, 'last');
    #index_cross = find(seg1_times == futures_start_time[0]);
    
    if not(index_cross != None):
        futures_start_1 = seg1_df[index_cross]
        futures_start_df_0 = 0 
        futures_end_df_0= futures_start_1/(1.0 + futures_rates_target*(futures_end_time_n - futures_start_time_n) );

    elif (seg1_times[len(seg1_times)] > futures_start_time_n): #% controllo della sovrapposizione: seg1 e futures

        t_target = futures_start_time_n
        df_n = seg1_times[index_seg1_last + 1] 
        df_o = seg1_times[index_seg1_last]
        t_n = seg1_times[index_seg1_last + 1]
        t_o = seg1_times[index_seg1_last]

        if (flag_interp1 == 1): #%----- SI SOVRAPPOSIZIONE
            
            #%----------intepolazione esponenziale sui fattori di scont ---------------
            futures_start_df_0 =  df_from_interp_df_exp(t_target, df_n, df_o, t_n, t_o)


        else: #%-------------- interpolazione lineare sui tassi -----------------------
            futures_start_df_0 = df_from_interp_R_lin(t_target, df_n, df_o, t_n, t_o)
            
        futures_end_df_0= futures_start_df_0/(1.0 + futures_rates_target*(futures_end_time_n - futures_start_time_n))

    else:  # ---- NO SOVRAPPOSIZIONE

        if any(seg1_val == 0):
            R_start1 = futures_rates_target
            futures_start_df_0 = np.exp(-R_start1*futures_start_time_n)
        
        elif (flag_futures_gap == 1) and (seg1_val != 0): #------------ gap colmato mantenedo costante il tasso spot
           
            R_start1 = -np.log(seg1_df(index_seg1_last))/seg1_times(index_seg1_last)
            futures_start_df_0 = np.exp(-R_start1*futures_start_time_n)
            
        else:  #------------ gap colmato mantenedo costante il tasso fwd ---------------------------------------------
            
            fwd_start1 = -np.log( seg1_df[index_seg1_last]/seg1_df[index_seg1_last - 1] )/( seg1_times[index_seg1_last] - seg1_times[index_seg1_last - 1])
            futures_start_df_0 = seg1_df[index_seg1_last]*np.exp(-fwd_start1*( futures_start_time_n- seg1_times[index_seg1_last] ))
            
        futures_end_df_0 = futures_start_df_0/(1 + futures_rates_target*(futures_end_time_n - futures_start_time_n) )

    return futures_start_df_0, futures_end_df_0 



def compute_convexity(par_convexity, futures_start_tx, futures_end_tx):
    
    nf = len(futures_start_tx) 
    conv_correction  = np.ndarray(len(futures_start_tx))
    mean_convexity   = par_convexity[0]
    vol_convexity    = par_convexity[1]
        
    for i in range(0, nf):

        t1 = futures_start_tx[i]
        t2 = futures_end_tx[i]

        if (mean_convexity == 0): #---------limite a-->0 (Hull& White) --> (Ho Lee)
            conv_corrTmp = 0.5*t1*t2*vol_convexity**2
        else:
            #%--------------- Calcolo della correzione mediante il modello di Hull&White
            
            B1  = ( 1 - np.exp(- mean_convexity*(t2 - t1)) )/mean_convexity
            B2  = ( 1 - np.exp(- mean_convexity*(t1) ) )/mean_convexity
            conv_corrTmp  = B1/(t2 - t1)*( B1*(1 - np.exp(-2*mean_convexity*t1)) + 2*mean_convexity*B2**2)*vol_convexity**2/(4*mean_convexity)

        conv_correction.append(conv_corrTmp)
            
            
    return conv_correction

def compute_rates(df_in, time_df_in, regime_output):
    
    
    n_dmerge = len(df_in)
    out_rates = np.zeros(n_dmerge)
    
    if (regime_output == 1):     #% regime composto
        #regime = 'composto';
        for i in range(0, n_dmerge):
            out_r_tmp    = ( (1.0/df_in[i])**(1.0/time_df_in[i]) - 1.0) #% regime composto
            out_rates.append(out_r_tmp)

    elif (regime_output == 0): #% regime semplice
        #regime = 'semplice';
        for i in range(0, n_dmerge):
            out_r_tmp    = ( (1.0/df_in[i]) - 1.0 )/time_df_in[i] #% semplice
            out_rates.append(out_r_tmp)

    elif (regime_output == 2): #% regime continuo
        #regime = 'continuo';
        for i in range(0, n_dmerge):
            out_r_tmp    = - np.log(df_in[i])/time_df_in[i] #% semplice
            out_rates.append(out_r_tmp)
    
    return out_rates 


def compute_first_df(seg1_dates, seg1_times, seg1_val, fix_time1, on_date, tn_date): 


    seg1_df_0 = df_simple(seg1_val[0], seg1_times[0])
    if ( (seg1_dates[0] != on_date) and (seg1_dates[0] != tn_date) ):

        spot_seg1_rate_1 = -math.log(seg1_df_0)/seg1_times[0];
        df_fix_seg1 = math.exp(-spot_seg1_rate_1*fix_time1);
        id_start = 1

    #---------------- caso 2) si ON, si TN ----------------------------------------------------
    elif (seg1_dates[0] == on_date) and (seg1_dates[1] == tn_date):

        seg1_df_1 = seg1_df_0*df_simple(seg1_val[1], seg1_times[1] - seg1_times[0])
        df_fix_seg1 = seg1_df_1
        id_start = 2

    #---------------- caso 3) no ON, si TN ---------------------------------------------------
    elif (seg1_dates[0] != on_date) and (seg1_dates[1] == tn_date):

        df_fix_seg1 = seg1_df_0
        id_start = 1

    #---------------- caso 4) si ON, no TN ---------------------------------------------------
    elif (seg1_dates[0] == on_date) and (seg1_dates[1]!= tn_date):

        #seg1_df(1) = 1./( 1 + seg1_val(1).*seg1_times(1) );
        seg1_df_1 = seg1_df_0*df_simple(seg1_val[1], seg1_times[1] - seg1_times[0])
        A = seg1_df_0;
        B = 1.0/(seg1_times[0] - seg1_times[1])*math.log(seg1_df_1/seg1_df_0)
        df_fix_seg1 = A*math.exp(-B*(fix_time1 - seg1_times[0]))
        id_start = 2
        
    else:
        
        df_fix_seg1 = seg1_df_0
        id_start = 1

    return df_fix_seg1, id_start 


def fwd_fun(x,n,m,f1, py_rate_m, time_new):

    #%--------------------------------------------------------------------------
    f0 = 0
    for k in range(n,m):    
        
        f0 = (time_new[k] - time_new[k-1])*np.exp(-x*(time_new[k]-time_new[n])) + f0
    
    f = py_rate_m*f0 + np.exp(-x*(time_new[m]-time_new[n])) - f1
    
    #%-------------------------------------------------------------------------
    return f

def make_graph_for_chk():


    """
    figure1 = figure('Name',['Tassi nel regime ', regime]);
    
    plot(seg1_times_g, seg1_rates_g, '-o')
    hold on
    plot(futures_times_g, futures_rates_g, 'r-*')
    hold on
    plot(swap_times_g, swap_rates_g, 'k-o')
    
    xlabel('Tempi [anni]')
    ylabel('Livello tassi')
    
    if     (  not(any(seg1_val == 0)) &&    not(any(futures_val == 0))  &&  not(any(swap_val == 0)) );
        
        if flag_1s == 1
        legend('Depositi', 'Futures', 'Swap', 0);
        elseif flag_1s == 2
        legend('Libor', 'Futures', 'Swap', 0);
        end
           
    elseif (  any(seg1_val == 0)      &&    not(any(futures_val == 0))  &&  not(any(swap_val == 0)) );
        legend('Futures', 'Swap', 0)
    elseif (  any(seg1_val == 0)      &&    any(futures_val == 0)       &&  not(any(swap_val == 0)) );
        legend('Swap', 0)
    elseif (  not(any(seg1_val == 0)) &&    any(futures_val == 0)       &&  any(swap_val == 0)      ); 
           
        if flag_1s == 1
        legend('Depositi', 0)
        elseif flag_1s == 2
        legend('Libor', 0);
        end
           
    elseif (  any(seg1_val == 0)      &&    not(any(futures_val == 0))  &&  any(swap_val == 0)      );
        legend('Futures', 0)
    elseif (not(any(seg1_val == 0))   &&    not(any(futures_val == 0))  && any(swap_val == 0)       );
    
        if flag_1s == 1
        legend('Depositi', 'Futures', 0)
        elseif flag_1s == 2
        legend('Libor', 'Futures', 0)
        end
    
    elseif (  not(any(seg1_val == 0)) &&    any(futures_val == 0)       && not(any(swap_val == 0)) );
     
        if flag_1s == 1
        legend('Depositi', 'Swap', 0)
        elseif flag_1s == 2
        legend('Libor', 'Swap', 0)
        end    
    end
    
    % legend('Depositi', 'Futures', 'Swap')
    
    title(['Tassi nel regime ', regime]);
    """
    

    g1 = plt.figure(1)
    return g1


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
        

def set_data_opt(): 
    
    data_opt['Basis'] = {}
    data_opt['Basis']['D'] = 'ACT/360'
    data_opt['Basis']['L'] = 'ACT/360'
    data_opt['Basis']['S'] = 'ACT/360'
    data_opt['Basis']['F'] = 'ACT/360'
    
    data_opt['BusConv'] = {}
    data_opt['ParConvexity'] = {}
    data_opt['ParConvexity']['A'] = 1.0
    data_opt['ParConvexity']['B'] = 0.25

    data_opt['TenorSwap'] = 6
    data_opt['Convexity'] = True
    data_opt['GapFutures'] = True
    data_opt['SwapGapMethod'] = True
    data_opt['InterpLinFutures'] = True

    data_opt['BusConv']['D'] = 'follow'
    data_opt['BusConv']['L'] = 'follow'
    data_opt['BusConv']['S'] = 'follow'
    data_opt['BusConv']['F'] = 'follow'
    data_opt['MKT']          = 'de'

    return data_opt

def boot3s_elab(data_opt, data_raw):


    '--------------------------------------------------------------------------------------'
    '------- retrive data -----------------------------------------------------------------'
    '--------------------------------------------------------------------------------------'


    futures_start_dates =[]
    swap_dates = []
    swap_val = []
    
    ref_date = swap_dates[0] # cambiare yearfrac perche punti alla prima data = ref_date
    seg1_dates = []
    seg1_times = []
    seg1_val   = []
    flag_regime_rate = 1
    futures_val = np.ndarray(10)


    par_a = data_opt['ParConvexity']['A']
    par_b = data_opt['ParConvexity']['B']

    par_convexity = [par_a,par_b]
    flag_convexity =  data_opt['Convexity']
    flag_futures_gap = data_opt['GapFutures']
    flag_interp1 = data_opt['InterpLinFutures']
    tenor_swap = data_opt['TenorSwap']
    flag_method_swap = data_opt['SwapGapMethod']
    
    nd  = len(seg1_dates);
    nf  = len(futures_start_dates);
    nsw = len(swap_dates);

    flag_1s = []
    flag_1s.append(1)
    flag_seg1          = flag_1s*np.ones(nd,1);
    
    #flag_futures       = 3*np.ones(2*nf,1);
    #flag_swap          = 4*np.ones(nsw,1);

    '--------------------------------------------------------------------------------------'
    '------- retrive data -----------------------------------------------------------------'
    '--------------------------------------------------------------------------------------'

    mkt_code = data_opt['MKT']
    mkt_ref  = holidays.get_calendar(mkt_code)
    
    map_s, n_s  = map_segment(data_raw)
    
    s1_code = map_s['1']
    s2_code = map_s['2']

    basis_1    = data_opt['Basis'][s1_code]
    basis_2    = data_opt['Basis'][s2_code]

    day_conv_1 = data_opt['BusConv'][s1_code]
    day_conv_2 = data_opt['BusConv'][s2_code]
    
    if n_s == 3:
         
        s3_code    = map_s['3']
        basis_3    = data_opt['Basis'][s3_code]
        day_conv_3 = data_opt['BusConv'][s3_code]
        
    
    basis_d = data_opt['Basis']['D']
    basis_l = data_opt['Basis']['L']
    basis_s = data_opt['Basis']['S']
    basis_f = data_opt['Basis']['F']

    day_conv_d = data_opt['BusConv']['D']
    day_conv_l = data_opt['BusConv']['L']
    day_conv_s = data_opt['BusConv']['S']
    day_conv_f = data_opt['BusConv']['F']
    
    dates_d, values_d = retrive_segment(data_raw, 'D')
    dates_l, values_l = retrive_segment(data_raw, 'L')
    dates_f, values_f = retrive_segment(data_raw, 'F')
    dates_s, values_s = retrive_segment(data_raw, 'S')
    
    basis_fix = basis_d
    basis_f_0 = basis_d
    basis_x = basis_s
    
    flag_make_graph = 1
    flag_save_graph = 1
    
    regime_output = 1
    
    '--------------------------------------------------------------------------'
    ' ------------------ inizializzazione data output -------------------------'
    '--------------------------------------------------------------------------'
    
    dtype  = [('Times', float), ('Dates', datetime.date), ('Df', float), ('Rates', float)]
    data_0 = np.array([], dtype=dtype)

    '-----------------------------------------------------------------------------------------'
    '------- generazione date di scadenza dei futures ----------------------------------------'
    '-----------------------------------------------------------------------------------------'
    
    
    future_tenor = 3*30
    
    if any(futures_start_dates):

        futures_end_dates  = busdayrule.rolldates(futures_start_dates + future_tenor, mkt_ref, day_conv_f)


    '-----------------------------------------------------------------------------------------'
    '------generazione dei tempi a partire dalle date ----------------------------------------'
    '-----------------------------------------------------------------------------------------'

    # ------------------ da rivedere ------------------------------------------
    if not(any(futures_val)):
        futures_end_times    = daycount.yearfractions(futures_end_dates, basis_f)
        futures_start_times  = daycount.yearfractions(futures_start_dates, basis_f)
        
        futures_start_tx = daycount.yearfractions(futures_start_dates, basis_f_0)
        futures_end_tx = daycount.yearfractions(futures_start_dates, basis_f_0)

    if not(any(swap_dates)):
        swap_times           = daycount.yearfractions(swap_dates, basis_s)
        swap_times_x         = daycount.yearfractions(swap_dates, basis_x)


    '-----------------------------------------------------------------------------------------'
    '------ generazione fattori di sconto per dati OverNight, Tomorrow Next ------------------'
    '-----------------------------------------------------------------------------------------'

    on_date = busdayrule.rolldates(ref_date + 1, mkt_ref, day_conv_s)
    tn_date = busdayrule.rolldates(on_date  + 1, mkt_ref, day_conv_s) 
    
    fix_date1      = tn_date;
    fix_time1      = daycount.yearfrac(ref_date, fix_date1, basis_fix)
    fix_swap_dates  = busdayrule.rolldates(on_date  + 3,  mkt_ref, day_conv_s) # nel codice matlab e' implementato in modo un po' differente
    fix_swap_times = daycount.yearfrac(ref_date, fix_swap_dates, basis_s)
    
    #------------------------------------------------------------------
    #e = ref_date;
    
    #for i = 1:n_fixing_days(3)
    #    e = busdate(e, 1, bank_holidays);
    #end
    #------------------------------------------------------------------
    
    '-----------------------------------------------------------------------------------------'
    '------ generazione fattori di sconto alle date ON/TN7FIX --------------------------------'
    '-----------------------------------------------------------------------------------------'

    if not(any(seg1_val == 0)):
            
        seg1_df = np.ndarray(len(seg1_val))
        df_fix_seg1, id_start = compute_first_df(seg1_dates, seg1_times, seg1_val, fix_time1, on_date, tn_date)
      
        seg1_df[0:id_start] = df_fix_seg1
      
        '-----------------------------------------------------------------------------------------'
        '------ generazione fattori di sconto depositi, Libor ------------------------------------'
        '-----------------------------------------------------------------------------------------'
    
        seg1_times_m_fix = seg1_times[id_start:] - fix_time1
        seg1_times_m_fix = seg1_times[id_start:] - fix_time1
      
        if (flag_regime_rate == 0): #% regime semplice
    
            seg1_df_ = df_fix_seg1*df_simple(seg1_val[id_start:], seg1_times_m_fix)
    
        elif (flag_regime_rate == 1): #% regime composto
    
            seg1_df_ = df_fix_seg1*df_cmp(seg1_val[id_start:], seg1_times_m_fix)
    
        elif (flag_regime_rate == 2): #% regime composto
    
                seg1_df_ = df_fix_seg1*df_cont(seg1_val[id_start:], seg1_times_m_fix)
        else:
            
            print 'Regime non gestito'
        
        seg1_df[id_start:] = seg1_df_
            
        
        
        #-------------------------------------------------------------    
        # ----------- set output -------------------------------------
        #-------------------------------------------------------------
        
        seg1_times   = [0, seg1_times]
        seg1_dates   = [ref_date, seg1_dates]
        seg1_df      = [1, seg1_df]
        flag_seg1    = [flag_1s[0], flag_seg1]
    
 
    else:

        df_fix_seg1 = 1
        seg1_times = 0
        seg1_dates = ref_date
        seg1_df    = 1

    '-----------------------------------------------------------------------------------------'
    '------ generazione fattori di sconto per i futures --------------------------------------'
    '-----------------------------------------------------------------------------------------'

    if not(any(futures_val == 0)):
  
        futures_rates = (100.0 - futures_val)/100.0 #-------- tasso future
    
        #-------- Correzione: "Convexity"------------------------------------------------
        
        if (flag_convexity == 1):
            conv_correction = compute_convexity(par_convexity, futures_start_tx, futures_end_tx)
            futures_rates = futures_rates - conv_correction;

        #----------------------------------------------------------------------------------------
        #------------------ generazione di Z alla prima data "start" del futures -----------------
        #-----------------------------------------------------------------------------------------
    
        nf = len(futures_end_times)
    
        futures_start_df = np.ndarray(nf)
        futures_end_df   = np.ndarray(nf)
    
        for i in range(0, nf):
            
            # ------------ SISTEMARE LE PRIME DATE -----------------
            futures_end_time_n   = futures_end_times[i]
            futures_start_time_n = futures_start_times[i]
            futures_rates_target = futures_rates[i]
        
            futures_start_df_tmp, futures_end_df_tmp = compute_df_future(seg1_times, 
                                                                       seg1_val, 
                                                                       seg1_df,  
                                                                       futures_rates_target, 
                                                                       flag_futures_gap, 
                                                                       futures_start_time_n,
                                                                       futures_end_time_n, 
                                                                       flag_interp1)
    
            futures_start_df[i] =  futures_start_df_tmp
            futures_end_df[i]   =  futures_end_df_tmp
        
        '------------------ SET OUTPUT  ----------------------'
    
    
        #index_zeros = find(futures_start_df == 0);
        #if not(isempty(index_zeros))
        #    futures_start_df_zeros    = futures_start_df(index_zeros);
        #    futures_start_time_zeros  = futures_start_time(index_zeros);
        #    futures_start_dates_zeros = futures_start_dates(index_zeros);
            
        #   futures_start_df = setxor(futures_start_df_zeros, futures_start_df);
        #   futures_start_df = sort(futures_start_df, 'descend');
    
        #   futures_start_time = setxor(futures_start_time, futures_start_time_zeros);
        #   futures_start_dates = setxor(futures_start_dates, futures_start_dates_zeros);        
        #   flag_futures = flag_futures(1:end - length(index_zeros));
        #end
       
        
    else:

        futures_start_time  = []
        futures_end_time    = []
        futures_start_dates = []
        futures_end_dates   = []
        futures_end_df      = []
        futures_start_df    = []
        flag_futures        = []
    
    #------------------- merge --------------------------------------------------------------
    #--------------------------------------------------------------

    futures_start_rates = futures_start_df
    futures_end_rates = futures_end_df

    df_array = futures_start_df
    r_array  = futures_start_rates
    t_array  = futures_start_times
    d_array  = futures_start_dates
    
    data_merge = merge_data(data_0, d_array, t_array, r_array, df_array)

    df_array = futures_end_df
    r_array  = futures_end_rates
    t_array  = futures_end_times
    d_array  = futures_end_dates

    data_merge = merge_data(data_merge, d_array, t_array, r_array, df_array)

    """
    [union_time_futures, id_start, id_end]  = np.union1d(futures_start_time, futures_end_time);
    
    futures_start_dates_x = futures_start_dates[id_start]
    futures_end_dates_x   = futures_end_dates[id_end]
    futures_start_time_x  = futures_start_time[id_start]
    futures_end_time_y    = futures_end_time[id_end]

    futures_start_df_x    = futures_start_df[id_start]
    futures_end_df_x      = futures_end_df[id_end]
    flag_futures_x        = flag_futures[id_end]
    flag_futures_y        = flag_futures[id_start]
    flag_futures          = [flag_futures_x, flag_futures_y]
   

    merge_dates_futures_tmp = [futures_start_dates_x, futures_end_dates_x]
    merge_df_futures_tmp    = [futures_start_df_x, futures_end_df_x]
    merge_time_futures_tmp  = [futures_start_time_x, futures_end_time_y]
    merge_futures_data_tmp  = [merge_time_futures_tmp, merge_dates_futures_tmp, merge_df_futures_tmp]

    if (merge_futures_data_tmp.size == 0):
        futures_time  = []
        futures_dates = []
        futures_df    = []
    else:

        #sort_future_data   = sortrows(merge_futures_data_tmp, 1)
        sort_future_data   = merge_futures_data_tmp
        
        futures_time       = sort_future_data[1][:]
        futures_dates      = sort_future_data[2][:]
        futures_df         = sort_future_data[3][:]

        merge_flag_tmp  = [flag_seg1, flag_futures]
        merge_time_tmp  = [seg1_times, futures_time]
        merge_dates_tmp = [seg1_dates, futures_dates]
        merge_df_tmp    = [seg1_df, futures_df]
        merge_data_tmp  = [merge_time_tmp, merge_dates_tmp, merge_df_tmp, merge_flag_tmp]
    
        #sort_data     = sortrows(merge_data_tmp, 1)
        sort_data     = merge_data_tmp

        merge_times   = sort_data[:,1]
        merge_dates   = sort_data[:,2]
        merge_df      = sort_data[:,3]
        merge_flag    = sort_data[:,4]
    

    merge_dates = []
    merge_df = []
    merge_times = []
    
    """

    merge_dates = data_merge['Dates']
    merge_times = data_merge['Times']
    merge_rates = data_merge['Rates']
    merge_df    = data_merge['Df']

    '-----------------------------------------------------------------------------------------'
    '------ generazione dei fattori di sconto a partire dai tassi swap  ----------------------'
    '-----------------------------------------------------------------------------------------'
    
    if not(any(swap_val == 0)):

        if (flag_future == False):
            df_fix_swap = df_fix_seg1
            
        else:
            
            index_last = find_indx(fix_swap_dates, merge_dates)
            A =  merge_df[index_last];
            B = 1./(merge_times[index_last] - merge_times[index_last+1])*np.log(merge_df[index_last + 1]/merge_df[index_last])
            df_fix_swap = A*np.exp(-B*(fix_swap_times - merge_times[index_last] ))
        #end
    
        #%--------------------------- Def. variabili ----------------------------------------------------------
    
        dt_time_swap     = swap_times[2:]- swap_times[1:len(swap_times)-2]  # differenze tempi
        index_condition  = (round(dt_time_swap) > tenor_swap/12)            # vettore indici date "irregolari"
        index_irregular  = index_condition.index[False]                     # vettore indici date "irregolari"

        index_irregular  = [index_irregular, index_irregular[len(index_irregular)-1] +  1] # % indici irregolari
    
        n_new_swap       = round(swap_times_x[len(swap_times_x)-1]/(tenor_swap/12));
        swap_dates_new   = np.zeros(0, n_new_swap);
        times_swap_new   = np.zeros(0, n_new_swap);
        times_swap_x_new = np.zeros(0, n_new_swap);
    
        for j in range(0, n_new_swap):
        
            swap_dates_tmp = add_months(fix_swap_dates, tenor_swap*j)
            swap_dates_tmp = busdayrule.rolldates(swap_dates_tmp, mkt_ref, day_conv_s) #% genero vettore delle date corrette DAY CONVENTION
            swap_dates_new.append(swap_dates_tmp)
            
            times_swap_tmp = daycount.yearfrac(ref_date, swap_dates_tmp, basis_s)
            times_swap_new.append(times_swap_tmp)
            
            times_swap_x_tmp = daycount.yearfrac(fix_swap_dates, swap_dates_tmp, basis_s);
            times_swap_x_new.append(times_swap_x_tmp)
    

        '-----------------------------------------------------------------------------------------'
        '------generazione dei tempi a partire dalle date ----------------------------------------'
        '-----------------------------------------------------------------------------------------'

        # ------------------ da rivedere ------------------------------------------


        tenors           = times_swap_new[0], times_swap_new[1:] - times_swap_new[0:len(times_swap_new)-2]
        index_start      = find_indx(swap_times[0], times_swap_new)
    
        times_swap_tmp1s  = times_swap_new[:index_start-1]   #% tempi fino al primo tasso swap
        times_swap_tmp2s  = times_swap_new[index_start:]  # % tempi oltre il primo tasso swap
        dates_swap_tmp1s  = swap_dates_new[:index_start-1]   #% tempi fino al primo tasso swap
        nsw1s             = len(times_swap_tmp1s)
        nsw2s             = len(times_swap_tmp2s)
    
        df_tmp1s = np.zeros(1,nsw1s) #; % fattori di sconto fino al primo tasso swap
    
    
        
        #%-------------------------------------------------------------------------------------
        #%----------------------------- metodo lsr --------------------------------------------
        #%-------------------------------------------------------------------------------------
    
        #%----------------------- Interoplazione lineare dei tassi swap alle date di interesse ----------------------
    
        #%        swap_val_new  = interp1(swap_times, swap_val, times_swap_new(index_start:end), 'linear','extrap'); % interpolazione
        swap_val_new  = np.interp(swap_dates_new[index_start:], swap_dates, swap_val)#; % interpolazione
    
        #%-------------------------------------------------------------------------------------
        #%----------------------- Generazione fattori di sconto alle date dei pagamenti "Swap"
        #%--------------------------------------------------------------------------------------
    

        for i in range(0, nsw1s):
            if (merge_dates[len(merge_dates)-1] < dates_swap_tmp1s[i]):  #%--- NO INTERPOLAZIONE ---------------
    
                #%------------ gap colmato mantenedo costante il tasso fwd -----------
    
                if all(merge_df == 1):
                    fwd = swap_val[1]
                    df_tmp1s[i] = np.exp(-fwd*(times_swap_tmp1s[i]- merge_times[len(merge_times)-1]))
                else:
                    fwd         = -np.log(merge_df[len(merge_df)-1]/merge_df[len(merge_df)-2])/(merge_times[len(merge_times)] - merge_times[len(merge_times)-2])
                    df_tmp1s[i] = np.exp(-fwd*(times_swap_tmp1s[i]- merge_times[len(merge_times)-1]))
    
            else: #%--- INTERPOLAZIONE --------------------------
    
                index_0 = (dates_swap_tmp1s[i] in merge_dates)
    
                if (index_0 == None): # se il fattore di sconto non e' stato calcolato
    
                    #----------intepolazione esponenziale sui fattori di sconto ------------
    
                    index_1     = find_indx(dates_swap_tmp1s[i], merge_dates)
                    A           =  merge_df[index_1]
                    B           =  1./(merge_times[index_1]-merge_times[index_1 + 1 ])*np.log(merge_df[index_1 + 1]/merge_df[index_1])
                    df_tmp1s[i] = A*np.exp(-B*(times_swap_tmp1s[i]- merge_times[index_1]))
                
                else:
    
                    #%---------- viene selezionato il fattore di sconto gia' calcolato -------
    
                    df_tmp1s[i]= merge_df[index_0]
    
        z_all        = [df_tmp1s, np.zeros(1, nsw2s)]  #% tutti i fattori di sconto
        z_swap2s     = np.zeros(1, nsw2s) #% fattori di sconto a partire dal primo tasso swap
    
        #%-------------------------------------------------------------------------------------
        #%--------------- Generazione dei fattori di sconto alle date dei tassi swap ----------
        #%-------------------------------------------------------------------------------------
    
    
        for i in range(0,nsw2s):
            sum_z            = np.sum(z_all*tenors);
            z_swap2s[i]      = (df_fix_swap - swap_val_new[i]*sum_z)/(1.0 + swap_val_new[i]*tenors[nsw1s + i])
            z_all[nsw1s + i] = z_swap2s[i]
    
        #%--------------------------------  set LSR method OUTPUT --------------------------------
    
        if (flag_method_swap == 1): #% considero il metodo: Linear Swap Rate:'LSR'
            swap_df_tmp    = [1, z_all]
            swap_times_tmp = [0, times_swap_new]
            swap_times_tmp = np.round(swap_times_tmp/(tenor_swap/12.0))*(tenor_swap/12)
    
        # inserire la parte CFR-----------
        

        else:
            
         flag_swap  = [];
         swap_times = [];
         swap_dates = [];
         swap_discount_factor = [];
        
         
        #%------------------------------------------------------------------------------------
        #%--------------------------------  CFR METHOD ----------------------------------------%
        #%------------------------------------------------------------------------------------
        if (flag_method_swap == 0): # Considero il metodo "Constant fwd rate" 'CFR'
    
            n_irr = len(index_irregular) # %n. date irregolari
            index_irregular2 = np.round(swap_times_x[index_irregular(1,1)]/(tenor_swap/12)) #; % indice prima data irregolare
            index_irregular3 = np.round(swap_times_x[index_irregular(len(index_irregular)-1)]/(tenor_swap/12)) #; % indice ultima data irregolare
    
            z_out           = [z_all[:index_irregular2], np.zeros(1, index_irregular3 - index_irregular2)]
            times_swap_new  = times_swap_new + fix_swap_times
            times_tmp       = np.round(times_swap_x_new/(tenor_swap/12))*(tenor_swap/12)
            #  times_tmp       = round(times_swap_new./(tenor_swap/12)).*(tenor_swap/12);
    
    
            for i in range(0, n_irr-1):     #% ciclo su tutte le date irregolari
    
                n_i = index_irregular[i]     #% indice dell'n-ma data "irregolare"
                m_i = index_irregular[i+1]   #% indice dell'n-ma data "irregolare"
    
                swap_times_n  = swap_times_x[n_i]  #% tempo dell'n-ma data "irregolare"
                swap_times_m  = swap_times_x[m_i]  #% tempo dell'n-ma data "irregolare"
    
                id_swt_n = np.round(swap_times_n)   #% n. anni dell'n-ma data "irregolare"
                id_swt_m = np.round(swap_times_m)   #% n. anni dell'n-ma data "irregolare"
    
                nz  = find_indx(times_tmp, id_swt_n) #% indice del vettore times_out corrispondente a n anni
                mz  = find_indx(times_tmp, id_swt_m) #% indice del vettore times_out corrispondente a n anni
    
                #!!!index_stop = find_index_p(z_out~=0, 1, 'last'); #% selezione dell'(n-1)-mo indice
                fid_cond = z_out != 0
                index_stop = find_indx_p(False, fid_cond); #% selezione dell'(n-1)-mo indice
    
                z_n = z_out*tenors
                z_n = z_n[0:index_stop-1]
                sum_z   = np.sum(z_n)          #% sommatoria di tutti gli z fino all'(n-1)-mo
    
                f1 = 1.0/z_out[nz]*(df_fix_swap - swap_val[m_i]*sum_z) #% Def. f1 per il calcolo di F
    
                #%f1 rappresenta il termine noto della equazione (2.10) pag.38 del
                #%documento modelli_tassi_v02.pdf
    
                
                #f2 = @(x)fwd_fun(x, nz, mz, f1, swap_val[m_i], times_swap_new);
                
                mthod_o ='Newton-CG' #-->ok
    
                x0 = 0.02
            
                #res = optimize.minimize(fwd_fun, x0, method = mthod_o,  args=(nz, mz, f1, swap_val[m_i], times_swap_new), bounds = x_bnd, constraints=cons)
                res = optimize.minimize(fwd_fun, x0, method = mthod_o,  args=(nz, mz, f1, swap_val[m_i], times_swap_new))
                F = res.x
    
                #% Def. f2 per il calcolo di F
                #% f2 rappresenta la funzione di cui vogliamo trovere lo "zero", data
                #%dalla differenza tra primo e secondo membro della equazione (2.10) pag.38 del
                #%documento modelli_tassi_v02.pdf
    
                #F = fzero(f2,[0 1000]) #;               % Cerco gli zeri di f2
    
                #% calcolo i fattori di sconto alle date intermedie non quotate (di pagamento del tasso swap)
    
                for k in range(nz+1,mz):
                    z_out[k] = z_out[nz]*np.exp(-F*(times_swap_new[k] - times_swap_new[nz])) #; % calcolo il fattore di sconto z_out
                #end
            #end
    
            #%------------------ set CFR Output --------------------------------
    
            swap_df_tmp   = z_out
            swap_times_tmp = times_tmp
        #end
        #%---------------------------------------------------------------------
        #%------------------ set swap output --------------------------------
        #%---------------------------------------------------------------------
    
        swap_times_adj     =  np.round(swap_times_x/(tenor_swap/12))*(tenor_swap/12)
        index_select       =  np.zeros(1, nsw)
    
        for i in range(0, nsw):
            index_select[i]= find_indx(swap_times_adj(i), swap_times_tmp)
    
        swap_discount_factor = swap_df_tmp[index_select]
    
        #%---------------------------------------------------------------------------
        #%------------------- merge all data ---------------------------------------
        #%---------------------------------------------------------------------------
     
    else:
        flag_swap  = []
        swap_times = []
        swap_dates = []
        swap_discount_factor = []
    
    """     
    merge_flag_tmp  = [merge_flag,  flag_swap]
    merge_time_tmp  = [merge_times, swap_times]
    merge_dates_tmp = [merge_dates, swap_dates]
    merge_df_tmp    = [merge_df,   swap_discount_factor]
    merge_data_tmp  = [merge_time_tmp, merge_dates_tmp, merge_df_tmp, merge_flag_tmp]
    """
    
    swap_rates = swap_discount_factor
    
    df_array = swap_discount_factor
    r_array  = swap_rates
    t_array  = swap_times
    d_array  = swap_dates

    data_merge = merge_data(data_merge, d_array, t_array, r_array, df_array)


    
    #sort_data   = sortrows(merge_data_tmp, 1)
    """
    sort_data   = (merge_data_tmp)

    #merge_times = sort_data(:,1);
    merge_dates = sort_data[2]
    merge_df    = sort_data[3]
    merge_flag  = sort_data[4]
    n_dmerge    = len(merge_df);
    merge_times = (merge_dates-merge_dates[0])/360; 
    """

    #%--------------------------------------------------------------------------
    #%-------------- CALCOLO Tassi OUT  ----------------------------------------
    #%--------------------------------------------------------------------------


    merge_rates = compute_rates(merge_df, merge_times, regime_output)
    
    
    if all(merge_df == 1):
        merge_rates[0] = 0;
    else:
        merge_rates[0] = merge_rates[1]
    
    #%-------------------------------------------------------------------------
    #%-------------- SET OUTPUT  ROUTINES ----------------------------------------------
    #%-------------------------------------------------------------------------
    
    
    discount_factors_out = merge_df;
    dates_out            = merge_dates;
    times_out            = merge_times;
    rendimenti_out       = merge_rates;
    #flag_data            = merge_flag;
    
    #%-------------------------------------------------------------------------
    #%-------------- GRAFICI DI CONTROLLO  ------------------------------------
    #%-------------------------------------------------------------------------
    
    """
    if (flag_1s == 1):
        index_seg1  = find_indx(flag_data == 1);
    elif (flag_1s == 2):
        index_seg1  = find_indx(flag_data == 2);
    
    index_futures  = find_indx(flag_data == 3);
    index_swap     = find_indx(flag_data == 4);
    
    seg1_df_g      = merge_df(index_seg1);
    seg1_times_g   = merge_times(index_seg1);
    seg1_rates_g   = merge_rates(index_seg1);
    
    futures_df_g    = merge_df(index_futures);
    futures_times_g = merge_times(index_futures);
    futures_rates_g = merge_rates(index_futures);
    
    swap_df_g    = merge_df(index_swap);
    swap_times_g = merge_times(index_swap);
    swap_rates_g = merge_rates(index_swap);
    
    """


    '-----------------------------------------------------------------------------------------'
    '--------------- make graph ------------------------------------------------------------'
    '-----------------------------------------------------------------------------------------'

    
    if (flag_make_graph == 1):

    
        g1 = make_graph_for_chk(seg1_times_g, seg1_rates_g, futures_times_g, futures_rates_g, swap_times_g, swap_rates_g, seg1_val, futures_val, swap_val)
        g2 = make_graph_for_chk(seg1_times_g, seg1_df_g, futures_times_g, futures_df_g, swap_times_g, swap_df_g, seg1_val, futures_val, swap_val)

    else:
        
        flag_save_graph = 0
    #%------------------- Salvataggio grafici ----------------------------------
    
    
    if (flag_save_graph == 1):
    
        data_ref_str = '31-03-2016'
        
        g1.savefig('fig/plot_tassi_%s.png' %(data_ref_str))
        g2.savefig('fig/plot_discount_%s.png' %(data_ref_str))


    '--------------- save graph ------------------------------------------------------------'
    
 
    date_out   = []
    tempi_out  = []
    valori_out = []

    nodi_out = convertTimeToNode(tempi_out)
    #nodi_out   = []

    data_elab_out = {}
    
    data_elab_out['Date scadenza'] = date_out        
    data_elab_out['Valori']        = valori_out
    data_elab_out['Tempi']         = tempi_out
    data_elab_out['Nodi']          = nodi_out
        
 
    
    
    return data_elab_out




if __name__ == "__main__":
 
    

    #FQ(223)
    
    data_opt = {}
    data_raw = {}
    
    data_opt = set_data_opt()    
    
    
    #------------------------------------------------------------
    

    inputCurveFile = r"input_test\raw_swap_curves.txt"
    data_raw = load_data_fromFile(inputCurveFile)

    ref_field = 'UsaNodo'
    val_not_allowed = 'N'

    data_raw_p = purge_data(data_raw, ref_field, val_not_allowed)
    
    dict1s, dict_f, dict_s = select_segments(data_raw_p)
    
    
    print 'dict1s: ', dict1s
    
    FQ(2233)
    #print data_raw_p.keys()
    dates_0 = data_raw_p['MatDate']
    
    
    
    print data_raw_p['UsaNodo']
    
    FQ(2229)
    time_n = daycount.yearfractions(dates_0, 'ACT/360')
    
    calendar_code = 'de'
    
    mkt_ref = holidays.get_calendar(calendar_code)
    
    day_conv = 'follow'
    
    dates_0_0 = (dates_0[0])
    
    print 'dates_0_0: ', dates_0_0    
    print 'type(dates_0_0): ', type(dates_0_0)
    #FQ(888)
    
    dates_0_1 = busdayrule.rolldates(dates_0, mkt_ref, day_conv)
    
    print 'dates_0: ', dates_0 
    print 'dates_0_1: ', dates_0_1 

    #------------------------------------------------------------------
    #FQ(111)   
    data_elab_out = boot3s_elab(data_opt, data_raw_p)
    
    FQ(223)

    
 
    dict_opt ={}
    dict_data = {}
    
    #time_ref = [1.0/52, 2.0/52, 5.0/52, 1.0/12, 6.0/12, 9/12, 1.0, 14.0/12.0, 24.0, 36]
    #time_ref = [3.0/52.0, 1.0/12.0, 1, 2, 3, 4, 1.2]
    #time_ref = [3.0/52.0, 1.0/12.0]

    time_ref = [1.0/24, 1.0/12, 1/2, 1, 10]
    
    node_out = convertTimeToNode(time_ref)
    
    for i in range(len(time_ref)):
        
        timeTmp = time_ref[i]
        nodeTmp = node_out[i]
        
        print 'Time tmp: ', timeTmp 
        print 'Node tmp: ', nodeTmp
        print '----------------------------' 
    
    #dict_curve_elab = boot3s_elab(dict_opt, dict_data)
    