import sys
#import datetime as dtime
import math
import numpy as np

#from sc_elab import numpy as np
from sc_elab.core.mdates import holidays
from sc_elab.core.mdates import daycount
from sc_elab.core.mdates import busdayrule

from scipy import optimize
from scipy.optimize import minimize
from scipy.optimize import fmin
from datetime import datetime as dtime
from dateutil.relativedelta import relativedelta

import matplotlib.pyplot as plt

import datetime
import calendar



def FQ(label):
    print ('------------- FIN QUI TUTTO OK  %s ----------' %(label))
    sys.exit()




def find_indx_last(x_target, x_list_ref):

    ln = len(x_list_ref) - 1

    for i in range(0, len(x_list_ref)):

        if (x_target < x_list_ref[0]): 
            i_ref = None
            break

        elif (x_target > x_list_ref[ln]) : 
            i_ref = None
            break

        elif (x_list_ref[i] == x_target) and (i == 0): 
            i_ref = i
            break

        elif (x_list_ref[i] > x_target):
            i_ref = i - 1
            break
            
        else:

            continue

    return i_ref



def zc_rate_by_SVE(dict_params, t_mat):

    t_mat   = np.asarray(t_mat)

    minTime = 0.0001
    
    for i in range(0, len(t_mat)):
        t_matTmp = t_mat[i]
        if (t_matTmp < minTime): t_matTmp = minTime
        t_mat[i] = float(t_matTmp)
        
    
    tau1   = dict_params['const1'][0]
    tau2   = dict_params['const2'][0]
    beta0  = dict_params['beta0'][0]
    beta1  = dict_params['beta1'][0]
    beta2  = dict_params['beta2'][0]
    beta3  = dict_params['beta3'][0]

    tmp1 = t_mat/tau1
    tmp2 = t_mat/tau2
    g1 = np.exp(-tmp1)
    g2 = np.exp(-tmp2)

    G1 = 1.0 - g1

    ns_model  =  beta0 + (beta1 + beta2)*(G1/tmp1) - beta2*g1
    zc_rate   =  ns_model - beta3*(g2 - (1.0 - g2)/tmp2)

    return zc_rate



def zc_rate_by_NS(dict_params, t_mat):

    t_mat   = np.asarray(t_mat)

    minTime = 0.0001
    
    for i in range(0, len(t_mat)):
        t_matTmp = t_mat[i]
        if (t_matTmp < minTime): t_matTmp = minTime
        t_mat[i] = float(t_matTmp)

    tau1   = dict_params['const1'][0]
    beta0  = dict_params['beta0'][0]
    beta1  = dict_params['beta1'][0]
    beta2  = dict_params['beta2'][0]

    tmp1 = t_mat/tau1
    g1 = np.exp(-tmp1)

    G1 = 1.0 - g1

    ns_model  =  beta0 + (beta1 + beta2)*(G1/tmp1) - beta2*g1
    zc_rate   =  ns_model

    return zc_rate


def zc_rate_by_CIR(dict_params, t_mat):

    t_mat = np.asarray(t_mat)
 
    minTime = 0.0001
    
    for i in range(0, len(t_mat)):
        t_matTmp = t_mat[i]
        if (t_matTmp < minTime): t_matTmp = minTime
        t_mat[i] = t_matTmp

    
    r0     = dict_params['r0'][0]
    kappa  = dict_params['kappa'][0]
    theta  = dict_params['theta'][0]
    sigma  = dict_params['sigma'][0]
    
    T = t_mat

    #%Lettura dei parametri
    
    h = (kappa*kappa + 2.0*sigma*sigma)**(0.5)

    g0 = 2*kappa*theta/(sigma*sigma)
    g1 = np.exp(T*h) - 1.0
    g2 = np.exp(T*(h + kappa)/2.0)

    A0 = (2*h*g2/(2.0*h + (kappa + h)*g1))
    B0 = (2.0*g1/(2.0*h + (kappa + h)*g1))
    
    zc_rate = -(g0*np.log(A0) - B0*r0)/T

    return zc_rate



def zc_rate_by_LIN(tempi,parameters, T_target):

    
    zc_rate = []
    
    type_ndarray = isinstance(T_target, np.ndarray)
    type_list = isinstance(T_target, list)
    
    if( type_list and type_ndarray) == 'FALSE':
            zc_rate = zc_rate_by_LIN_s(tempi, parameters, T_target)
            zc_rate = [zc_rate]
    
    else:

        for i in range(0, len(T_target)):    
            t_tmp = T_target[i]
            zc_outTmp = zc_rate_by_LIN_s(tempi,  parameters, t_tmp)            
            zc_rate.append(zc_outTmp)
    
    return zc_rate


def zc_rate_by_LIN_s(tempi, parameters, T_target):

    t = tempi
    n_knots  = len(t)

    # Se la maturity T e' negativa, pongo Z=-1
    if (T_target < 0):
        zc_rate = -1
        return  zc_rate  

    elif (T_target == 0): #se il tempo T e' nullo Z=1
        
        zc_rate =  parameters['b'][0]
        return zc_rate[0]

    elif (T_target > 0) and (T_target < t[n_knots-1]):
        index = find_indx_last(T_target, tempi)        
        zc_rate = parameters['a'][index]/T_target + parameters['b'][index]
        zc_rate = zc_rate[0] 

    elif T_target >= t[n_knots-1]:
        index = n_knots - 2
        
        zc_rate = parameters['a'][index]/T_target + parameters['b'][index]
        zc_rate = zc_rate[0] 
    
    return zc_rate

def zc_rate_by_AVD(tempi, zc_values, prms, T_target):

    zc_rate = []
    df_values = []
    
    for i in range(0, len(tempi)):
        
        t_tmp = tempi[i]
        df_tmp = np.exp(-zc_values[i]*t_tmp)
        df_values.append(df_tmp)
    
    type_ndarray = isinstance(T_target, np.ndarray)
    type_list = isinstance(T_target, list)
    
    if( type_list and type_ndarray) == 'FALSE':

            zc_rate = zc_rate_by_AVD_s(tempi, df_values, prms, T_target)
            zc_rate = [zc_rate]
    else:

        for i in range(0, len(T_target)):    

            t_tmp = T_target[i]          
            zc_outTmp = zc_rate_by_AVD_s(tempi, df_values, prms, t_tmp) 
            
            zc_rate.append(zc_outTmp)
        
    return zc_rate

def zc_rate_by_AVD_s(tempi, df_v, prms, T_target):
    
    if (T_target < 0.001): T_target = 0.001
    
    df_value = df_by_AVD_s(tempi, df_v, prms, T_target)
    
    
    zc_rate_out = -np.log(df_value)/T_target

    return zc_rate_out 

def df_by_AVD_s(tempi, df_v, prms, T_target):

    #Dichiarazione variabili
    t = tempi
    T = T_target
    
    n_knots = len(t)
    
    #%Se la maturity T e' negativa, pongo Z=-1

    if (T < 0):
        zc_rate = -1
        return zc_rate

        #Se il tempo T e' nullo Z=1
    elif (T==0):
        df_value =1.0
        return df_value
    
        #Se e' compreso tra 0 e il tempo del primo nodo
    elif ( T > 0) and  ( T < t[0]):
        i = 0
        Tp = 0
        Zp = 1
        
    #Se e' compreso tra il tempo del primo nodo ed il tempo dell%ultimo nodo
    #%individuo la fascia temporale di appartenenza
    elif (T >= t[0]) and (T < t[n_knots-1]):
    
        i = find_indx_last(T,t)
        Tp = t[i]
        Zp = df_v[i]
    
    elif (T >= t[n_knots-1]):
        i = n_knots-2
        Tp = t[n_knots-1]
        Zp = df_v[n_knots-1]
    
    #%--------------------------------------------------------------------------
    #%Calcolo il prezzo ZC usando i parametri
    #%secondo l'algoritmo di Adams & Van Denventer (1994)
    #%--------------------------------------------------------------------------
   
    #zc_rate = parameters['a'][index]/T_target + parameters['b'][index] 
    #dummy = prms['a'][i+1]*(T) + (prms['b'][i+1]/2)*(T^2)+(prms['c'][i+1]/3)*(T ^3) + (prms['d'][i+1]/4)*(T^4) + (prms['e'][i+1]/5)*(T^5)
    #zc_rate = dummy/T
    zc_rate_t = prms['a'][i]*(T-Tp) + (prms['b'][i]/2)*(T**2-Tp**2)+(prms['c'][i]/3)*(T**3-Tp**3) + (prms['d'][i]/4)*(T**4-Tp**4) + (prms['e'][i]/5)*(T**5-Tp**5)
    
    df_out = Zp*np.exp(-zc_rate_t)
    
    return df_out

def df_linear_s(tempi, discount_factors, parameters, T_target):

    t = tempi;
    n_knots  = len(t);
    
    # Se la maturity T e' negativa, pongo Z=-1
    if (T_target < 0):
        z_out = -1
        return    

    elif (T_target == 0): #se il tempo T e' nullo Z=1
        z_out = 1
        return

    elif (T_target > 0) and (T_target < t[n_knots-1]):
        index = find_indx_last(T_target, tempi)
        dummy = parameters[index,0] + parameters[index,1]*T_target 
        z_out     = np.exp(-dummy)

    elif T_target >= t[n_knots-1]:
        index = n_knots
        dummy = parameters[index,0] + parameters[index,1]*T_target 
        z_out    = np.exp(-dummy)
    
    return z_out





def estimate_linear_params(t_times, zc_rates):


    zc_rates = np.asarray(zc_rates)
    t_times = np.asarray(t_times)

    discount_factors = np.exp(-t_times*zc_rates)

    t = t_times 
    # Calcolo numero nodi
    n_knot = len(t_times) 
    # Numero funzioni da stimare
    n_functions = n_knot - 1
    # numero parameteri per funzione
    n_parameters = 2
        
    # Dimensionamento vettori e matrice
    m = n_functions*n_parameters

    #lin_parameters = np.zeros((n_knot-1,n_parameters))    
    lin_parameters = {}
    lin_parameters['a'] = {}
    lin_parameters['b'] = {}

    system_coeff   = np.zeros((m,m))
    known_term     = np.zeros((m,1))
    
    # Riempimento matrice dei coefficienti del sistema e del vettore termine noto
    # Le condizioni usate sono a pag.42-43 del documento 'modelli_tassi'
    for k in range(0, n_functions-1):
        for j in range(0, n_parameters):
            # 1ma condizione: f=f
            system_coeff[k][k*n_parameters+j] = 1.0*t[k+1]**(j)
            system_coeff[k][k*n_parameters+n_parameters+j] = -1.0*t[k+1]**(j)
    
    for k in range(0, n_functions-1):
        for j in range(0, n_parameters):
            # 2 condizione che sul nodo il fattore di sconto si a pari
            # a quello quotato sul mercato Z_lin = Z__mkt
            system_coeff[(n_functions-1)+k][(k)*n_parameters+j] = 1*t[k+1]**(j)
            known_term[(n_functions-1)+k] = -1.0*np.log(discount_factors[k+1])
    
    #  condizione: sul nodo zero (0)
    system_coeff[n_parameters*(n_functions-1),0] = 1
    known_term[n_parameters*(n_functions-1)] = 0
    
    #  condizioni: sull'ultimo nodo (n_knot)
    for j in range(0,n_parameters):
        system_coeff[n_parameters*(n_functions-1)+1][(n_functions-1)*n_parameters+j] = 1.0*t[n_knot-1]**(j)

    known_term[n_parameters*n_functions-1] = -1.0*np.log(discount_factors[n_knot-1])
        
    # Risoluzione del sistema
    
    solution = np.linalg.solve(system_coeff, known_term)

    # Disposizione della soluzione sottoforma di matrice
    for i in range(0, n_knot-1):
        lin_parameters['a'][i] = solution[(i)*n_parameters]
        lin_parameters['b'][i] = solution[(i)*n_parameters + 1]

    return lin_parameters 



def estimate_avd_params(curve_times, zc_rate):

    zc_rate = np.asarray(zc_rate)
    curve_times = np.asarray(curve_times)

    discount_factors = np.exp(-curve_times*zc_rate)

  
    t = curve_times
    
    # Calcolo numero nodi
    n_knot = len(curve_times) 
    # Numero funzioni da stimare
    n_functions = n_knot - 1
    # numero parameteri per funzione
    n_parameters = 5
        
    # Dimensionamento vettori e matrice
    m = n_functions*n_parameters
    
    system_coeff     =  np.zeros((m,m))
    known_term       =  np.zeros((m,1))
    lin_parameters   =  np.zeros((n_knot-1,n_parameters))
    
    
    
    # Calcolo tassi ZC continui sui nodi
    r_cont_t0 = -1/t[1]*np.log(discount_factors[1])
    
    
    # Riempimento matrice dei coefficienti del sistema e del vettore termine noto
    # Le condizioni usate sono a pag.42-43 del documento 'modelli_tassi'
    
    
    avd_parameters = {}
    avd_parameters['a'] = {}
    avd_parameters['b'] = {}
    avd_parameters['c'] = {}
    avd_parameters['d'] = {}
    avd_parameters['e'] = {}


    for k in range(0, n_functions-1):
        for j in range(0, n_parameters):
            # 1 condizione: f=f 
            system_coeff[k][k*n_parameters+j] = 1.0*t[k+1]**(j)
            system_coeff[k][k*n_parameters+n_parameters+j] = -1.0*t[k+1]**(j)
            
            # 2 condizione f'=f' 
            
            system_coeff[(n_functions-1)+k][(k)*n_parameters+j] =(j)*t[k+1]**(j-1)
            system_coeff[(n_functions-1)+k][(k)*n_parameters+n_parameters+j] =-(j)*t[k+1]**(j-1)

        
        # 3 condizione f''=f''
        
        system_coeff[2*(n_functions-1)+k][(k)*n_parameters+0] = 0
        system_coeff[2*(n_functions-1)+k][(k)*n_parameters+1] = 0
        system_coeff[2*(n_functions-1)+k][(k)*n_parameters+2] = 2*t[k+1]**0
        system_coeff[2*(n_functions-1)+k][(k)*n_parameters+3] = 6*t[k+1]**1
        system_coeff[2*(n_functions-1)+k][(k)*n_parameters+4] = 12*t[k+1]**2
        system_coeff[2*(n_functions-1)+k][(k)*n_parameters+n_parameters+0] = 0
        system_coeff[2*(n_functions-1)+k][(k)*n_parameters+n_parameters+1] = 0
        system_coeff[2*(n_functions-1)+k][(k)*n_parameters+n_parameters+2] = -2*t[k+1]**0
        system_coeff[2*(n_functions-1)+k][(k)*n_parameters+n_parameters+3] = -6*t[k+1]**1
        system_coeff[2*(n_functions-1)+k][(k)*n_parameters+n_parameters+4] = -12*t[k+1]**2
    
        # 4 condizione f'''=f''' 
        system_coeff[3*(n_functions-1)+k][(k)*n_parameters+0] = 0
        system_coeff[3*(n_functions-1)+k][(k)*n_parameters+1] = 0
        system_coeff[3*(n_functions-1)+k][(k)*n_parameters+2] = 0
        system_coeff[3*(n_functions-1)+k][(k)*n_parameters+3] = 6*t[k+1]**0
        system_coeff[3*(n_functions-1)+k][(k)*n_parameters+4] = 24*t[k+1]**1
        system_coeff[3*(n_functions-1)+k][(k)*n_parameters+n_parameters+0] = 0
        system_coeff[3*(n_functions-1)+k][(k)*n_parameters+n_parameters+1] = 0
        system_coeff[3*(n_functions-1)+k][(k)*n_parameters+n_parameters+2] = 0
        system_coeff[3*(n_functions-1)+k][(k)*n_parameters+n_parameters+3] =-6*t[k+1]**0
        system_coeff[3*(n_functions-1)+k][(k)*n_parameters+n_parameters+4] = -24*t[k+1]**1
    
        # 5 condizione  -ln(discount_factors(k)/discount_factors(k-1)), uguaglianza con fattori sconto

        for k in range(0, n_functions):
            for j in range(0, n_parameters):
                    
                t_x = (1.0/(j+1))*(t[k+1]**(j+1)-t[k]**(j+1))
                t_y = -1.0*np.log(discount_factors[k+1]/discount_factors[k])
                indx_r = 4*(n_functions-1)+k
                indx_c = k*n_parameters+j
                
                """
                print 'indx_r: ', indx_r
                print 'indx_c: ', indx_c
                print 't_x: ', t_x
                """
                
                system_coeff[indx_r][indx_c] = t_x
                known_term[indx_r] = t_y
                
        # 6 e 7 condizioni: sul nodo zero (0)
        system_coeff[n_parameters*n_functions-4][0] = 1
        
        known_term[n_parameters*n_functions-4] = r_cont_t0;
        
        system_coeff[n_parameters*n_functions-3][0] = 0
        system_coeff[n_parameters*n_functions-3][1] = 0
        system_coeff[n_parameters*n_functions-3][2] = 2*t[0]**0
        system_coeff[n_parameters*n_functions-3][3] = 6*t[0]**1
        system_coeff[n_parameters*n_functions-3][4] = 12*t[0]**2
        
        # 8 e 9 condizioni: sull'ultimo nodo (n_knot)
        for j in range(1, n_parameters):
            
            indx_r = n_parameters*n_functions-2
            indx_c = (n_functions-1)*n_parameters+j
            t_x = (j)*t[n_knot-1]**(j)

            """
            print 'indx_r: ', indx_r
            print 'indx_c: ', indx_c
            print 't_x: ', t_x
            """
            system_coeff[indx_r][indx_c] = t_x
         
        system_coeff[n_parameters*n_functions-1][(n_functions-1)*n_parameters] = 0
        system_coeff[n_parameters*n_functions-1][(n_functions-1)*n_parameters+1] = 0
        system_coeff[n_parameters*n_functions-1][(n_functions-1)*n_parameters+2] = 2*t[n_knot-1]**0
        system_coeff[n_parameters*n_functions-1][(n_functions-1)*n_parameters+3] = 6*t[n_knot-1]**1
        system_coeff[n_parameters*n_functions-1][(n_functions-1)*n_parameters+4] = 12*t[n_knot-1]**2
         
        
        # Risoluzione del sistema
    solution = np.linalg.solve(system_coeff, known_term)
    

    # Disposizione della soluzione sottoforma di matrice
    
    for i in range(0, n_knot-1):
        
        avd_parameters['a'][i] = solution[(i)*n_parameters]
        avd_parameters['b'][i] = solution[(i)*n_parameters + 1]
        avd_parameters['c'][i] = solution[(i)*n_parameters + 2]
        avd_parameters['d'][i] = solution[(i)*n_parameters + 3]
        avd_parameters['e'][i] = solution[(i)*n_parameters + 4]

    return avd_parameters 
    



def loss_fun_for_fitting(par, t_mkt, zc_mkt, model_type):
        
        
        
        if (model_type == '2'): # modello SVE
            dict_params = {}

            dict_params['const1'] = [par[0]]
            dict_params['const2'] = [par[1]]
            dict_params['beta0']  = [par[2]]
            dict_params['beta1']  = [par[3]]
            dict_params['beta2']  = [par[4]]
            dict_params['beta3']  = [par[5]]
                        
            zc_model = zc_rate_by_SVE(dict_params, t_mkt)
            
        elif (model_type == '3'): # modello CIR

            dict_params = {}
            
            dict_params['r0']     = [par[0]]
            dict_params['kappa']  = [par[1]]
            dict_params['theta']  = [par[2]]
            dict_params['sigma']  = [par[3]]
                        
            zc_model = zc_rate_by_CIR(dict_params, t_mkt)

        elif (model_type == '4'): # modello NS

            dict_params = {}
            
            dict_params['const1'] = [par[0]]
            dict_params['beta0']  = [par[1]]
            dict_params['beta1']  = [par[2]]
            dict_params['beta2']  = [par[3]]
                        
            zc_model = zc_rate_by_NS(dict_params, t_mkt)
        
        else:
            
            dict_params = {}
            
            dict_params['const1'] = [par[0]]
            dict_params['const2'] = [par[1]]
            dict_params['beta0']  = [par[2]]
            dict_params['beta1']  = [par[3]]
            dict_params['beta2']  = [par[4]]
            dict_params['beta3']  = [par[5]]
            
            zc_model = zc_rate_by_SVE(dict_params, t_mkt)
            
        zc_mkt = np.array(zc_mkt)
        
        #print 'dict_params: ', dict_params
        

        zc_mkt = zc_mkt*10000.0
        zc_model = zc_model*10000.0

        #zc_mkt = zc_mkt
        #zc_model = zc_model
        #n_pt = len(zc_mkt)

        diff = abs((zc_mkt - zc_model))
        
        diff2 = diff*diff

        
        val_sum = np.sum(diff2)
        val_sum = val_sum
        
        return val_sum


def makeRatesFromModel(mkt_times, mkt_values, dict_model_par, target_times, model_type):



    if (model_type == '0'):
        mdl_values = zc_rate_by_LIN(mkt_times, dict_model_par, target_times)

    elif (model_type == '1'):
        mdl_values = zc_rate_by_AVD(mkt_times, mkt_values, dict_model_par, target_times)
        
    elif (model_type == '2'):
        mdl_values = zc_rate_by_SVE(dict_model_par, target_times)

    elif (model_type == '3'):
        mdl_values = zc_rate_by_CIR(dict_model_par, target_times)

    elif (model_type == '4'):
        mdl_values = zc_rate_by_NS(dict_model_par, target_times)

    else:
        mdl_values = zc_rate_by_CIR(dict_model_par, target_times)



    return mdl_values 



def chk_graph(mkt_times, mkt_values, model_type, dict_model_par):


    model_dict = {}
    model_dict['0'] = 'LIN'
    model_dict['1'] = 'AVD'
    model_dict['2'] = 'SVE'
    model_dict['3'] = 'CIR'
    model_dict['4'] = 'NS'
    
    modelRef = model_dict[model_type]

    mkt_times = np.asarray(mkt_times)
    target_times = mkt_times

    mdl_values = makeRatesFromModel(mkt_times, mkt_values, dict_model_par, target_times, model_type)
        
    
    #print 'mdl_values: ', mdl_values
    
    plt.plot(mkt_times, mkt_values, 'o', label='market')
    plt.plot(target_times, mdl_values, '-', label='model')

    plt.xlabel('Tempi [anni]')
    plt.ylabel('Livello tassi zc')
    plt.title('Fitting via %s model' %(modelRef))

    plt.legend()
    plt.show()




def  fitting(c_dates, c_values, opt_dict):
    
    #print 'c_dates: ', c_dates
    #print 'c_values: ', c_values
    #print 'opt_dict: ', opt_dict
    
    c_values = np.array(c_values)
    #c_values = c_values/100.0
    
    
    
    try:
        opt_dict['MakeGraph']
    except:
        opt_dict['MakeGraph'] = False

        
    #opt_dict['MakeGraph'] = True
    
    t_mkt = []
    for i in range(0, len(c_dates)):
        
        dateDiff = c_dates[i] - c_dates[0]
        timeTmp = float(dateDiff.days/365.2425)
        t_mkt.append(timeTmp)
    
    model_type = opt_dict['interp']
    
    
    #mthod_o ='SLSQP' #ok
    mthod_o ='TNC' #ok

    zc_mkt = c_values
    #zc_mkt = data_raw['ValoreZC']
    
    #t_mkt  = data_raw['MatTimes']
    t_mkt  = np.asarray(t_mkt)


    
    if (model_type == '0'): # caso lineare

        prms_lin = estimate_linear_params(t_mkt, zc_mkt)
        
    elif (model_type == '1'): # caso AVD
        
        prms_avd = estimate_avd_params(t_mkt, zc_mkt)
        
    else:

        x_bnd = []
        x0 = []
        
        if (model_type == '2'): # caso SVE
            
            ln_prms = len(opt_dict['bound_min_sve'])

            for i in range(0, ln_prms):
                
                x_bnd_minTmp = opt_dict['bound_min_sve'][i]
                x_bnd_maxTmp = opt_dict['bound_max_sve'][i]
                #x0Tmp        = opt_dict['x0'][i]
                
                x0Tmp        = float((x_bnd_maxTmp + x_bnd_minTmp)/2.0)
                
                x_bnd.append([x_bnd_minTmp, x_bnd_maxTmp])
                x0.append(x0Tmp)


        if (model_type == '3'): # caso CIR
        
            ln_prms = len(opt_dict['bound_min_cir'])
            
            for i in range(0, ln_prms):
                
                x_bnd_minTmp = opt_dict['bound_min_cir'][i]
                x_bnd_maxTmp = opt_dict['bound_max_cir'][i]
                #x0Tmp = opt_dict['x0'][i]
                x0Tmp        = float((x_bnd_maxTmp + x_bnd_minTmp)/2.0)
                
                x_bnd.append([x_bnd_minTmp, x_bnd_maxTmp])
                x0.append(x0Tmp)
        
        if (model_type == '4'): # caso NS
        
            ln_prms = len(opt_dict['bound_min_ns'])
            
            for i in range(0, ln_prms):
                
                x_bnd_minTmp = opt_dict['bound_min_ns'][i]
                x_bnd_maxTmp = opt_dict['bound_max_ns'][i]
                x_bnd.append([x_bnd_minTmp, x_bnd_maxTmp])

                x0Tmp        = float((x_bnd_maxTmp + x_bnd_minTmp)/2.0)

                #x0Tmp = opt_dict['x0'][i]
                
                x0.append(x0Tmp)


        res = optimize.minimize(loss_fun_for_fitting, x0, method = mthod_o,  args=(t_mkt, zc_mkt, model_type), bounds = x_bnd)

    dict_model_par = {}
    
    """
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
    """
    

    if (model_type == '0'): # Linear

        date_prms = []
        a_prms = []
        b_prms = []
        
        values_list = prms_lin['a'].keys()
        ln = len(values_list)
        
        #date_prms.append(c_dates[0])
        #a_prms.append(None)
        #b_prms.append(None)
        
        for i in range(0, ln):
        
            
            #dateTmp = data_raw['MatDate'][i]
            dateTmp =c_dates[i+1]

            a_prmsTmp = prms_lin['a'][i]
            b_prmsTmp = prms_lin['b'][i]

            a_prms.append(a_prmsTmp)
            b_prms.append(b_prmsTmp)
            date_prms.append(dateTmp)


        dict_model_par['a']      = a_prms
        dict_model_par['b']      = b_prms
        dict_model_par['Dates']  = date_prms

    elif (model_type == '1'): # AVD

        date_prms = []
        a_prms = []
        b_prms = []
        c_prms = []
        d_prms = []
        e_prms = []

        values_list = prms_avd['a'].keys()
        ln = len(values_list)
        
        for i in range(0, ln):
        
            #dateTmp = data_raw['MatDate'][i]
            dateTmp = c_dates[i+1]

            a_prmsTmp = prms_avd['a'][i]
            b_prmsTmp = prms_avd['b'][i]
            c_prmsTmp = prms_avd['c'][i]
            d_prmsTmp = prms_avd['d'][i]
            e_prmsTmp = prms_avd['e'][i]


            date_prms.append(dateTmp)

            a_prms.append(a_prmsTmp)
            b_prms.append(b_prmsTmp)
            c_prms.append(c_prmsTmp)
            d_prms.append(d_prmsTmp)
            e_prms.append(e_prmsTmp)

        dict_model_par['a']      = a_prms
        dict_model_par['b']      = b_prms
        dict_model_par['c']      = c_prms
        dict_model_par['d']      = d_prms
        dict_model_par['e']      = e_prms

        dict_model_par['Dates']  = date_prms

    elif (model_type == '2'): # Svensson

        #dataRef = data_raw['MatDate'][0]
        dataRef = c_dates[0]

        dict_model_par['Dates']  = [dataRef]
        dict_model_par['const1'] = [res.x[0]]
        dict_model_par['const2'] = [res.x[1]]
        dict_model_par['beta0']  = [res.x[2]]
        dict_model_par['beta1']  = [res.x[3]]
        dict_model_par['beta2']  = [res.x[4]]
        dict_model_par['beta3']  = [res.x[5]]
        
    elif (model_type == '3'): # CIR

        #dataRef = data_raw['MatDate'][0]
        dataRef = c_dates[0]

        dict_model_par['Dates'] = [dataRef]
        dict_model_par['r0']    = [res.x[0]]
        dict_model_par['kappa'] = [res.x[1]]
        dict_model_par['theta'] = [res.x[2]]
        dict_model_par['sigma'] = [res.x[3]]

    elif (model_type == '4'): # NS

        #dataRef = data_raw['MatDate'][0]
        dataRef = c_dates[0]

        dict_model_par['Dates'] = [dataRef]
        dict_model_par['const1']    = [res.x[0]]
        dict_model_par['beta0'] = [res.x[1]]
        dict_model_par['beta1'] = [res.x[2]]
        dict_model_par['beta2'] = [res.x[3]]



    else:

        #dataRef = data_raw['MatDate'][0]
        dataRef = c_dates[0]
        
        dict_model_par['r0']    = [res.x[0]]
        dict_model_par['Dates'] = [dataRef ]      
        dict_model_par['kappa'] = [res.x[1]]
        dict_model_par['theta'] = [res.x[2]]
        dict_model_par['sigma'] = [res.x[3]]
    
    #--------------- MAKE GRAPH -------------------
    make_graph = opt_dict['MakeGraph']
    
    if (make_graph == True):

        dataRef = c_dates[0]

        model_type = opt_dict['interp']
        
        mkt_times  = t_mkt
        mkt_values = c_values
        
        chk_graph(mkt_times, mkt_values, model_type, dict_model_par)


    return dict_model_par



def packModelPrms(tipo_modello, bench_dates, mdl_prms_list):
    
        prms_dict = {}
        par1List  = []
        par2List  = []
        
        bench_dates_list = []
        if (tipo_modello == 'LIN'):
        
            ln = len(bench_dates)
            for i in range(1, ln):
                
                dateTmp = bench_dates[i]
                par1Tmp    = mdl_prms_list[2*(i-1)]
                par2Tmp    = mdl_prms_list[2*(i-1) + 1]
                bench_dates_list.append(dateTmp)
                
                par1Tmp = np.array([par1Tmp])
                par2Tmp = np.array([par2Tmp])
                
                par1List.append(par1Tmp)
                par2List.append(par2Tmp)
                
            prms_dict['Dates'] = bench_dates_list
            prms_dict['a'] = par1List
            prms_dict['b'] = par2List

        elif (tipo_modello == 'SVE'):
            
            prms_dict['Dates']  = [bench_dates[0]]
            prms_dict['const1'] = [mdl_prms_list[0]]
            prms_dict['const2'] = [mdl_prms_list[1]]
            prms_dict['beta0']  = [mdl_prms_list[2]]
            prms_dict['beta1']  = [mdl_prms_list[3]]
            prms_dict['beta2']  = [mdl_prms_list[4]]
            prms_dict['beta3']  = [mdl_prms_list[5]]

        elif (tipo_modello == 'NS'):
            
            prms_dict['Dates']  = [bench_dates[0]]
            prms_dict['const1'] = [mdl_prms_list[0]]
            prms_dict['beta0']  = [mdl_prms_list[2]]
            prms_dict['beta1']  = [mdl_prms_list[3]]
            prms_dict['beta2']  = [mdl_prms_list[4]]

        elif (tipo_modello == 'CIR'):
            
            prms_dict['Dates']  = [bench_dates[0]]
            prms_dict['r0']     = [mdl_prms_list[0]]
            prms_dict['kappa']  = [mdl_prms_list[1]]
            prms_dict['theta']  = [mdl_prms_list[2]]
            prms_dict['sigma']  = [mdl_prms_list[3]]
            
        else:
            
            print 'MODELLO NON DISOPONIBILE!!!'


        return prms_dict





def  fitting_with_plot(c_dates, c_values, opt_dict, flag_plot):
    
    opt_dict['MakeGraph'] = flag_plot
    #opt_dict['MakeGraph'] = True
    
    t_mkt = []
    for i in range(0, len(c_dates)):
        
        dateDiff = c_dates[i] - c_dates[0]
        timeTmp = float(dateDiff.days/365.2425)
        t_mkt.append(timeTmp)
    
    model_type = opt_dict['interp']
    
    
    #mthod_o ='SLSQP' #ok
    mthod_o ='TNC' #ok

    zc_mkt = c_values
    #zc_mkt = data_raw['ValoreZC']
    
    #t_mkt  = data_raw['MatTimes']
    t_mkt  = np.asarray(t_mkt)


    
    if (model_type == '0'): # caso lineare

        prms_lin = estimate_linear_params(t_mkt, zc_mkt)
        
    elif (model_type == '1'): # caso AVD
        
        prms_avd = estimate_avd_params(t_mkt, zc_mkt)
        
    else:

        x_bnd = []
        x0 = []
        
        if (model_type == '2'): # caso SVE
            
            ln_prms = len(opt_dict['bound_min_sve'])

            for i in range(0, ln_prms):
                
                x_bnd_minTmp = opt_dict['bound_min_sve'][i]
                x_bnd_maxTmp = opt_dict['bound_max_sve'][i]
                x0Tmp        = float((x_bnd_maxTmp + x_bnd_minTmp)/2.0)
        
                x_bnd.append([x_bnd_minTmp, x_bnd_maxTmp])
                x0.append(x0Tmp)


        if (model_type == '3'): # caso CIR
        
            ln_prms = len(opt_dict['bound_min_cir'])
            
            for i in range(0, ln_prms):
                
                x_bnd_minTmp = opt_dict['bound_min_cir'][i]
                x_bnd_maxTmp = opt_dict['bound_max_cir'][i]
                x0Tmp        = float((x_bnd_maxTmp + x_bnd_minTmp)/2.0)
                
                x_bnd.append([x_bnd_minTmp, x_bnd_maxTmp])
                x0.append(x0Tmp)
        
        if (model_type == '4'): # caso NS
        
            ln_prms = len(opt_dict['bound_min_ns'])
            
            for i in range(0, ln_prms):
                
                x_bnd_minTmp = opt_dict['bound_min_ns'][i]
                x_bnd_maxTmp = opt_dict['bound_max_ns'][i]
                x0Tmp        = float((x_bnd_maxTmp + x_bnd_minTmp)/2.0)
                
                x_bnd.append([x_bnd_minTmp, x_bnd_maxTmp])
                x0.append(x0Tmp)

            
        res = optimize.minimize(loss_fun_for_fitting, x0, method = mthod_o,  args=(t_mkt, zc_mkt, model_type), bounds = x_bnd)

    dict_model_par = {}
    
    """
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
    """
    

    if (model_type == '0'): # Linear

        date_prms = []
        a_prms = []
        b_prms = []
        
        values_list = prms_lin['a'].keys()
        ln = len(values_list)
        
        for i in range(0, ln):
        
            
            #dateTmp = data_raw['MatDate'][i]
            dateTmp =c_dates[i]

            a_prmsTmp = prms_lin['a'][i]
            b_prmsTmp = prms_lin['b'][i]

            a_prms.append(a_prmsTmp)
            b_prms.append(b_prmsTmp)
            date_prms.append(dateTmp)


        dict_model_par['a']      = a_prms
        dict_model_par['b']      = b_prms
        dict_model_par['Dates']  = date_prms

    elif (model_type == '1'): # AVD

        date_prms = []
        a_prms = []
        b_prms = []
        c_prms = []
        d_prms = []
        e_prms = []

        values_list = prms_avd['a'].keys()
        ln = len(values_list)
        
        for i in range(0, ln):
        
            #dateTmp = data_raw['MatDate'][i]
            dateTmp = c_dates[i]

            a_prmsTmp = prms_avd['a'][i]
            b_prmsTmp = prms_avd['b'][i]
            c_prmsTmp = prms_avd['c'][i]
            d_prmsTmp = prms_avd['d'][i]
            e_prmsTmp = prms_avd['e'][i]

            date_prms.append(dateTmp)

            a_prms.append(a_prmsTmp)
            b_prms.append(b_prmsTmp)
            c_prms.append(c_prmsTmp)
            d_prms.append(d_prmsTmp)
            e_prms.append(e_prmsTmp)

        dict_model_par['a']      = a_prms
        dict_model_par['b']      = b_prms
        dict_model_par['c']      = c_prms
        dict_model_par['d']      = d_prms
        dict_model_par['e']      = e_prms

        dict_model_par['Dates']  = date_prms

    elif (model_type == '2'): # Svensson

        #dataRef = data_raw['MatDate'][0]
        dataRef = c_dates[0]

        dict_model_par['Dates']  = [dataRef]
        dict_model_par['const1'] = [res.x[0]]
        dict_model_par['const2'] = [res.x[1]]
        dict_model_par['beta0']  = [res.x[2]]
        dict_model_par['beta1']  = [res.x[3]]
        dict_model_par['beta2']  = [res.x[4]]
        dict_model_par['beta3']  = [res.x[5]]
        
    elif (model_type == '3'): # CIR

        #dataRef = data_raw['MatDate'][0]
        dataRef = c_dates[0]

        dict_model_par['Dates'] = [dataRef]
        dict_model_par['r0']    = [res.x[0]]
        dict_model_par['kappa'] = [res.x[1]]
        dict_model_par['theta'] = [res.x[2]]
        dict_model_par['sigma'] = [res.x[3]]

    else:

        #dataRef = data_raw['MatDate'][0]
        dataRef = c_dates[0]
        
        dict_model_par['r0']    = [res.x[0]]
        dict_model_par['Dates'] = [dataRef ]      
        dict_model_par['kappa'] = [res.x[1]]
        dict_model_par['theta'] = [res.x[2]]
        dict_model_par['sigma'] = [res.x[3]]
    
    #--------------- MAKE GRAPH -------------------
    make_graph = opt_dict['MakeGraph']
    
    if (make_graph == True):

        dataRef = c_dates[0]

        model_type = opt_dict['interp']
        
        mkt_times  = t_mkt
        mkt_values = c_values
        
        chk_graph(mkt_times, mkt_values, model_type, dict_model_par)


    return dict_model_par

def convertNodeToMnth(nodeRef):

    dict_N2T = {} 
    dict_N2T['1M'] = 1 
    dict_N2T['3M'] = 3
    dict_N2T['6M'] = 6
    dict_N2T['1Y'] = 12
    
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

        elif (ts_tmp == 'S') or (ts_tmp == '3G'):      

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

    ln = len(x_list_ref) 

    for i in range(0, ln):

        if (x_target < x_list_ref[0]): 
            i_ref = None
            break

        elif (x_target > x_list_ref[ln-1]) : 
            i_ref = None
            break

         
        elif (x_list_ref[i] == x_target) and (i == 0): 
            i_ref = i
            break

        
        elif (x_list_ref[i] >= x_target):
            i_ref = i - 1
            break
            
        else:

            continue

    return i_ref


def find_indx_next(x_target, x_list_ref):

    ln = len(x_list_ref) 

    for i in range(0, ln):

        if (x_target < x_list_ref[0]): 
            i_ref = 0
            break

        elif (x_target > x_list_ref[ln-1]) : 
            i_ref = ln -1
            break

         
        elif (x_list_ref[i] == x_target) and (i == 0): 
            i_ref = i
            break

        
        elif (x_list_ref[i] >= x_target):
            i_ref = i
            break
            
        else:

            continue

    return i_ref


def find_indx_toll(x_target, x_list_ref, toll):

    ln = len(x_list_ref) - 1

    if (x_target < x_list_ref[0]): 
        i_ref = 0

    elif (x_target > x_list_ref[ln]): 
        i_ref = ln

    else:

        for i in range(0, len(x_list_ref)):
    
            if (abs(x_list_ref[i] - x_target) < toll ): 
                i_ref = i
                break
    
            elif (x_target < x_list_ref[i] - toll/2.0): 
                i_ref = i
    
            else:
    
                continue

    return i_ref



def find_indx_n(x_target, x_list_ref):

    ln = len(x_list_ref) 

    for i in range(0, ln):

        if (x_list_ref[0] > x_target): 
            i_ref = None
            break

        elif (x_list_ref[ln-1] < x_target): 
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




def interp_df(t_s_target, df_n_s, df_o_s, t_n_s, t_o_s, flag_interp1):
        
    if (flag_interp1 == 1):  #%----------intepolazione esponenziale sui fattori di sconto ---------------
    
        futures_start_df_0 =  df_from_interp_df_exp(t_s_target, df_n_s, df_o_s, t_n_s, t_o_s)
    
    else: #%-------------- interpolazione lineare sui tassi -----------------------
        
        futures_start_df_0 = df_from_interp_R_lin(t_s_target, df_n_s, df_o_s, t_n_s, t_o_s)

    return futures_start_df_0


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



def compute_df_future(seg1_times, seg1_val, seg1_df,  futures_rates, flag_futures_gap, futures_start_time, futures_end_time,  futures_end_df_x, flag_interp1, indx_ref):

    
    t_s_target   = futures_start_time[indx_ref]
    t_e_target   = futures_end_time[indx_ref]
    
    indx_seg1_last_s = find_indx(t_s_target, seg1_times)
    indx_seg1_last_e = find_indx(t_e_target, seg1_times)


    
    index_cross     = find_indx_equal(futures_start_time[indx_ref], seg1_times)
    

    if (indx_ref != 0):
        df_o = futures_end_df_x[indx_ref-1]
        
    t_o  = futures_end_time[indx_ref-1]
    
    futures_rates_old = futures_rates[indx_ref - 1]

    
    
    if (index_cross != None): # caso sovrapposizione
        futures_start_1 = seg1_df[index_cross]

        futures_start_df_0 = futures_start_1
        futures_rates_target = futures_rates[indx_ref]

        futures_end_df_0 = futures_start_1/(1.0 + futures_rates_target*(futures_end_time[indx_ref] - t_s_target) )
        
        

    elif (t_s_target < seg1_times[len(seg1_times)-1]): #% controllo della sovrapposizione: seg1 e futures



        df_o_s = seg1_df[indx_seg1_last_s]        
        df_n_s = seg1_df[indx_seg1_last_s + 1] 

        t_o_s = seg1_times[indx_seg1_last_s]        
        t_n_s = seg1_times[indx_seg1_last_s + 1] 


        futures_start_df_0 = interp_df(t_s_target, df_n_s, df_o_s, t_n_s, t_o_s, flag_interp1)
        
        if (t_e_target < seg1_times[len(seg1_times)-1]):

            df_o_e = seg1_df[indx_seg1_last_e]        
            df_n_e = seg1_df[indx_seg1_last_e + 1] 

            t_o_e = seg1_times[indx_seg1_last_e]        
            t_n_e = seg1_times[indx_seg1_last_e + 1] 

            futures_end_df_0 = interp_df(t_e_target, df_n_e, df_o_e, t_n_e, t_o_e, flag_interp1)
            
        else:
            
            futures_rates_target = futures_rates[indx_ref]                    
            futures_end_df_0 = futures_start_df_0/(1.0 + futures_rates_target*(futures_end_time[indx_ref] - t_s_target))
    

    else:  # ---- NO SOVRAPPOSIZIONE

        #t_o  = futures_start_time[indx_ref - 1]
        t_o  = futures_end_time[indx_ref - 1]
        df_o = futures_end_df_x[indx_ref-1]

        #--------------------------

        futures_rates_target = futures_rates[indx_ref]
        
        if (indx_ref == 0) and (len(seg1_val) == 0):
                
            R_start1 = futures_rates_target
            futures_start_df_0 = np.exp(-R_start1*futures_start_time[indx_ref])
        
        elif (int(flag_futures_gap) == 1) and (len(seg1_val) > 0): #------------ gap colmato mantenedo costante il tasso spot


            if (indx_ref == 0):
                 
                ln_s = len(seg1_df)-1
                df_o = seg1_df[ln_s]
                t_o   = seg1_times[ln_s]


            R_start1            = -np.log(df_o)/t_o #CHK indici
            
            
            futures_start_df_0 =  np.exp(-R_start1*t_s_target)
            
        else:  #------------ gap colmato mantenedo costante il tasso fwd ---------------------------------------------
            
            
            if (indx_ref == 0):
                
                
                ln_s = len(seg1_times)-1
                t_o   = seg1_times[ln_s]
                t_oo  = seg1_times[ln_s-1]
        
                df_o = seg1_df[ln_s]
                df_oo = seg1_df[ln_s-1]
                
                fwd_start1         = -np.log(df_o/df_oo)/(t_o - t_oo)
                futures_start_df_0 = df_o*np.exp(-fwd_start1*(t_s_target - t_o ))
                

            else:
                fwd_start1         = futures_rates_old
                futures_start_df_0 = df_o/(1.0 + fwd_start1*(t_s_target - t_o ))
            
        


        futures_end_df_0 = futures_start_df_0/(1 + futures_rates_target*(t_e_target - t_s_target) )
        #dr = 1.0/(t_e_target - t_s_target)*(futures_start_df_0/futures_end_df_0 - 1.0)



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

    f0 = 0
    for k in range(n,m):    
        
        f0 = (time_new[k] - time_new[k-1])*np.exp(-x*(time_new[k]-time_new[n])) + f0
    
    f = f1 - x
    
    return f



def fwd_fun_0(x, n, m, f1, py_rate_m, time_new):

    #--------------------------------------------------------------------------
    f0 = 0
    
    for k in range(n+1,m+1):    
        
        dt_tmp = (time_new[k] - time_new[k-1])
        n_tmp = k-n
        z_Tmp = (1.0 + x*dt_tmp)**-(n_tmp)
        f0 = dt_tmp*z_Tmp + f0

    X = py_rate_m*f0
    Y = (1.0 + x*dt_tmp)**-(m-n)
    
    f_out = X  + Y - f1

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
        
            if (valRefTmp in val_not_allowed):
                continue
            else:
                new_vec.append(data_out[kk][i])
         
        data_out_n[kk] = new_vec
    
    return data_out_n
        



def graphrates(dep_times, dep_rates, fu_times, fu_rates, sw_times, sw_rates, boot_tims, boot_rates):


    import matplotlib
    matplotlib.use('TkAgg')
    
    from numpy import arange, sin, pi
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
    from matplotlib.backend_bases import key_press_handler
    
    
    from matplotlib.figure import Figure
    
    import Tkinter as Tk
    
    root = Tk.Tk()
    root.wm_title("Plot Bond fitting results")
    
    f = Figure(figsize=(5, 4), dpi=100)
    a = f.add_subplot(111)

    a.plot(dep_times, dep_rates, '-o', label='depositi')

    if (fu_times[0] != 999): a.plot(fu_times, fu_rates, '-o', label='futures')   
    if (sw_times[0] != 999): a.plot(sw_times, sw_rates,'-o', label='swap')
    
    a.plot(boot_tims, boot_rates,'--o', label='bootstrap')
    
    a.set_title('Bootstrap curve')
    a.set_xlabel('Tempo [anni]')
    a.set_ylabel('Livello tassi')
    legend = a.legend(loc='upper left', shadow=False)
    a.grid()
    
    #------------------------------------------------------------
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
    
    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

    def _quit():
        root.quit()     # stops mainloop
        root.destroy()  # this is necessary on Windows to prevent

    def on_key_event(event):
        print('you pressed %s' % event.key)
        key_press_handler(event, canvas, toolbar)

    canvas.mpl_connect('key_press_event', on_key_event)

    button = Tk.Button(master=root, text='Quit', command=_quit)
    button.pack(side=Tk.BOTTOM)
    
    Tk.mainloop()


    
    """
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
    """
    



def graphdf(dep_times, dep_df, fu_times, fu_df, sw_times, sw_df):


    import matplotlib
    matplotlib.use('TkAgg')
    
    from numpy import arange, sin, pi
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
    from matplotlib.backend_bases import key_press_handler
    
    
    from matplotlib.figure import Figure
    
    import Tkinter as Tk
    
    root = Tk.Tk()
    root.wm_title("Plot Bond fitting results")
    
    f = Figure(figsize=(5, 4), dpi=100)
    a = f.add_subplot(111)



    #g1 = plt.figure()

    a.plot(dep_times, dep_df, 'o-', label='depositi')

    if (fu_times[0] != 999): a.plot(fu_times, fu_df, 'o-', label='futures')
    if (sw_times[0] != 999): a.plot(sw_times, sw_df,'o-', label='swap')
    
    a.set_title('Bootstrap curve')
    a.set_xlabel('Tempo [anni]')
    a.set_ylabel('Livello fattore di sconto')
    legend = a.legend(loc='upper right', shadow=False)
    a.grid()
    
    """
    plt.xlabel('Tempi [anni]')
    plt.ylabel('Livello fattore di sconto')
    plt.grid()
    plt.legend()
    plt.show()
    return g1
    """

    #------------------------------------------------------------
    #------------------------------------------------------------

    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
    
    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

    def _quit():
        root.quit()     # stops mainloop
        root.destroy()  # this is necessary on Windows to prevent

    def on_key_event(event):
        print('you pressed %s' % event.key)
        key_press_handler(event, canvas, toolbar)

    canvas.mpl_connect('key_press_event', on_key_event)

    button = Tk.Button(master=root, text='Quit', command=_quit)
    button.pack(side=Tk.BOTTOM)
    
    Tk.mainloop()







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



def chkDataCoherence(dict1s, dict_f, dict_s):
    
    s1 = list(set(dict1s['TipoSegmento']))
    
    try:
        s1.remove('3G')
    except:
        
        pass
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

        '--------------- save graph ------------------------------------------------------------'

    flag_save_graph = 0
    if (flag_save_graph == 1):
    
        print 'AAAAAA'
        #g1.savefig('fig/plot_tassi_%s.png' %(ref_date))
        #g2.savefig('fig/plot_discount_%s.png' %(ref_date))
 
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


def set_data_opt(refDate):
    
    
    
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

    data_opt['RefDate'] = refDate
    #data_opt['RefDate'] = datetime.date(2017, 02, 15)

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
    opt_dict['hr_bootMethod']  = 0 #1 = LCS, 0 = CHR
    opt_dict['Basis']          = 'ACT/360' 
    opt_dict['BusConv']        = 'follow'


    opt_dict['opt_path_graph']  =  'C:\\'

    #data_opt['hr_bootMethod']  = opt_boot_meth #0 = LCS, 1 = CHR
    opt_dict['bench_interp']   ='2' #
    opt_dict['hr_interp']      ='2'
    

    return opt_dict


def set_data_opt_for_fit(model_fit, flag_make_graph, data_ref):
    
    
    opt_dict = {}
    
    opt_dict['MakeGraph'] = flag_make_graph 
    opt_dict['interp']    = model_fit 
    opt_dict['DataRef']   = data_ref

    
    """
    opt_dict['interp'] == '0' # 'LIN'
    opt_dict['interp'] == '1' # 'AVD'
    opt_dict['interp'] == '2' # 'SVE'
    opt_dict['interp'] == '3' # 'CIR'
    opt_dict['interp'] == '4' # 'NS'
    """

    opt_dict['opt_fwd_tenor']   = "1M"
    
    """ 
    opt_dict['opt_fwd_tenor']   = "3M"
    opt_dict['opt_fwd_tenor']   = "6M"
    opt_dict['opt_fwd_tenor']   = "12M"
    """
    
    opt_dict['opt_path_graph']  =  'C:\\'
    
    opt_dict['fit_type']        = "boot"
    #opt_dict['fit_type']        = "py"
    

    
    if (opt_dict['interp'] == '2'): #->SVE
    
        #bound_min_sve = [0.0001,  0.0001, -10.00, -10.050, -10.00, -10.00]

        bound_min_sve = [0.0001,  0.0001, -10.00, -10.00, -10.00, -10.00]
        bound_max_sve = [10.0,     50.0,  10.03,   10.0,   10.5,   10.0]
        bound_x0_sve  = [1.0,       10.0,   0.03,   0.03,   0.03,   0.03]

        opt_dict['bound_min_sve'] = bound_min_sve
        opt_dict['bound_max_sve'] = bound_max_sve
        opt_dict['x0']            = bound_x0_sve

    elif (opt_dict['interp'] == '3'): #->CIR

        bound_min_cir = [  -0.1,  0.1, 0.001, 0.001]
        bound_max_cir = [ 10.00, 10.0, 10.00, 1.000]
        bound_x0_cir  = [0.0001,  1.0, 0.015,  0.01]

        opt_dict['bound_min_cir'] = bound_min_cir
        opt_dict['bound_max_cir'] = bound_max_cir
        opt_dict['x0']            = bound_x0_cir

    elif (opt_dict['interp'] == '4'): #->NS

        bound_min_ns = [0.0001,  -10.0, -10.0, -10.0]
        bound_max_ns = [100,     +10.0, +10.0, +10.0]
        bound_x0_ns  = [5.0001,  0.03, 0.03, 0.03]

        opt_dict['bound_min_ns'] = bound_min_ns
        opt_dict['bound_max_ns'] = bound_max_ns
        opt_dict['x0']            = bound_x0_ns


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


def computeTimesFromDates(dateList):

    matTimes = []
    dateRef  = dateList[0]
    
    ln = len(dateList)
    
    matTimes.append(0.0)
    
    for i in range(1, ln):
        
        timeTmp = (dateList[i] - dateRef).days/365.2425

        matTimes.append(timeTmp)

    return matTimes


def fromCurveToSpread(df_bench_values, zc_bench_dates, prms_bench, bench_model, py_risky_val, py_risky_dates, risky_model, targetDates, targetTimes):
    
    
    dataRef = zc_bench_dates[0]
    flagPlot = False
    
    fitting_opt_dict = set_data_opt_for_fit(risky_model, flagPlot, dataRef)
    
    #fitting_opt_dict = set_opt_for_fit(risky_model, dataRef)    
    zc_bench_times   = fromDates2Times(zc_bench_dates)
    py_risky_times   = fromDates2Times(py_risky_dates)

    refDates = zc_bench_dates[0]
    
    prms_risky_for_py      = fitting(py_risky_dates, py_risky_val, fitting_opt_dict)

    py_risky_val_fitted    = makeRatesFromModel(py_risky_times, py_risky_val, prms_risky_for_py, targetTimes, risky_model)
    zc_risk_free           = makeRatesFromModel(zc_bench_times, df_bench_values, prms_bench, targetTimes, bench_model)
    
    pyFreq = 0.25
    py_risk_free = computePYRates(zc_bench_times, df_bench_values, pyFreq, 0.0, targetTimes)
    
    py_spread = py_risky_val_fitted - py_risk_free
    
    
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
    
    data_opt_for_boot = set_data_opt(dataRef)
    data_opt_for_boot['RefDate'] = dataRef
    data_opt_for_boot['SwapGapMethod'] = 1
    
    #n_rates = len(py_risky_val_n)
    
    target_dates = fromTimes2Dates(refDates, targetTimes)

    #data_raw_risky = buildDataRawForBoot(target_dates, targetTimes, py_risky_val_n)
    
    
    #--------------- IMPLEMENTA LA PARTE DI BOOT ----------------------
    
    freq_pay = 0.5
    
    time_ref, zc_risky_val_n = formPytoZC(targetTimes, py_risky_val_n, freq_pay)
    

    #zc_risky_val_n = py_risky_val_n
        
    zc_spread_n    = zc_risky_val_n - zc_risk_free
    
        
        
    return py_spread_n, zc_spread_n

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


    return surv, hr
    




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
    k =  find_indx_next(maturities_cds[0], times)

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
        m =  find_indx_next(times[i], maturities_cds)

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
            j =  find_indx_next(maturities_cds[m], times)
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

def func_h2(hr,cds_rr,t,z,cds_spread,delta_t):

    dummy1 = z[1:]*np.exp(-hr*(t[0:-1]-t[0]))*(1-np.exp(-hr*(t[1:]-t[0:-1])))
    dummy2 = z[1:]*delta_t[1:]*np.exp(-hr*(t[1:]-t[0]))
    
    f = cds_spread - (1-cds_rr) * sum(dummy1)/sum(dummy2)
    
    return f

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

def func_chr(hr,b,cds_rr,s_n,t,z,cds_m,delta_t):
    
    dummy1 = 0
    dummy2 = 0
    
    for i in range(0,len(t)-1):
        dummy1 = dummy1 + ( np.exp(-hr*(t[i]-t[0])) - np.exp(-hr*(t[i+1]-t[0])) )*z[i]  
        dummy2 = dummy2 + delta_t[i]*z[i]*np.exp(-hr*(t[i+1]-t[0]))
    
    f = b - (1-cds_rr)*s_n*dummy1 + s_n*cds_m*dummy2
    
    return f


def fromTimes2Dates(refDates, times_list):

    target_dates = []
    
    n_times = len(times_list)
    
    
    for i in range(0, n_times):
    
        timeTmp = times_list[i]
        daysTmp = int(365.2425*timeTmp)
        
        target_datesTmp = refDates  + datetime.timedelta(days=daysTmp)
        target_dates.append(target_datesTmp)
        
    return target_dates

def mapStr2Time(str_in):

    
    unit   = str_in[-1:]
    number = int(str_in[:-1])

    if (unit == 'M'):

        time_out = float(number*1.0/12.0)
    
    elif (unit == 'W'):

        time_out = float(number*1.0/52.0)

    elif (unit == 'G'):
    
        time_out = float(number/365.2425)

    else:
    
        print 'Unita non riconosciuta'
        
    return time_out



def mapTime2Str(time_in):


    n_mnth = time_in*12.0
    r_mnth = n_mnth%1.0

    n_week = time_in*52.0
    r_week = n_week%1.0
    
    n_gg = time_in*365
    
    if (n_mnth >= 1.0):
    
        if (r_mnth > 0.5):
        
            n_mnth = int(n_mnth) + 1
        else:
            n_mnth = int(n_mnth)

        tag = str(n_mnth) + 'M'
    else:

        if (n_week >= 1):
    
            if (r_week > 0.5):
                n_week = int(n_week) + 1
            else:
                n_week = int(n_week)
                
            tag = str(n_week) + 'W'

        else:

            n_gg = int(n_gg)
                
            tag = str(n_gg) + 'G'

        
    return tag


def    formPytoZC(T, py_values, pay_rate):


    ln = len(T)
    T_max = T[ln-1]
    n_step = int (T_max/pay_rate) +1
    
    T_new = []
    T_new.append(0.0)
    
    T_tmp = 0.0
    for j in range(1, n_step):
        T_tmp = T_tmp + pay_rate 
        T_new.append(T_tmp)
    

    pyNew = np.interp(T_new, T, py_values)

    Z   = []

    ZDict   = {}
    rZcDict = {}

    Z.append(1.0)
    Z1 = (1.0/(1.0 + pyNew[1]))**pay_rate
    Z.append(Z1)
    
    pay_rate_str = mapTime2Str(pay_rate)

    ZDict['0G']     = 1.0
    ZDict[pay_rate_str]   = Z1

    rZcDict['0G']   = pyNew[1]
    rZcDict[pay_rate_str] = pyNew[1]
    
        
    for i in range(2, len(pyNew)):
    
        
        D_j = 0.0

        for j in range(1, i):

            D_tmp = pay_rate*Z[j]
            D_j = D_j + D_tmp
        
        T_i    = T_new[i]
        sr_i   = pyNew[i]

        Z_i = (1.0 - D_j*sr_i)/(1.0 + sr_i*pay_rate)
        Z.append(Z_i)

        rZc_tmp = (1.0/Z_i)**(1.0/T_i) - 1.0 
        
        T_i_str = mapTime2Str(T_i)
        
        rZcDict[T_i_str] = rZc_tmp
        ZDict[T_i_str]   = Z_i

    rZcDictNew = {}
    
    zr_out = []
    time_out = []
    
    for j in range(0, len(py_values)):
    
        T_i = T[j]
        T_i_str = mapTime2Str(T_i)
        
        time_out.append(T_i)
    
        if (T_i <= 0.5):
            rZcDictNew[T_i_str] = py_values[j]
            
        else:

            rZcDictNew[T_i_str] = rZcDict[T_i_str]
        
        zr_out.append(float(rZcDictNew[T_i_str]))
        
    
    time_out = np.asarray(time_out)
    zr_out = np.asarray(zr_out)
    
    return time_out, zr_out
    
    
    


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
    
    #model_risky  = bench_data['Model']
    model_rf     = bench_data['Model']
    model_risky  = int(opt_dict['hr_interp'])

    boot_method   = opt_dict['hr_bootMethod']
    #model_risky   = '0'
    
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
    


    boot_method = 1
    
    if (boot_method == 0):

                        
        surv_prob_all, hr_values_all = bootstrap_lcs(cds_times_n, cds_val_n, cds_pay_times, RR, df_risk_free)

        log_surv = np.log(surv_prob_all)

        hr_values   = np.interp(cds_times_n, cds_pay_times, hr_values_all)
        log_surv    = np.interp(cds_times_n, cds_pay_times, log_surv)
        
        surv_prob = np.exp(log_surv) 

        """
        fileOutSurv = 'dataSurv_lcs_v0.txt'
        fileOutHR   = 'dataaH_lcs_v0.txt'
        
        dumpDataOnFile(cds_times_n, surv_prob, fileOutSurv)
        dumpDataOnFile(cds_times_n, hr_values, fileOutHR)
        """
            
    elif (boot_method == 1):
        
        
        cds_pay_times = np.array(cds_pay_times)

        #df_risk_free = np.exp(-cds_pay_times*0.01)

        ln = len(cds_pay_times)
        surv_prob_all, hr_values_all, times_n = bootstrap_chr(cds_times, cds_values, RR, cds_pay_times, df_risk_free)
        
        hr_values   = np.interp(cds_times_n, times_n, hr_values_all)
        surv_prob   =  np.interp(cds_times_n, times_n, surv_prob_all)

        
        """
        fileOutHR         = 'output_test/dataH_chr_v2.txt'
        fileOutSurv       = 'output_test/dataSurv_chr_v2.txt'
        fileOutZCSpread   = 'output_test/dataSpread_zc_v2.txt'
        fileOutPYSpread   = 'output_test/dataSpread_py_v2.txt'
        fileOutMDefault   = 'output_test/dataMDefault_py_v2.txt'
        """


        """     
        fileOutHR         = 'output_test/dataH_chr_flat_v0.txt'
        fileOutSurv       = 'output_test/dataSurv_chr_flat_v0.txt'
        fileOutZCSpread   = 'output_test/dataSpread_zc_flat_v0.txt'
        fileOutPYSpread   = 'output_test/dataSpread_py_flat_v0.txt'
        fileOutMDefault   = 'output_test/dataMDefault_py_flat_v0.txt'
        """
        
        """
        dumpDataOnFile(cds_times_n, surv_prob, fileOutSurv)
        dumpDataOnFile(cds_times_n, hr_values, fileOutHR)
        dumpDataOnFile(cds_times_n, zcSpread, fileOutZCSpread)
        dumpDataOnFile(cds_times_n, pySpread, fileOutPYSpread)
        """
        
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




