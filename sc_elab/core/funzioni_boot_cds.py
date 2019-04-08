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
from datetime import timedelta

import funzioni_base as fb


def set_data_opt_for_cds(data_ref):

    # bench_boot_method 
    # interp_method 
    
    opt_dict = {}

    #opt_dict['DataRef']        = datetime.date(2005, 12, 31)
    opt_dict['DataRef']        = data_ref

    opt_dict['MKT']            = 'de'
    opt_dict['tenor']          = 3 # 
    opt_dict['Basis']          = 'ACT/365' #ACT/365
    opt_dict['interp']         = '0' # 
    opt_dict['BusConv']        = 'modfollow'
    opt_dict['fixingDays']     = 2
    opt_dict['compounding']    = 0   #0 = semplice, 1 = composto, 2 = continuo
    opt_dict['ReocveryRate']   = 0.4
    opt_dict['hr_bootMethod']  = 0 #1 = LCS, 0 = CHR
    opt_dict['Basis']          = 'ACT/360' 
    opt_dict['BusConv']        = 'follow'


    opt_dict['opt_path_graph']  =  'C:\\'

    #data_opt['hr_bootMethod']  = opt_boot_meth #0 = LCS, 1 = CHR
    opt_dict['bench_interp']   ='0' #
    opt_dict['hr_interp']      ='0'
    

    return opt_dict


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
    


def func_h2(hr,cds_rr,t,z,cds_spread,delta_t):

    dummy1 = z[1:]*np.exp(-hr*(t[0:-1]-t[0]))*(1-np.exp(-hr*(t[1:]-t[0:-1])))
    dummy2 = z[1:]*delta_t[1:]*np.exp(-hr*(t[1:]-t[0]))
    
    f = cds_spread - (1-cds_rr) * sum(dummy1)/sum(dummy2)
    
    return f


def bootstrap_chr(maturities_cds, adjusted_cds, recovery_rate, times, discount_factors, dt):
                                            
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
    #dt     = 0.25

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

def func_chr(hr,b,cds_rr,s_n,t,z,cds_m,delta_t):
    
    dummy1 = 0
    dummy2 = 0
    
    for i in range(0,len(t)-1):
        dummy1 = dummy1 + ( np.exp(-hr*(t[i]-t[0])) - np.exp(-hr*(t[i+1]-t[0])) )*z[i]  
        dummy2 = dummy2 + delta_t[i]*z[i]*np.exp(-hr*(t[i+1]-t[0]))
    
    f = b - (1-cds_rr)*s_n*dummy1 + s_n*cds_m*dummy2
    
    return f



def  plotCDS(cds_times, cds_values,timeOut, pySpread, cds_title):

    import matplotlib
    matplotlib.use('TkAgg')
    
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
    from matplotlib.backend_bases import key_press_handler
    
    
    from matplotlib.figure import Figure
    
    import Tkinter as Tk
    
    root = Tk.Tk()
    root.wm_title("Plot Bond fitting results")
    
    f = Figure(figsize=(5, 4), dpi=100)
    a = f.add_subplot(111)

    
    a.plot(cds_times, cds_values,'o', label='CDS Mkt')
    a.plot(timeOut, pySpread,'-.k', label='Py Spread')
    
    a.set_title('CDS bootstrap results: %s'%cds_title)
    a.set_xlabel('Tempo [anni]')
    a.set_ylabel('Livello CDS')
    a.legend(loc='lower right', shadow=False)
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



def boot_cds(opt_dict, raw_data, bench_data, swap_data):
    

    output_data = {}
    
    #--> definisco struttura output: lista date
    #outputDates
    
    dateOut, timeOut = computeDateOutputCDS(raw_data, opt_dict)
    benchTimes       = fb.computeTimesFromDates(bench_data['MatDate'])
    
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
    
    py_bench_rates = fb.computePYRates(benchTimes, benchDf, cds_pay_freq, time_0, timeOut)

    
    benchDf = np.asarray(benchDf)
    benchTimes = np.asarray(benchTimes)

    benchZcRates = -1.0/benchTimes*np.log(benchDf)
    benchZcRates[0] = benchZcRates[1]
    
    tipo_curva_bench = bench_data['Type']
    
    if(tipo_curva_bench != 'swap'):

        swapTimes_risk_free = fb.computeTimesFromDates(swap_data['MatDate'])
        swapDf_risk_free    = swap_data['DiscountFactor']

        py_swap_rates   = fb.computePYRates(swapTimes_risk_free, swapDf_risk_free, pyFreq, time_0, timeOut)

        py_spread_adj   = py_swap_rates - py_bench_rates
        cds_adj         = cds_val_n + py_spread_adj
        cds_val_n       = cds_adj
        

    py_risky = cds_val_n + py_bench_rates
    

    py_risky_times = timeOut
    py_risky_dates = dateOut
    
    
    py_risky = np.insert(py_risky, 0, py_risky[0])
    py_risky_dates = np.insert(py_risky_dates, 0, py_risky_dates[0])
    py_risky_dates[1] = py_risky_dates[0] +  timedelta(days=1)
    

    
    pySpread, zcSpread = fb.fromCurveToSpread(benchDf, zcBenchDates, prmsBench, model_rf, py_risky, py_risky_dates, model_risky, dateOut, timeOut, cds_pay_freq)
    

      
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
        cds_pay_dateTmp = busdayrule.rolldate(cds_pay_dateTmp, mkt_ref, day_conv_cds)

        cds_pay_dates.append(cds_pay_dateTmp)
        
    cds_pay_times = daycount.yearfractions(cds_pay_dates, basis_cds)
    cds_pay_times = np.asarray(cds_pay_times)

    
    
    benchZcRates_n = np.interp(cds_pay_times, benchTimes, benchZcRates)
    benchDf_n = np.exp(-benchZcRates_n*cds_pay_times)
    df_risk_free = benchDf_n
    
    #dt = 0.25
    times_n = []
    t_last = np.around(cds_pay_times[len(cds_pay_times)-1], decimals = 0)
    #n_t = int(t_last/dt)

    n_t = int(t_last/cds_pay_freq)    
    for i in range(0, n_t+1):
        
        tTmp = cds_pay_freq*i
        times_n.append(tTmp)
        
    
    log_z        = -np.log(df_risk_free[1:])/cds_pay_times[1:]
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

            
    elif (boot_method == 1):
        
        
        cds_pay_times = np.array(cds_pay_times)

        #df_risk_free = np.exp(-cds_pay_times*0.01)

        ln = len(cds_pay_times)
        surv_prob_all, hr_values_all, times_n = bootstrap_chr(cds_times, cds_values, RR, cds_pay_times, df_risk_free, cds_pay_freq)
        
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
    
    cds_title = opt_dict['Emittente']

    plotCDS(cds_times, cds_values,timeOut, pySpread, cds_title) 

    #------------------------------------------
    #------ set up Output----------------------
    #------------------------------------------

    output_data['outputTimes']     = timeOut
    output_data['outputDates']     = dateOut
    output_data['zcSpread']        = zcSpread
    output_data['pySpread']        = pySpread
    output_data['hazardRate']      = hr_values
    output_data['survProbCum']     = surv_prob
    output_data['marginalDefault'] = marg_d_n
    
    return output_data     
    




