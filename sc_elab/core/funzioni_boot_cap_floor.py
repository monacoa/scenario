import numpy as np
import pandas as pd
from scipy.optimize import root
from scipy.stats import norm

# -------------- Valore atteso di Black con shift ------------------------------
def BlackExpectedValue(fwdPrice,strike,TimeToReset,vol,shift,type):
    if TimeToReset == 0.0:
        diff = fwdPrice - strike
        if type=='Call':
            ExpectedValue=np.maximum(0.0,diff)
        elif type=='Put':
            ExpectedValue=np.maximum(0.0,-diff)
        else:
            ExpectedValue = None
            print('tipo opzione non valido')
    else:
        d1=(np.log((fwdPrice+shift)/(strike+shift))+(vol*vol/2)*TimeToReset)/(vol*np.sqrt(TimeToReset))
        d2=d1-vol*np.sqrt(TimeToReset)
        if type=='Call':
            ExpectedValue=(fwdPrice+shift)*norm.cdf(d1)-(strike+shift)*norm.cdf(d2)
        elif type=='Put':
            ExpectedValue=(strike+shift)*norm.cdf(-d2)-(fwdPrice+shift)*norm.cdf(-d1)
        else:
            ExpectedValue=None
            print('tipo opzione non valido')
    return ExpectedValue

# --------- Inversione della formula di Black ------------------------------------
def BlackExpectedValue_tobe_inverted(vol,fwdPrice,strike,TimeToReset,shift,type,price):
    diff = BlackExpectedValue(fwdPrice,strike,TimeToReset,vol,shift,type)-price
    return diff

def BlackExpectedValue_inverse(fwdPrice,strike,TimeToReset,shift,type,price):
    vol=root(fun=BlackExpectedValue_tobe_inverted,args=(fwdPrice,strike,TimeToReset,shift,type,price),x0=0.10)
    return vol.x[0]

# Procedura di Bootstrap delle volatilita' dei Caplets dalle volatilita' flat implicite nei Cap quotati
# tasso di riferimento Euribor 6 mesi

def Bootstrap_CapFloor_ATM(shift, discount_curve, volsdata):

    # ----- import delle volatilita' su cui eseguire il bootstrap -----
    # strike_list  = volsdata['Strikes']
    vol_mat_list = volsdata['Maturities']
    vols_data    = volsdata['Volatilities']

    # ----- import della curva dei fattori di sconto ----------
    t_zc_list = discount_curve['discount times']
    discount_vec = discount_curve['discount factors']

    # ================================================================
    # Calcolo tempi di esercizio singoli Caplets e strike ATM mancanti
    # ================================================================
    time_max      = float(vol_mat_list[-1])+0.5
    times_list    = np.arange(time_max,step=0.5)

    # fattori di sconto per i tempi di esercizio dei Caplets
    zc_discount_list = np.exp(np.interp(times_list, t_zc_list, np.log(discount_vec)))

    # Tenors
    Tenors = []
    for i in range(1, len(times_list)):
        tenor_tmp = times_list[i] - times_list[i - 1]
        Tenors.append(tenor_tmp)

    # Strike ATM
    StrikeATM=[]
    for i in range(1, len(times_list)):
        numerator = zc_discount_list[0] - zc_discount_list[i]
        denominator = 0
        for j in range(0, i):
            den_tmp = Tenors[j] * zc_discount_list[j + 1]
            denominator += den_tmp
        StrikeATM_tmp = numerator / denominator
        StrikeATM.append(StrikeATM_tmp)

    #==========================================================================
    # Interpolazione lineare delle volatilita' flat relative ai singoli Caplets
    # =========================================================================
    vols = np.interp(times_list[1:],vol_mat_list,vols_data)

    # ============================
    # Calcolo variabili intermedie
    # ============================

    # tassi forward
    forward_rates=[]
    for i in range(0,len(Tenors)):
        P1=zc_discount_list[i]
        P2=zc_discount_list[i+1]
        fwd_tmp=(1/Tenors[i])*((P1/P2)-1)
        forward_rates.append(fwd_tmp)


    # prezzi teorici Cap

    Cap_prices = []
    for i in range(1, len(times_list)):
        Cap_tmp = BlackExpectedValue(forward_rates[0],StrikeATM[i-1],times_list[0],vols[i-1],shift,'Call')
        for j in range(1, i):
            BlackValue = BlackExpectedValue(forward_rates[j], StrikeATM[i-1], times_list[j], vols[i-1], shift, 'Call')
            caplet_tmp = zc_discount_list[j + 1] * Tenors[j] * BlackValue
            Cap_tmp += caplet_tmp
        Cap_prices.append(Cap_tmp)

    # ==================
    # Ciclo di Bootstrap
    # ==================

    Caplet_prices=[]
    vol_boot=[]
    Caplet_prices.append(Cap_prices[0])
    vol_boot.append(vols[0])
    for i in range(1,len(Tenors)):
        Caplet_tmp = Cap_prices[i]-Cap_prices[i-1]
        Black_tmp  = Caplet_tmp/(Tenors[i]*zc_discount_list[i+1])

        vol_boot_tmp = BlackExpectedValue_inverse(forward_rates[i - 1], StrikeATM[i], times_list[i], shift, 'Call',Black_tmp)

        Caplet_prices.append(Caplet_tmp)
        vol_boot.append(vol_boot_tmp)

    # moltiplico per 100 le volatilita' e gli strike
    vol_boot  = np.multiply(vol_boot,100)
    StrikeATM = np.multiply(StrikeATM,100)

    bootstrap_result = pd.DataFrame()
    bootstrap_result['Time']= times_list[1:]
    bootstrap_result['Strike (x100)']= StrikeATM
    bootstrap_result['Volatility (x100)']= vol_boot

    return bootstrap_result