import numpy as np
import pandas as pd
from scipy.optimize import root
from scipy.stats import norm
import matplotlib.pyplot as plt

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

# --------- Funzione di calcolo di diversi Caplet sullo stesso strike ---------------
def Multiple_Caplets_Price(zc_df_list, fwd_list, strike, times_list, sigma, shift):
    if len(fwd_list) != len(times_list): raise Exception('lunghezza lista tempi e lunghezza lista forward diverse!')
    price_sum = 0.
    for i in range(0,len(fwd_list)):
        BlackValueTmp = BlackExpectedValue(fwd_list[i], strike, times_list[i], sigma, shift,
                                        'Call')
        caplet_tmp = zc_df_list[i] * 0.5 * BlackValueTmp
        price_sum += caplet_tmp
    return price_sum


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

    # ============================
    # Calcolo variabili intermedie
    # ============================

    # tassi forward
    forward_rates=[]
    for i in range(0,len(Tenors)):
        P1=zc_discount_list[i]
        P2=zc_discount_list[i+1]
        fwd_tmp=(1./Tenors[i])*((P1/P2)-1.)
        forward_rates.append(fwd_tmp)


    # prezzi teorici Cap

    Cap_prices = []
    idx = 1
    for i in range(0, len(vol_mat_list)):
        while times_list[idx] != vol_mat_list[i]:
            idx+=1
        Cap_tmp = 0
        for j in range(0, idx):
            BlackValue = BlackExpectedValue(forward_rates[j], StrikeATM[idx-1], times_list[j], vols_data[i], shift, 'Call')
            caplet_tmp = zc_discount_list[j + 1] * Tenors[j] * BlackValue
            Cap_tmp += caplet_tmp
        Cap_prices.append(Cap_tmp)

    # ==================
    # Ciclo di Bootstrap
    # ==================

    Caplet_prices=[]
    vol_boot=[]
    idx = 1
    vol_boot.append(vols_data[0])
    while vol_mat_list[0]!= times_list[idx]:
        vol_boot.append(vols_data[0])
        idx += 1

    idx_n = idx
    for i in range(1,len(Cap_prices)):
        diff_tmp = Cap_prices[i]-Cap_prices[i-1]
        while times_list[idx_n] != vol_mat_list[i]:
            idx_n += 1
        caplets_prices_tmp = lambda sigma : Multiple_Caplets_Price(zc_discount_list[idx+1:idx_n+1], forward_rates[idx:idx_n], StrikeATM[idx_n-1], times_list[idx:idx_n], sigma, shift) - diff_tmp
        vol_boot_tmp = root(fun=caplets_prices_tmp,x0=0.10).x[0]
        for j in range(0,idx_n-idx):
            vol_boot.append(vol_boot_tmp)
        idx = idx_n

    # moltiplico per 100 le volatilita' e i forward
    vol_boot  = np.multiply(vol_boot,100)
    forward_rates = np.multiply(forward_rates,100)

    bootstrap_result = pd.DataFrame()
    bootstrap_result['Time']= times_list[1:]
    bootstrap_result['Fwd rate (x100)']= forward_rates
    bootstrap_result['Volatility (x100)']= vol_boot

    return bootstrap_result


def Bootstrap_CapFloor_Surface(shift, discount, volatilities):

    df_input_data=pd.DataFrame()
    df_input_data['Maturity'] = volatilities['Maturities']
    df_input_data['Strike']   = volatilities['Strikes']
    df_input_data['Volatility'] = volatilities['Volatilities']

    vol_mat_list = df_input_data['Maturity'].drop_duplicates().values.astype(float)
    strike_list = df_input_data['Strike'].drop_duplicates().values.astype(float)
    vols_data = df_input_data.pivot_table(index='Maturity', columns='Strike', values='Volatility',fill_value = None).values.astype(float)
    vols_data = vols_data.transpose()

    # --------- Calcolo tempi di esercizio singoli Caplets -----
    time_max = float(vol_mat_list[-1]) + 0.5
    times_list = np.arange(time_max, step=0.5)

    # ----------------------------------------------------------

    # ------- Interpolazione delle volatilita' flat relative ai singoli Caplets ---------
    # uso una semplice interpolazione lineare sulle maturita' per ogni strike (riga per riga)
    number_strike = len(strike_list)
    vols = []
    for i in range(0, number_strike):
        vols_tmp = np.interp(times_list[1:], vol_mat_list, vols_data[i])
        vols.append(vols_tmp)

    # ----------------------------------------------------------------------------------

    # ----- Calcolo variabili intermedie ------------------------
    # Tenors
    Tenors = []
    for i in range(1, len(times_list)):
        tenor_tmp = times_list[i] - times_list[i - 1]
        Tenors.append(tenor_tmp)

    # fattori di sconto per i tempi di esercizio dei Caplets
    zc_discount_list = np.exp(np.interp(times_list, discount['discount times'], np.log(discount['discount factors'])))

    # tassi forward
    forward_rates = []
    for i in range(0, len(Tenors)):
        P1 = zc_discount_list[i]
        P2 = zc_discount_list[i + 1]
        fwd_tmp = (1. / Tenors[i]) * ((P1 / P2) - 1.)
        forward_rates.append(fwd_tmp)

    # prezzi teorici Cap
    Cap_prices = []
    for s in range(0, number_strike):
        Cap_prices.append([])
        for i in range(2, len(times_list)):
            Cap_tmp = 0
            for j in range(0, i):
                BlackValue = BlackExpectedValue(forward_rates[j], strike_list[s], times_list[j], vols[s][i - 1], shift,
                                                'Call')
                caplet_tmp = zc_discount_list[j + 1] * Tenors[j] * BlackValue
                Cap_tmp += caplet_tmp
            Cap_prices[s].append(Cap_tmp)

    # --------------------- Ciclo di Bootstrap ---------------------
    Caplet_prices = []
    vol_boot = []
    for s in range(0, number_strike):
        Caplet_prices.append([])
        vol_boot.append([])
        vol_boot[s].append(vols[s][1])
        for i in range(1, len(Cap_prices[s])):
            Caplet_tmp = Cap_prices[s][i] - Cap_prices[s][i - 1]
            Black_tmp = Caplet_tmp / (Tenors[i+1] * zc_discount_list[i + 2])

            vol_boot_tmp = BlackExpectedValue_inverse(forward_rates[i+1], strike_list[s], times_list[i+1], shift,
                                                      'Call', Black_tmp)

            Caplet_prices[s].append(Caplet_tmp)
            vol_boot[s].append(vol_boot_tmp)
    vol_boot = np.array(vol_boot)

    for s in range(len(vol_boot)):
        plt.plot(times_list[2:], vol_boot[s])
    plt.show()

    # moltiplico per 100 le volatilita' e gli strike e ridimensiono gli output per la scrittura
    vol_boot = np.multiply(vol_boot, 100.).flatten()
    strike_list = np.repeat(np.multiply(strike_list, 100.),len(times_list[2:]))

    bootstrap_result = pd.DataFrame()
    bootstrap_result['Time'] = np.tile(times_list[2:],number_strike)
    bootstrap_result['Strike (x100)'] = strike_list
    bootstrap_result['Volatility (x100)'] = vol_boot

    return bootstrap_result