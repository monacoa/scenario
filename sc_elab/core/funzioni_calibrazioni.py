from Tkinter import *
import tkMessageBox
import numpy as np
import pandas as pd
from scipy.special import ive
from scipy.stats import norm
from scipy import optimize
import datetime
from anagrafica_dati import *
from sc_elab.excel_hook.W_calibration import W_dividends


def preProcessingCurve(df, rate_time_zero=False, out_type='rate', curve_nom=None):

    # verifico la capitalizzazione dei tassi
    capitalization = df.loc[df.loc[:, 0] == 'Interest rate Type',1].values[0]

    # leggo la data di riferimento
    ref_date = df.loc[df.loc[:, 0] == 'Date Ref',1].values[0]


    # elimino le intestazioni: in base alle date o ai tempi
    tmp = df.loc[:, 0] == 'Date'
    if tmp.any(axis=0):
        row_sep = df.loc[df.loc[:, 0] == 'Date',].index[0]

    tmp = df.loc[:, 0] == 'Times'
    if tmp.any(axis=0):
        row_sep = df.loc[df.loc[:, 0] == 'Times',].index[0]

    out = df.loc[(row_sep+1):df.shape[0],:]

    # seleziono solo le righe da considerare
    out = out.loc[(out.loc[:,2] == 'Y'),0:1]
    out.columns = ['TIME','VALUE']


    if isinstance(out.iloc[0,0],datetime.datetime):
        out.iloc[:, 0] = (out.iloc[:,0] - ref_date)
        out.iloc[:, 0] = out.iloc[:, 0].apply(lambda x: x.days) / 365.2425

    # verifico che non ci sia tempo pari a 0
    if rate_time_zero==False:
        out = out.loc[out.iloc[:, 0] != 0.,:]

    out.reset_index(drop=True, inplace=True)
    out = out.astype('float')

    out_origin = out.copy()

    if capitalization == 'SMP':
        timeTmp = out.iloc[:, 0]
        out.iloc[:, 1] = -1.0 / timeTmp * np.log(1.0 / (1.0 + out.iloc[:, 1] * timeTmp))

    elif capitalization == 'CMP':
        out.iloc[:, 1] = np.log(1.0 + out.iloc[:, 1])

    elif capitalization == 'CNT':
        print 'Nothing to do, right capitalization to calibrate'

    else:
        tkMessageBox.showinfo("Error", "Non e' presente la capitalizzazione")
        capitalization = ""

    out_to_fit = out.copy()
    out_to_fit['VALUE'] = np.divide(out_to_fit['VALUE'],100.) # si assume che i tassi siano passati moltiplicati per 100
    out_to_fit.sort_values(by='TIME',inplace=True)

    if out_type=='discount':
        out_to_fit['VALUE'] = np.exp(-out_to_fit['VALUE']*out_to_fit['TIME'])

    if df.loc[df[0]==u'CurveType',1].values[0]==u'Inflation':
        out_to_fit = createRealCurve(curve_nom,out_to_fit)

    return out_origin ,out_to_fit , capitalization


def fromContinuousToCompost(r):
    rate = np.exp(r) - 1.
    return rate



def preProcessignTimeSeries(df,dt_min,dt_max):

    # elimino le intestazioni
    row_sep = df.loc[df.loc[:, 0] == 'Date',].index[0]
    out = df.loc[(row_sep+1):df.shape[0],:]

    # seleziono solo le righe da considerare
    out = out.loc[(out.loc[:,2] == 'Y'),0:1]
    out.columns = ['DATE','VALUE']

    out = out.loc[(out['DATE'] >= dt_min) & (out['DATE'] <= dt_max),]

    out['TIME'] = (out.iloc[:,0] - out.iloc[0,0])
    out['TIME'] = out.iloc[:, 2].apply(lambda x: x.days) / 365.2425


    out.reset_index(drop=True, inplace=True)
    out['VALUE'].astype('float')
    out['TIME'].astype('float')

    out_to_fit = out.copy()

    return out_to_fit


def preProcessingOptions(W_calib, curve):

    optiontype = W_calib.OptionChosen.loc[(W_calib.OptionChosen[0] == 'OptionType'), 1].values[0]

    if optiontype == 'Swaption':
        # leggo le opzioni
        volsdata_noint = W_calib.OptionChosen.loc[(W_calib.OptionChosen[3] == 'Y'), [0, 1, 2]]
        market_data = pd.DataFrame()
        market_data['expiry'] = volsdata_noint.loc[:, 0].map(MaturityFromStringToYear).values.astype(float)
        market_data['maturity'] = volsdata_noint.loc[:, 1].map(MaturityFromStringToYear).values.astype(float)
        market_data['value'] = np.divide(volsdata_noint.loc[:, 2].values.astype(float), 100.)

        # leggo il tipo di contratto e il tenor
        if W_calib.OptionChosen.loc[W_calib.OptionChosen[0] == 'SwaptionType', 1].values[0] == 'Payer':
            call_flag = 1.
        elif W_calib.OptionChosen.loc[W_calib.OptionChosen[0] == 'SwaptionType', 1].values[0] == 'Receiver':
            call_flag = -1.
        else:
            root = Tk()
            tkMessageBox.showerror(title='Errore', message='Tipo di opzione non valido')
            root.destroy()
            return
        tenr = float(MaturityFromStringToYear[W_calib.OptionChosen.loc[W_calib.OptionChosen[0] == 'Tenor swap', 1].values[0]])

        curve_times = curve['TIME'].values
        curve_values = curve['VALUE'].values
        # calcolo i prezzi delle Swaptions
        market_data.reset_index(drop=True, inplace=True)
        for i in xrange(0, int(market_data.shape[0])):
            t = market_data.at[i, "expiry"]
            T = market_data.at[i, "maturity"]
            srate, annuity = forwardSwap(tenr, curve_times, curve_values, 0, t, T)
            market_data.at[i, "swap"] = srate
            if market_data.at[i, "swap"] > 0:
                market_data.at[i, "market price"] = fromVolaToPrice(t, T, tenr, market_data.at[i, "value"], curve_times,
                                                                    curve_values, call_flag)

        # per calibrare seleziono i dati con swap positivo
        market_data = market_data.loc[market_data['swap'] >= 0]
        market_data.reset_index(drop=True, inplace=True)
        print market_data

        return market_data, tenr, call_flag

    elif optiontype == 'Vol Cap Floor':

        # leggo le opzioni
        volsdata_noint = W_calib.OptionChosen.loc[(W_calib.OptionChosen[3] == 'Y'), [0, 1, 2]]
        market_data = pd.DataFrame()
        market_data['time'] = volsdata_noint.loc[:, 0].values.astype(float)
        market_data['strike'] = np.divide(volsdata_noint.loc[:, 1].values.astype(float), 100.)

        shift = float(W_calib.OptionChosen.loc[W_calib.OptionChosen[0] == 'Shift', 1])
        market_data['vols_data'] = np.divide(volsdata_noint.loc[:, 2].values.astype(float), 100.)

        # Calcolo i prezzi di mercato dei Caplet
        market_data['market price'] = np.array(compute_Black_prices(curve, market_data, 1, shift))

        return market_data

    elif optiontype == 'Caplets':

        # leggo le opzioni
        volsdata_noint = W_calib.OptionChosen.loc[(W_calib.OptionChosen[3] == 'Y'), [0, 1, 2]]
        market_data = pd.DataFrame()
        market_data['time'] = volsdata_noint.loc[:, 0].values.astype(float)
        market_data['strike'] = np.divide(volsdata_noint.loc[:, 1].values.astype(float), 100.)

        market_data['market price'] = volsdata_noint.loc[:, 2].values.astype(float)

        return market_data

    elif optiontype == 'PUT / CALL':

        # leggo le opzioni il prezzo iniziale e strike e maturity da cui estrarre la volatilita' di Black-Scholes
        options_noint = W_calib.OptionChosen.loc[(W_calib.OptionChosen[4] == 'Y'), [0, 1, 2, 3]]
        market_data = pd.DataFrame()
        market_data['maturity'] = options_noint.loc[:, 0].values.astype(float)
        market_data['strike'] = options_noint.loc[:, 1].values.astype(float)
        market_data['market price'] = options_noint.loc[:, 2].values.astype(float)
        market_data['type'] = options_noint.loc[:, 3].values.astype(str)

        S0 = W_calib.OptionChosen.loc[(W_calib.OptionChosen.loc[:, 0] == 'Initial Price'), 1].values.astype(float)[0]

        strike = W_calib.OptionChosen.loc[(W_calib.OptionChosen[0] == 'Strike'), 1].values.astype(float)[0]
        maturity = W_calib.OptionChosen.loc[(W_calib.OptionChosen[0] == 'Mat'), 1].values.astype(float)[0]
        # preprocessing dati
        market_data = market_data.sort_values(by=['maturity', 'type', 'strike'])

        # calcolo dei dividendi impliciti nei prezzi delle opzioni
        dividends_data, dividends = implicit_dividends(S0, market_data, curve)
        # gestisco il caso in cui non ci siano dati disponibili per calcolare i dividendi impliciti
        if len(dividends['VALUE']) == 0:

            root = Tk()
            W_dvd = W_dividends(root)
            root.mainloop()

            if W_dvd.res == 0:
                return
            elif W_dvd.res == 1:
                dividends = pd.DataFrame()
                dividends['TIME'] = [0., 10.]
                dividends['VALUE'] = [float(W_dvd.dvd.get()), float(W_dvd.dvd.get())]
                dividends_data = dividends

        return market_data, S0, strike, maturity, dividends_data, dividends

    else:
        root = Tk()
        tkMessageBox.showerror(message='Tipo opzione non valido')
        root.destroy()
        return

############################################################
#       CIR
############################################################

def loss_zc_model_cir(list_model_params, mkt_prices,power, absrel):

    model_price_tmp = compute_zc_cir_rate(list_model_params, mkt_prices["TIME"])
    diff = np.absolute(model_price_tmp - mkt_prices["VALUE"])
    if absrel == 'rel':
        diff = diff /np.absolute(mkt_prices["VALUE"])

    diff = np.power(diff,power)

    return diff.sum()


def compute_zc_cir_rate(p,t):
    r0 = p[0]
    kappa = p[1]
    theta = p[2]
    sigma = p[3]

    h = (kappa * kappa + 2.0 * sigma * sigma) ** (0.5)

    g0 = 2 * kappa * theta / (sigma * sigma)
    g1 = np.exp(t * h) - 1.0
    g2 = np.exp(t * (h + kappa) / 2.0)

    A0 = (2 * h * g2 / (2.0 * h + (kappa + h) * g1))
    B0 = (2.0 * g1 / (2.0 * h + (kappa + h) * g1))

    rate = -(g0 * np.log(A0) - B0 * r0) / t

    return rate


def mle_cir(list_model_params, ts_data):
    kappa = list_model_params[1]
    theta = list_model_params[2]
    sigma = list_model_params[3]

    n_obs = ts_data.shape[0]

    # calcolo il delta t, eliminando la prima riga
    dt = ts_data['TIME'].diff(periods=1).dropna()
    dt.reset_index(drop=True, inplace=True)

    r = ts_data['VALUE']

    dt= dt.astype('float')
    r = r.astype('float')

    # calcolo r_i come r_old
    r_old = r.drop(r.tail(1).index)

    # calcolo r_i+1 come r_new
    r_new = r.drop(r.head(1).index)
    r_new.reset_index(drop=True, inplace=True)

    c = (2.0 * kappa) / (np.power(sigma,2) * (1.0 - np.exp(-kappa * dt)))
    q = ((2.0 * kappa * theta) / (np.power(sigma,2))) - 1.0

    u_i = c * r_old * (np.exp(-kappa * dt))
    v_i = c * r_new

    z_i = 2.0 * np.sqrt(u_i * v_i)
    b_i = ive(q, z_i)

    # controllo sulla grandezza di b_i
    b_i[b_i < 0.0000001] = 0.0000001

    tmp = (u_i + v_i - 0.5 * q * np.log(v_i / u_i) - np.log(b_i) - z_i)

    lnL = -(n_obs - 1) * np.log(c.mean()) + tmp.sum()

    return lnL


def generate_cir_perc(params, data_to_fit):

    r0 = params[0]
    kappa = params[1]
    theta = params[2]
    sigma = params[3]

    time = data_to_fit['TIME']

    g1 = np.exp(-kappa * time)
    g2 = np.exp(-2.0 * kappa * time)

    G1 = 1.0 - g1
    gamma = (sigma * sigma) / kappa

    mid_tmp = r0 * g1 + theta * G1
    var_tmp = r0 * gamma * (g1 - g2) + theta * gamma / 2.0 * G1 * G1

    s1_tmp = mid_tmp + np.sqrt(var_tmp)
    s2_tmp = mid_tmp - np.sqrt(var_tmp)

    model_value = pd.DataFrame()

    model_value['DATE'] = data_to_fit ['DATE']
    model_value['VALUE'] = data_to_fit ['VALUE']
    model_value['MODEL VALUE MEAN'] = mid_tmp
    model_value['MODEL VALUE DOWN'] = s2_tmp
    model_value['MODEL VALUE UP'] = s1_tmp

    return model_value


############################################################
#       VSCK
############################################################

def loss_zc_model_vsck(list_model_params, mkt_prices,power, absrel):

    model_price_tmp = compute_zc_vsck_rate(list_model_params, mkt_prices["TIME"])
    diff = np.absolute(model_price_tmp - mkt_prices["VALUE"])
    if absrel == 'rel':
        diff = diff / np.absolute(mkt_prices["VALUE"])

    diff = np.power(diff,power)

    return diff.sum()


# output='rate' restituisce il tasso, output='discount' restituisce il fattore di sconto
def compute_zc_vsck_rate(list_model_params, T,output='rate'):
    r0 = list_model_params[0]
    kappa = list_model_params[1]
    theta = list_model_params[2]
    sigma = list_model_params[3]

    B0 = (1.0 / kappa) * (1.0 - np.exp(-kappa * T))
    g0 = (sigma * sigma) / (4.0 * kappa)
    G0 = (theta - sigma * sigma / (2.0 * kappa * kappa))
    A0 = np.exp(G0 * (B0 - T) - g0 * B0 * B0)

    if output=='rate':
        model_rate = -(np.log(A0) - B0 * r0) / T
    elif output=='discount':
        model_rate = A0*np.exp(-B0*r0)

    return model_rate


def mle_vsck(list_model_params, ts_data):

    kappa = list_model_params[1]
    theta = list_model_params[2]
    sigma = list_model_params[3]

    # calcolo il delta t, eliminando la prima riga
    dt = ts_data['TIME'].diff(periods=1).dropna()
    dt.reset_index(drop=True, inplace=True)

    r = ts_data['VALUE']

    dt= dt.astype('float')
    r = r.astype('float')

    # calcolo r_i come r_old
    r_old = r.drop(r.tail(1).index)

    # calcolo r_i+1 come r_new
    r_new = r.drop(r.head(1).index)
    r_new.reset_index(drop=True, inplace=True)

    mu_i = theta + (r_old - theta) * np.exp(-kappa * dt)
    var_i = np.power(sigma, 2) / (2 * kappa) * (1.0 - np.exp(-2 * kappa * dt))
    sigma_i = np.sqrt(var_i)
    norm_i = 1.0 / (sigma_i * np.sqrt(2.0 * np.pi))
    mu_tmp = np.power((r_new - mu_i) / sigma_i , 2)
    sum_mu = mu_tmp.sum()
    sum_norm = np.log(norm_i).sum()

    mll = sum_norm - 0.5 * sum_mu
    mll = -mll

    return mll


def generate_vsck_perc(params, data_to_fit):

    r0 = params[0]
    kappa = params[1]
    theta = params[2]
    sigma = params[3]

    time = data_to_fit['TIME']

    g1 = np.exp(-kappa * time)
    g2 = np.exp(-2.0 * kappa * time)

    G1 = (1.0 - g1)
    G2 = (1.0 - g2)

    gamma = (sigma * sigma) / (2.0 * kappa)

    mid_tmp = r0 * g1 + theta * G1
    var_tmp = gamma * G2

    s1_tmp = mid_tmp + np.sqrt(var_tmp)
    s2_tmp = mid_tmp - np.sqrt(var_tmp)

    model_value = pd.DataFrame()

    model_value['DATE'] = data_to_fit ['DATE']
    model_value['VALUE'] = data_to_fit ['VALUE']
    model_value['MODEL VALUE MEAN'] = mid_tmp
    model_value['MODEL VALUE DOWN'] = s2_tmp
    model_value['MODEL VALUE UP'] = s1_tmp

    return model_value

#======= Calibrazione sulle opzioni Cap Floor ================
# Calcolo del singolo Caplet
# I parametri del Vasicek vanno passati come dizionario {r0, k, sigma, theta}
# Va passata la curva dei tassi come dizionario {TIME, VALUE}
def Caplet_Vasicek(parameters, reset_time, exercise_time, nominal_amount, strike):
    r0    = float(parameters['r0'])
    k     = float(parameters['k'])
    theta = float(parameters['theta'])
    sigma = float(parameters['sigma'])

    T = float(reset_time)
    S = float(exercise_time)
    N = float(nominal_amount)

    if T == 0.:
        P6m = compute_zc_vsck_rate([r0,k,theta,sigma],S,'discount')
        forward_rate = ((1. / P6m)-1) / S
        Caplet = P6m*nominal_amount*S*np.maximum(0., forward_rate-strike)
    else:
        PT = compute_zc_vsck_rate([r0, k, theta, sigma], T, 'discount')
        PS = compute_zc_vsck_rate([r0, k, theta, sigma], S, 'discount')
        Nmod = N * (1 + float(strike) * (S - T))
        BTS = (1/k)*(1-np.exp(-k*(S-T)))
        SIGMAP = np.sqrt((1-np.exp(-2*k*T))/(2*k))*BTS
        h = (1/SIGMAP)*np.log((PS*Nmod)/(PT*N))+(SIGMAP/2)

        Caplet = N*PT*norm.cdf(SIGMAP-h)-Nmod*PS*norm.cdf(-h)

    return Caplet



# ---- Lista prezzi dei Caplet da modello Vasicek -----------------------------
# i parametri del modello vengono passati tramite una lista [r0, k, theta, sigma]
def compute_Vasicek_prices(parameters_list, times, strikes):
    parameters_dict={}
    parameters_dict['r0']=parameters_list[0]
    parameters_dict['k']=parameters_list[1]
    parameters_dict['theta'] = parameters_list[2]
    parameters_dict['sigma'] = parameters_list[3]

    caplet_prices=[]

    for i in range(0,len(times)):
        time_tmp   = times[i]
        strike_tmp = strikes[i]
        caplet_tmp = Caplet_Vasicek(parameters_dict, time_tmp - 0.5, time_tmp, 1., strike_tmp)
        caplet_prices.append(caplet_tmp)

    return caplet_prices




# ---------- Funzione di calcolo norma da ottimizzare ------------------------------------
# power=1 metrica di Manhattan
# power=2 metrica euclidea
# se absrel='rel' viene calcolata la norma delle differenze relative
# mkt_data e' un DataFrame contenente tre series 'time', 'strike' e 'market price'
def loss_caplets_Vasicek(list_model_params, mkt_prices, power, absrel):

    model_price_tmp = np.array(compute_Vasicek_prices(list_model_params, mkt_prices['time'], mkt_prices['strike']))
    diff = np.absolute(model_price_tmp - mkt_prices['market price'])
    if absrel == 'rel':
        diff = diff /np.absolute(mkt_prices['market price'])

    diff = np.power(diff,power)

    return diff.sum()

def computeCHI2(mkt, mdl,type_calib='CURVE'):
    if type_calib=='CURVE':
        tmp = np.power(mkt - mdl,2)/np.absolute(mkt)
    elif type_calib=='CURVE_OPT':
        tmp=[]
        for i in range(0,len(mkt)):
            tmp.append(np.power(mkt[i]-mdl[i],2))
            if mkt[i]!=0: tmp[i]=tmp[i]/mkt[i]
        tmp=np.array(tmp)

    return tmp.sum()


####################################################
#       G2++ su Cap e Floor
####################################################

# ====================================================================
#  Funzioni per il calcolo dei Caplet col modello di Black con shift
# ====================================================================

# valore atteso di Black con shift
# si considera come valuation_time il tempo 0
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
        d1=(np.log((fwdPrice+shift)/(strike+shift))+((vol*vol)/2)*TimeToReset)/(vol*np.sqrt(TimeToReset))
        d2=d1-vol*np.sqrt(TimeToReset)
        if type=='Call':
            ExpectedValue=(fwdPrice+shift)*norm.cdf(d1)-(strike+shift)*norm.cdf(d2)
        elif type=='Put':
            ExpectedValue=(strike+shift)*norm.cdf(-d2)-(fwdPrice+shift)*norm.cdf(-d1)
        else:
            ExpectedValue=None
            print('tipo opzione non valido')
    return ExpectedValue

# ------------------- Prezzo Caplet --------------------------------------------------
# Va passata la curva dei fattori di sconto come dizionario {TIME, VALUE}
def Black_shifted_Caplet(curve, reset_time, exercise_time, nominal_amount, strike, shift, vol):

    P0 = np.exp(np.interp(reset_time, curve['TIME'], np.log(curve['VALUE'])))
    P1 = np.exp(np.interp(exercise_time, curve['TIME'], np.log(curve['VALUE'])))
    tenor = float(exercise_time-reset_time)
    forward_rate = ((P0/P1)-1)/tenor

    Caplet = nominal_amount*P1*tenor*BlackExpectedValue(forward_rate,strike,reset_time,vol,shift,'Call')

    return  Caplet


# ---------- Lista prezzi di mercato dei Caplet---------------------------
# mkt data va dato come DataFrame
def compute_Black_prices(curve, mkt_data, nominal_amount, shift):

    caplet_prices=[]
    for i in range(0,len(mkt_data['time'])):
        time = mkt_data['time'][i]
        strike = mkt_data['strike'][i]
        volatility = mkt_data['vols_data'][i]
        caplet_tmp=Black_shifted_Caplet(curve,time-0.5,time,nominal_amount,strike,shift,volatility)
        caplet_prices.append(caplet_tmp)

    return caplet_prices


# ===========================================================
#     Funzioni di calcolo prezzo dei Caplet nel modello G2++
# ===========================================================
# Calcolo del singolo Caplet
# I parametri del G2++ vanno passati come un dizionario {a, b, sigma, eta, rho}
# Va passata la curva dei fattori di sconto come dataframe {TIME, VALUE}
def Caplet_G2pp(curve,parameters,reset_time,exercise_time,nominal_amount,strike):
    a     = float(parameters['a'])
    b     = float(parameters['b'])
    sigma = float(parameters['sigma'])
    eta   = float(parameters['eta'])
    rho   = float(parameters['rho'])

    T= float(reset_time)
    S= float(exercise_time)
    N= float(nominal_amount)
    Nmod=N*(1+float(strike)*(S-T))

    P1  = np.exp(np.interp(T, curve['TIME'], np.log(curve['VALUE'])))
    P2  = np.exp(np.interp(S, curve['TIME'], np.log(curve['VALUE'])))

    if T == 0.:
        forward_rate=((P1/P2)-1)/(S-T)
        Caplet = P2*nominal_amount*(S-T)*np.maximum(0.0,forward_rate-strike)
    else:
        SIGMA = np.sqrt((np.power(sigma,2)/(2*np.power(a,3)))*np.power(1-np.exp(-a*(S-T)),2)*(1-np.exp(-2*a*(T))) \
                        + (np.power(eta,2)/(2*np.power(b,3)))*np.power(1-np.exp(-b*(S-T)),2)*(1-np.exp(-2*b*(T))) \
                        + 2*rho*((sigma*eta)/(a*b*(a+b)))*(1-np.exp(-a*(S-T)))*(1-np.exp(-b*(S-T)))*(1-np.exp(-(a+b)*(T))))

        x1= (np.log((N*P1)/(Nmod*P2))/SIGMA) - ((1/2)*SIGMA)
        x2 = x1 + SIGMA
        Caplet = -Nmod*P2*norm.cdf(x1)+P1*N*norm.cdf(x2)

    return Caplet



# ---- Lista prezzi dei Caplet da modello G2++ -----------------------------
# i parametri del modello vengono passati tramite una lista [a,b,sigma,eta,rho]
def compute_G2pp_prices(parameters_list, curve, times, strikes):
    parameters_dict={}
    parameters_dict['a']=parameters_list[0]
    parameters_dict['sigma'] = parameters_list[1]
    parameters_dict['b'] = parameters_list[2]
    parameters_dict['eta'] = parameters_list[3]
    parameters_dict['rho'] = parameters_list[4]

    caplet_prices=[]
    # for elem in strikes_and_times:
    #     caplet_tmp=Caplet_G2pp(curve,parameters_dict,0.0,elem[1]-0.5,elem[1],1.0,elem[0])
    #     caplet_prices.append(caplet_tmp)

    for i in range(0,len(times)):
        time_tmp   = times[i]
        strike_tmp = strikes[i]
        caplet_tmp = Caplet_G2pp(curve, parameters_dict, time_tmp - 0.5, time_tmp, 1.0, strike_tmp)
        caplet_prices.append(caplet_tmp)

    return caplet_prices

# ---------- Funzione di calcolo norma da ottimizzare ------------------------------------
# power=1 metrica di Manhattan
# power=2 metrica euclidea
# se absrel='rel' viene calcolata la norma delle differenze relative
# mkt_data e' un DataFrame contenente tre series 'time', 'strike' e 'market price'
def loss_G2pp(list_model_params, curve, mkt_prices, power, absrel):

    model_price_tmp = np.array(compute_G2pp_prices(list_model_params, curve, mkt_prices['time'], mkt_prices['strike']))
    diff = np.absolute(model_price_tmp - mkt_prices['market price'])
    if absrel == 'rel':
        diff = diff /np.absolute(mkt_prices['market price'])

    diff = np.power(diff,power)

    return diff.sum()

#################################################
#   G2++ su Swaptions
#################################################
from utils_g2pp_newton import Pt_MKT,found_opt,price_swaption

def forwardSwap(tenor, rf_times,rf_values, t, t_exp, t_mat):
    T = t_exp + t_mat
    P_0_t = Pt_MKT(t, rf_times,rf_values)
    P_0_exp = Pt_MKT(t_exp, rf_times,rf_values)
    P_0_T = Pt_MKT(T, rf_times,rf_values)

    P_t_exp = P_0_exp / P_0_t
    P_t_T = P_0_T / P_0_t

    n = int(t_mat / tenor)

    ForwardAnnuityPrice = 0.0

    for i in range(1, n + 1):
        ti = t_exp + i * tenor
        P_0_ti = Pt_MKT(ti, rf_times,rf_values)
        ForwardAnnuityPrice = ForwardAnnuityPrice + (P_0_ti / P_0_t) * tenor

    swap = (P_t_exp - P_t_T) / ForwardAnnuityPrice
    return swap, ForwardAnnuityPrice


def CurveFromDictToList(rf_Curve):
    rf_times = rf_Curve.keys()
    rf_times.sort()
    rf_values = []

    for j in xrange(0, len(rf_times)):
        rf_values.append(rf_Curve[rf_times[j]])
    return np.array(rf_times), np.array(rf_values)


def fromVolaToPrice(t_exp, t_mat, tenor, vol, rf_times,rf_values,call_type):
    if (vol <= 0.0) or (vol >= 2.0):
        print "in swaption: it takes 0 < vol < 2; found ", vol

    srate, ForwardAnnuityPrice = forwardSwap(tenor, rf_times,rf_values, 0, t_exp, t_mat)
    call = srate * call_type * (2.0 * norm.cdf(call_type * 0.5 * vol * np.sqrt(t_exp)) - 1.0)

    price = ForwardAnnuityPrice * call

    return price


def fromPriceToVola(t_exp, t_mat, tenor,price, curve_dict,call_type):

    srate, ForwardAnnuityPrice = forwardSwap(tenor, curve_dict, 0, t_exp, t_mat)

    z = 0.5* ( price / (srate * call_type * ForwardAnnuityPrice) + 1)
    p = norm.ppf(z)
    v = 2*p/call_type
    vol = v / np.sqrt(t_exp)

    return vol

#################################################
#              Variance Gamma
#################################################

# --------- Fourier Transform per il prezzo Call nel Variance Gamma -----
# parameters e' un dizionario contenenti le chiavi {sigma, nu, theta}
def PhiCaratt(parameters, S0, r, q, T, u):

    sigma = parameters['sigma']
    nu = parameters['nu']
    theta = parameters['theta']

    omega = (1. / nu) * np.log(1 - theta * nu - 0.5 * sigma * sigma * nu)
    base = 1-1j*theta*nu*u+0.5*sigma*sigma*nu*u*u
    fatt2 = np.power(base,-T/nu)
    fatt1 = np.exp(1j*u*(np.log(S0)+(r-q+omega)*T))

    return fatt1*fatt2

def Psi(parameters, S0, r, q, T, alpha, v):

    arg_phi = v-(alpha+1)*1j
    numeratore = np.exp(-r*T)*PhiCaratt(parameters, S0, r, q, T, arg_phi)
    denominatore = alpha*alpha + alpha -v*v +1j*(2*alpha+1)*v

    return numeratore/denominatore

def Call_Fourier_VG(parameters, S0, r, q, T, alpha, eta, N):

    lambd = 2*np.pi/(N*eta)
    # print 'log strike spacing: %.4f' %lambd
    vect = np.arange(N)
    v_m = eta*vect
    b = np.log(S0)-N*lambd/2
    nodes_1 = np.exp(-1j*b*v_m)
    nodes_2 = Psi(parameters, S0, r, q, T, alpha, v_m)
    w = np.ones(N)*4.
    w[0:N:2] = 2.
    w[0] = 1.
    w[N-1] = 1.
    w = eta*w/3.
    nodes = nodes_1*nodes_2*w
    FFT = np.real(np.fft.fft(nodes))
    k_vect = b + lambd*vect

    prices = np.exp(-alpha*k_vect)*FFT/np.pi
    strikes = np.exp(k_vect)

    return strikes, prices

# --------------- funzioni per la calibrazione ----------------------
def alpha_Madan(list_model_params):

    sigma = list_model_params[0]
    nu = list_model_params[1]
    theta = list_model_params[2]

    alphaUpBound = np.sqrt(np.power(theta, 2) / np.power(sigma, 4) + 2 / (np.power(sigma, 2) * nu)) - theta / np.power(
        sigma, 2) - 1
    flag = 0
    while alphaUpBound > 2.1:
        flag = 1
        alphaUpBound = alphaUpBound / 2
    if alphaUpBound > 1.1:
        alphaUpBound -= 0.1
    alpha = alphaUpBound

    return alpha


def compute_VG_prices(list_model_params, S0, curve, dividends, market_data):

    parameters = {}
    parameters['sigma'] = list_model_params[0]
    parameters['nu'] = list_model_params[1]
    parameters['theta'] = list_model_params[2]

    times = market_data['maturity'].drop_duplicates().tolist()
    rates = np.interp(times,curve['TIME'],curve['VALUE'])
    dividends = np.interp(times, dividends['TIME'],dividends['VALUE'])

    alpha = alpha_Madan(list_model_params)

    VG_prices = []
    for i in range(0,len(times)):
        strikesTmp, pricesTmp = Call_Fourier_VG(parameters, S0, rates[i], dividends[i], times[i], alpha, 0.25, 4096)
        strikeMkt = market_data.loc[market_data['maturity']==times[i],['strike']].values.flatten()
        typeMkt = market_data.loc[market_data['maturity']==times[i],['type']].values.flatten()
        price_interp = np.interp(strikeMkt,strikesTmp,pricesTmp).flatten().tolist()
        for j in range(0,len(strikeMkt)):
            if typeMkt[j]=='PUT':
                price_interp[j] = price_interp[j]-S0*np.exp(-dividends[i]*times[i])+strikeMkt[j]*np.exp(-rates[i]*times[i])
        VG_prices = VG_prices + price_interp

    return VG_prices

# ---------- Funzione di calcolo norma da ottimizzare ------------------------------------
# power=1 metrica di Manhattan
# power=2 metrica euclidea
# se absrel='rel' viene calcolata la norma delle differenze relative
# mkt_data e' un DataFrame contenente tre series 'maturity', 'strike' e 'market price'
# curve e' un DataFrame contenente due series: 'time' e 'rate'
def loss_Call_VG(list_model_params, S0, mkt_data, curve, dividends, power, absrel):

    model_price_tmp = np.array(compute_VG_prices(list_model_params, S0, curve, dividends, mkt_data[['maturity','strike','type']]))
    diff = np.absolute(model_price_tmp - mkt_data['market price'])
    if absrel == 'rel':
        diff = diff /np.absolute(mkt_data['market price'])

    diff = np.power(diff,power)

    return diff.sum()

# ---- Prezzo analitico Call nel Variance Gamma --------
# non viene utilizzato nella calibrazione perche' presenta dei problemi numerici quando il
# rapporto T/nu e' molto grande. Tuttavia la formula funziona per la maggiorparte delle maturita' a breve.
def funcPsi(a,b,gam):

    c = np.absolute(a)*np.sqrt(2.+np.power(b,2))
    u = b/np.sqrt(2.+np.power(b,2))
    v1 = np.sign(a)*c
    v2 = np.power(c,gam+0.5)*np.exp(v1)/(np.sqrt(2*np.pi)*special.gamma(gam))

    integrand1 = lambda s: np.power(s, gam - 1.) * np.power(1. - 0.5 * s * (1. + u), gam - 1.) * np.exp(
        -s * v1 * (1. + u))
    x_values = np.arange(0,1.01,0.01)
    integral1 = integrate.simps(integrand1(x_values),x_values)
    Phi1 = gam*integral1
    integrand2 = lambda s: np.power(s, gam) * np.power(1. - 0.5 * s * (1. + u), gam - 1.) * np.exp(-s * v1 * (1. + u))
    integral2 = integrate.simps(integrand2(x_values),x_values)
    Phi2 = (1.+gam)*integral2

    Psi = v2*np.power(1.+u,gam)/gam*special.kv(gam+0.5,c)*Phi1 - \
          np.sign(a)*v2*np.power(1.+u,1.+gam)/(1.+gam)*special.kv(gam-0.5,c)*Phi2 +\
          np.sign(a)*v2*np.power(1.+u,gam)/gam*special.kv(gam-0.5,c)*Phi1

    return Psi

def Call_VG(S0,K,T,r,sigma,nu,theta):

    gamm = T/nu
    omega = (1/nu)*np.log(1-nu*theta-0.5*sigma*sigma*nu)
    zeta = (np.log(S0/K)+omega*T)/sigma
    vartheta = 1-nu*(theta+0.5*sigma*sigma)
    a_1 = zeta*np.sqrt(vartheta/nu)
    b_1 = (1./sigma)*(theta + sigma*sigma)*np.sqrt(nu/vartheta)
    a_2 = zeta*np.sqrt(1./nu)
    b_2 = (1./sigma)*theta*np.sqrt(nu)

    fattPsi_1 = funcPsi(a_1,b_1,gamm)
    fattPsi_2 = funcPsi(a_2,b_2,gamm)
    Call_price = S0*np.exp(-r*T)*fattPsi_1 \
                 -K*np.exp(-r*T)*fattPsi_2

    return Call_price

# -------- funzione di Calcolo dei dividendi implciti nei prezzi delle opzioni ---------
def implicit_dividends(S0,market_data,curve):

    # eseguo un reshape dei dati per aggregare i prezzi di Call e Put corrispondenti agli stessi strike e maturita'
    base_c = market_data.loc[market_data['type'] == "CALL",['maturity','strike','market price']].rename(columns={'market price':'price CALL'})
    base_p = market_data.loc[market_data['type'] == "PUT",['maturity','strike','market price']].rename(columns={'market price':'price PUT'})
    reshaped_data = pd.merge(base_c,base_p, how='inner', on=['maturity','strike'])

    # controlli sui dati disponibili
    maturity_list = market_data['maturity'].drop_duplicates().tolist()
    reshaped_maturity_list = reshaped_data['maturity'].drop_duplicates().tolist()
    # controllo di avere dati su cui eseguire l'estrazione dei dividendi
    if len(reshaped_maturity_list)==0:
        print 'No data available to compute implcit dividends'
        dividend = pd.DataFrame()
        dividend['TIME']  = []
        dividend['VALUE'] = []
        return reshaped_data , dividend
    # controllo per quali maturita' non ho dati su cui calcolare il dividendo implicito
    for t in maturity_list:
        if t not in reshaped_maturity_list:
            print 'No Call and Put options available for the maturity %.4f year \n' \
                  'No implicit dividend available at that maturity.' %t

    # aggiungo al df reshaped una colonna con il tasso risk free
    rates = np.interp(reshaped_data['maturity'],curve['TIME'],curve['VALUE'])
    reshaped_data['discount factor'] = np.exp(-rates*reshaped_data['maturity'])

    # calcolo i dividendi impliciti
    log_argument = np.divide(reshaped_data['price CALL']-reshaped_data['price PUT']+reshaped_data['strike']*reshaped_data['discount factor'],S0)
    reshaped_data['dividend rate'] = -np.divide(np.log(log_argument),reshaped_data['maturity'])

    # calcolo la media dei dividendi impliciti per maturita' e restituisco la curva dei dividendi impliciti
    res = pd.DataFrame()
    res['TIME'] = reshaped_maturity_list
    res['VALUE'] = reshaped_data.groupby(['maturity']).mean()['dividend rate'].tolist()

    return reshaped_data, res

# -------- Call Black-Scholes e inversione -------------------------------------------
def Call_BS(S0,K,T,r,q,sigma):
    d1 = (np.log(S0/K)+(r-q+0.5*np.power(sigma,2))*T)/(sigma*np.sqrt(T))
    d2 = d1 - sigma*np.sqrt(T)
    Call = np.exp(-q*T)*S0*norm.cdf(d1)-K*np.exp(-r*T)*norm.cdf(d2)
    return Call

def fromPriceBSToVol(Price, S0, K, T, r, q):
    diff = lambda s: Call_BS(S0,K,T,r,q,s)-Price
    vol = optimize.root(diff,x0=0.10)
    return vol.x[0]

# Trovo la volatilita' Black-Scholes associata a strike e maturita' date in input tale che il prezzo
# Black-Scholes e quello Variance Gamma coincidano
def fromPriceVGtoVolBS(parameters_list,S0,strike,maturity,curve,dividends):

    parameters = {'sigma': parameters_list[0], 'nu': parameters_list[1], 'theta': parameters_list[2]}
    r = np.interp(maturity, curve['TIME'], curve['VALUE'])
    q = np.interp(maturity, dividends['TIME'], dividends['VALUE'])
    alpha = alpha_Madan(parameters_list)
    strike_list, price_list = Call_Fourier_VG(parameters, S0, r, q, maturity, alpha, 0.25, np.power(2, 12))
    priceFourier = np.interp(strike, strike_list, price_list)
    vol = fromPriceBSToVol(priceFourier, S0, strike, maturity, r, q)
    priceBS = Call_BS(S0, strike, maturity, r, q, vol)
    if np.abs(priceFourier-priceBS) > 2.*priceFourier/100. :
        root = Tk()
        tkMessageBox.showwarning(title='Attenzione',message='Inversione del prezzo non riuscita, BS volatility non attendibile')
        root.mainloop()

    return vol

###############################################
#  Jarrow Yildirim
###############################################

# Preprocessing data
def preProcessingDataJY(W_calib,curve_nom,real_curve):

    flag_error = False
    data = W_calib.OptionChosen

    param = {}
    param['I0']=float(W_calib.setting_I0.get())
    param['Strike_scale']=float(W_calib.setting_strike_scale.get())
    param['K']=float(W_calib.setting_K.get())
    param['Noz']=float(W_calib.setting_Noz.get())
    param['curve_n']=curve_nom
    param['curve_r']=real_curve
    param['tenorOplet']=data.loc[data[0]=='tenorOplet',1].values.astype(float)[0]

    if data.loc[data[0]=='OptionType',1].values[0]=='Floor':
        param['w']='floor'
    elif data.loc[data[0]=='OptionType',1].values[0]=='Cap':
        param['w']='cap'
    else:
        flag_error = True

    if data.loc[data[0]=='Type value',1].values[0]=='Price':
        flag_optim = True # Se True la calibrazione avviene sui prezzi di mercato delle opzioni
    elif data.loc[data[0]=='Type value',1].values[0]=='Volatility':
        flag_optim = False # Se False la calibrazione avviene sulla volatilita' delle opzioni
    else:
        flag_error = True

    market_data = data.loc[data[2]=='Y',[0,1]].sort_values(by=0).astype(float).reset_index(drop=True)
    market_data.rename(columns={0:'time',1:'market value'},inplace=True)

    if flag_error == True:
        root=Tk()
        tkMessageBox.showwarning(title='Error',message='Dati opzioni non validi, ricontrollare.')
        root.destroy()
        return

    return param, flag_optim, market_data

# -----------------------------------------------------
# Funzioni per curva reale
# -----------------------------------------------------

# fattori di sconto CONTINUI
def computeDfCurve(rf_Curve, time):
    rf_times = rf_Curve['TIME'].tolist()
    rf_values = rf_Curve['VALUE'].tolist()

    RateTime = np.interp(time, rf_times, rf_values)
    dfCurveTime = np.exp(- RateTime * time)

    return dfCurveTime

# Curva reale al tasso di interesse CONTINUO
def createRealCurve(curve_n, curve_i):

    times_n = curve_n['TIME'].tolist()
    times_i = curve_i['TIME'].tolist()

    if len(times_n)>0:
        times_r = times_n
    elif len(times_i)>0:
        times_r = times_i
    else:
        root = Tk()
        tkMessageBox.showwarning(title = 'Curve error', message = 'Lunghezza curve in input non valida')
        root.destroy()
        return

    rates_n = np.interp(times_r,curve_n['TIME'],curve_n['VALUE'])
    rates_i = np.interp(times_r,curve_i['TIME'],curve_i['VALUE'])

    curve_r = pd.DataFrame()
    curve_r['TIME']  = times_r
    curve_r['VALUE'] = rates_n - rates_i

    return curve_r

# -----------------------------------------------------
# no fixing Index
# -----------------------------------------------------

def A(ts, te, alpha):
    return (1. - np.exp(-alpha * (te - ts))) / alpha


# ----
def C(t0, ts, te, aN, aR, sigmaR, sigmaN, sigmaI, rhoNR, rhoNI, rhoRI):
    cc = sigmaR * A(ts, te, aR)
    tmp = rhoRI * sigmaI - .5 * sigmaR * A(t0, ts, aR)
    tmp += (rhoNR * sigmaN) / (aN + aR) * (1. + aR * A(t0, ts, aN))
    tmp *= A(t0, ts, aR)
    tmp -= A(t0, ts, aN) * (rhoNR * sigmaN) / (aN + aR)
    cc *= tmp

    return cc


# ----

def VarYoY(t0, ts, te, aN, aR, sigmaR, sigmaN, sigmaI, rhoNR, rhoNI, rhoRI):
    Dt1 = ts - t0
    Dt2 = te - ts
    exp1N = np.exp(- aN * Dt2)
    exp2N = np.exp(- 2. * aN * Dt2)
    num1N = 1. - exp1N

    exp1R = np.exp(- aR * Dt2)
    exp2R = np.exp(- 2. * aR * Dt2)
    num1R = 1. - exp1R

    vv = sigmaN * sigmaN / (2. * pow(aN, 3.)) * num1N * num1N * (1. - np.exp(-2. * aN * Dt1)) + sigmaI * sigmaI * Dt2

    vv += sigmaR * sigmaR / (2. * pow(aR, 3.)) * num1R * num1R * (1. - np.exp(-2. * aR * Dt1))
    vv += - 2. * rhoNR * (sigmaN * sigmaR) / (aR * aN * (aN + aR)) * num1N * num1R * (1 - np.exp(-(aN + aR) * Dt1))

    vv += sigmaN * sigmaN / (aN * aN) * (Dt2 + 2. / aN * exp1N - .5 / aN * exp2N - 3. / (2. * aN))

    vv += sigmaR * sigmaR / (aR * aR) * (Dt2 + 2. / aR * exp1R - .5 / aR * exp2R - 3. / (2. * aR))
    vv += - 2. * rhoNR * (sigmaN * sigmaR) / (aR * aN) * (
                Dt2 - num1N / aN - num1R / aR + (1. - np.exp(-(aN + aR) * Dt2)) / (aN + aR))
    vv += 2. * rhoNI * (sigmaN * sigmaI) / aN * (Dt2 - num1N / aN) - 2. * rhoRI * (sigmaR * sigmaI) / aR * (
                Dt2 - num1R / aR)
    return vv


# ----

def VarInfIndex(t, aN, aR, sigmaR, sigmaN, sigmaI, rhoNR, rhoNI, rhoRI):
    dN = (np.exp(- aN * t) - 1.) / aN
    dR = (np.exp(- aR * t) - 1.) / aR

    eta2N = (sigmaN * sigmaN) / (aN * aN) * (t + 2. * dN + (1 - np.exp(- 2. * aN * t)) / (2. * aN))
    eta2R = (sigmaR * sigmaR) / (aR * aR) * (t + 2. * dR + (1 - np.exp(- 2. * aR * t)) / (2. * aR))

    eta2NR = (sigmaN * sigmaR * rhoNR * 2.) / (aN * aR)
    eta2NR = eta2NR * (t + dN + dR + (1. - np.exp(-(aN + aR) * t)) / (aN + aR))

    eta2NI = (sigmaN * sigmaI * rhoNI * 2.) / aN * (t + dN)
    eta2RI = (sigmaR * sigmaI * rhoRI * 2.) / aR * (t + dR)
    eta2I = (sigmaI * sigmaI) * t
    eta2 = eta2N + eta2R - eta2NR + eta2NI - eta2RI + eta2I
    return eta2


# ----

# -----------------------------------------------------
# no fixing Index
# -----------------------------------------------------

def analitical_oplet_jy(vals, curve_n, curve_r, Verbose=False):
    aN = vals['aN']
    aR = vals['aR']
    sigmaR = vals['sigmaR']
    sigmaN = vals['sigmaN']
    sigmaI = vals['sigmaI']
    rhoNR = vals['rhoNR']
    rhoNI = vals['rhoNI']
    rhoRI = vals['rhoRI']

    Pn_s = computeDfCurve(curve_n, vals['ts'])
    Pn_e = computeDfCurve(curve_n, vals['te'])
    Pr_s = computeDfCurve(curve_r, vals['ts'])
    Pr_e = computeDfCurve(curve_r, vals['te'])

    # Calcolo media
    cc = C(vals['t0'], vals['ts'], vals['te'], aN, aR, sigmaR, sigmaN, sigmaI, rhoNR, rhoNI, rhoRI)
    m = (Pn_s * Pr_e) / (Pn_e * Pr_s) * np.exp(cc)
    # print m

    # Calcolo varianza
    vv = VarYoY(vals['t0'], vals['ts'], vals['te'], aN, aR, sigmaR, sigmaN, sigmaI, rhoNR, rhoNI, rhoRI)

    if vv < 1.e-8:
        print "@" * 20
        return max(vals['omega'] * (m - vals['Strike']), 0.)

    # Calcolo oplet
    d1 = vals['omega'] * (np.log(m / vals['Strike']) + .5 * vv) / np.sqrt(vv)
    d2 = vals['omega'] * (np.log(m / vals['Strike']) - .5 * vv) / np.sqrt(vv)
    expValue = vals['omega'] * (m * norm.cdf(d1) - vals['Strike'] * norm.cdf(d2))

    if expValue < 0.:
        print "@" * 20
        if expValue < -1.e-8:
            raise Exception("%s bruttissimo valore per un'opzione" % expValue)
        else:
            expValue = 0.

    if Verbose:
        print "%12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f" % (
        cc, vals['ts'], vals['te'], m, vv, d1, d2, expValue, Pn_e, (Pn_e * expValue))

    return expValue


# -----------------------------------------------------
# fixing Index a ts
# -----------------------------------------------------
def analytical_option_jy(vals, curve_n, curve_r):
    I0 = vals['I0']
    K = vals['Strike']
    t = vals['te']

    Pn = computeDfCurve(curve_n, t)
    Pr = computeDfCurve(curve_r, t)

    # Calcolo il forward dell'indice d'inflazione
    fwdI = I0 * Pr / Pn

    # Calcolo la volatility flat
    eta2 = VarInfIndex(t, vals['aN'], vals['aR'], vals['sigmaR'], vals['sigmaN'], vals['sigmaI'], vals['rhoNR'],
                       vals['rhoNI'], vals['rhoRI'])

    if eta2 < 1.e-8:
        print "@" * 20
        return max(vals['omega'] * (fwdI - K), 0.)

    d1 = vals['omega'] * (np.log(fwdI / K) + eta2 * .5) / np.sqrt(eta2)
    d2 = vals['omega'] * (np.log(fwdI / K) - eta2 * .5) / np.sqrt(eta2)

    expValue = vals['omega'] * (I0 * Pr * norm.cdf(d1) - K * Pn * norm.cdf(d2)) / Pn
    # print Pr, I0, Pn, eta2
    if expValue < 0.:
        print "@" * 20
        if expValue < -1.e-8:
            raise Exception("%s bruttissimo valore per un'opzione" % expValue)
        else:
            expValue = 0.
    return expValue


def inflationOptionJY(list_model_params, param, schedule):
    curve_n = param['curve_n']
    curve_r = param['curve_r']

    noz = param['Noz']
    I0 = param['I0']

    if param['w'] == 'cap':
        w = 1.0
    elif param['w'] == 'floor':
        w = -1.0
    else:
        print "skipping option type %s" % param['w']

    # scale = param['Strike_scale']
    strike = param['K'] / param['Strike_scale']

    vals = {"t0": 0.
        , "I0": I0
        , "aN": list_model_params[0]
        , "aR": list_model_params[1]
        , "sigmaR": list_model_params[2]
        , "sigmaN": list_model_params[3]
        , "sigmaI": list_model_params[4]
        , "rhoNR": list_model_params[5]
        , "rhoNI": list_model_params[6]
        , "rhoRI": list_model_params[7]
        , "omega": w
            }

    first = True
    price = 0.0

    if vals['sigmaR'] < 1.e-8 or vals['sigmaN'] < 1.e-8 or vals['sigmaI'] < 1.e-8:
        print '*********** bruttissimo ***************'
        return 0.

    for i in range(len(schedule)):
        if i == 0:
            continue
        if first:
            # richiamo analytical_call_jy(I0, K, t, Pn, Pr, param)
            t = schedule[i]
            vals['te'] = t
            Dt = t - schedule[i - 1]
            vals['Strike'] = I0 * (1.0 + strike)

            price = noz * computeDfCurve(curve_n, t) * Dt * analytical_option_jy(vals, curve_n, curve_r) / I0
            first = False
        # return price
        else:
            # richiamo analitical_caplet_jy
            ts = schedule[i - 1]
            te = schedule[i]
            Dt = te - ts
            vals['ts'] = ts
            vals['te'] = te
            vals['Strike'] = (1.0 + strike)
            jy = noz * computeDfCurve(curve_n, te) * Dt * analitical_oplet_jy(vals, curve_n, curve_r)
            price += jy
    return price

#--------- funzioni loss -------------------------------------------
def loss_jy_model(list_model_params , param , mkt_prices):
    diff_sum = 0.0
    time_mkt_list = mkt_prices['time'].tolist()

    for t in range(0,len(time_mkt_list)):
        schedule = np.arange(0.0, time_mkt_list[t],  param["tenorOplet"]/12.).tolist()
        schedule.append(time_mkt_list[t])
        model_price_tmp = inflationOptionJY(list_model_params , param , schedule)
        mkt_price_tmp = mkt_prices['market value'][t]
        diff = abs(model_price_tmp - mkt_price_tmp)
        diff = diff * diff
        diff_sum = diff_sum + diff

    return diff_sum

def loss_jy_model_var(list_model_params , wishVol):

    aN	    = list_model_params[0]
    aR	    = list_model_params[1]
    sigmaR  = list_model_params[2]
    sigmaN  = list_model_params[3]
    sigmaI  = list_model_params[4]
    rhoNR   = list_model_params[5]
    rhoNI   = list_model_params[6]
    rhoRI   = list_model_params[7]

    diff_sum = 0.0
    time_mkt_list = wishVol['time'].tolist()

    for t in range(0,len(time_mkt_list)):
        model_vol_tmp = np.sqrt(VarInfIndex(time_mkt_list[t], aN, aR, sigmaR, sigmaN, sigmaI, rhoNR, rhoNI, rhoRI))
        mkt_vol_tmp = wishVol['market value'][t]
        diff = abs(model_vol_tmp - mkt_vol_tmp)
        diff = diff * diff
        diff_sum = diff_sum + diff

    return diff_sum

# ---------- funzione post calibrazione ---------------------
def compute_values_post_calib_JY(flag_optim, param, market_data, calib_result_list):

    mdl_opt_list = []

    if flag_optim:
        t_ref_list = market_data['time'].tolist()

        for tt in range(0,len(t_ref_list)):
            schedule = np.arange(0.0, t_ref_list[tt], param["tenorOplet"] / 12.).tolist()
            schedule.append(t_ref_list[tt])
            mdl_tmp = inflationOptionJY(calib_result_list, param, schedule)

            mdl_opt_list.append(float(mdl_tmp))

    else:
        t_ref_list = market_data['time'].tolist()

        for tt in t_ref_list:
            mdl_tmp = VarInfIndex(float(tt), calib_result_list[0], calib_result_list[1], calib_result_list[2],
                                  calib_result_list[3], calib_result_list[4], calib_result_list[5],
                                  calib_result_list[6], calib_result_list[7])

            mdl_opt_list.append(float(mdl_tmp))

    return mdl_opt_list