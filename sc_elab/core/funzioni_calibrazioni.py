import tkMessageBox
import numpy as np
import pandas as pd
from scipy.special import ive
from scipy.stats import norm
import datetime


def preProcessignCurve(df):

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
# Va passata la curva dei tassi come dizionario {t_zc_list, rf_rate_list}
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
#       G2++
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
# Va passata la curva dei fattori di sconto come dizionario {t_zc_list, discount_vec}
def Black_shifted_Caplet(curve, reset_time, exercise_time, nominal_amount, strike, shift, vol):

    P0 = np.exp(np.interp(reset_time, curve['t_zc_list'], np.log(curve['discount_vec'])))
    P1 = np.exp(np.interp(exercise_time, curve['t_zc_list'], np.log(curve['discount_vec'])))
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
# Va passata la curva dei fattori di sconto come dizionario {t_zc_list, discount_vec}
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

    P1  = np.exp(np.interp(T, curve['t_zc_list'], np.log(curve['discount_vec'])))
    P2  = np.exp(np.interp(S, curve['t_zc_list'], np.log(curve['discount_vec'])))

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