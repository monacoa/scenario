import tkMessageBox
import numpy as np
import pandas as pd
from scipy.special import ive
from scipy.stats import norm
import datetime


def preProcessignCurve(df,rate_time_zero=False, out_type='rate'):

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

    if out_type=='discount':
        out_to_fit['VALUE'] = np.exp(-out_to_fit['VALUE']*out_to_fit['TIME'])

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
#              Variance Gamma
#################################################

# --------- Fourier Transform per il prezzo Call nel Variance Gamma -----
# parameters e' un dizionario contenenti le chiavi {sigma, nu, theta}
def PhiCaratt(parameters, S0, r, T, u):

    sigma = parameters['sigma']
    nu = parameters['nu']
    theta = parameters['theta']

    omega = (1. / nu) * np.log(1 - theta * nu - 0.5 * sigma * sigma * nu)
    base = 1-1j*theta*nu*u+0.5*sigma*sigma*nu*u*u
    fatt2 = np.power(base,-T/nu)
    fatt1 = np.exp(1j*u*(np.log(S0)+(r+omega)*T))

    return fatt1*fatt2

def Psi(parameters, S0, r, T, alpha, v):

    arg_phi = v-(alpha+1)*1j
    numeratore = np.exp(-r*T)*PhiCaratt(parameters, S0, r, T, arg_phi)
    denominatore = alpha*alpha + alpha -v*v +1j*(2*alpha+1)*v

    return numeratore/denominatore

def Call_Fourier_VG(parameters, S0, r, T, alpha, eta, N):

    lambd = 2*np.pi/(N*eta)
    # print 'log strike spacing: %.4f' %lambd
    vect = np.arange(N)
    v_m = eta*vect
    b = np.log(S0)-N*lambd/2
    nodes_1 = np.exp(-1j*b*v_m)
    nodes_2 = Psi(parameters, S0, r, T, alpha, v_m)
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


def compute_VG_prices(list_model_params, S0, curve, market_data):

    parameters = {}
    parameters['sigma'] = list_model_params[0]
    parameters['nu'] = list_model_params[1]
    parameters['theta'] = list_model_params[2]

    times = market_data['maturity'].drop_duplicates().tolist()
    rates = np.interp(times,curve['TIME'],curve['VALUE'])

    alpha = alpha_Madan(list_model_params)

    VG_prices = []
    for i in range(0,len(times)):
        strikesTmp, pricesTmp = Call_Fourier_VG(parameters, S0, rates[i], times[i], alpha, 0.25, 4096)
        strikeMkt = market_data.loc[market_data['maturity']==times[i],['strike']].values.flatten()
        typeMkt = market_data.loc[market_data['maturity']==times[i],['type']].values.flatten()
        price_interp = np.interp(strikeMkt,strikesTmp,pricesTmp).flatten().tolist()
        for j in range(0,len(strikeMkt)):
            if typeMkt[j]=='PUT':
                price_interp[j] = price_interp[j]-S0+strikeMkt[j]*np.exp(-rates[i]*times[i])
        VG_prices = VG_prices + price_interp

    return VG_prices

# ---------- Funzione di calcolo norma da ottimizzare ------------------------------------
# power=1 metrica di Manhattan
# power=2 metrica euclidea
# se absrel='rel' viene calcolata la norma delle differenze relative
# mkt_data e' un DataFrame contenente tre series 'maturity', 'strike' e 'market price'
# curve e' un DataFrame contenente due series: 'time' e 'rate'
def loss_Call_VG(list_model_params, S0, mkt_data, curve, power, absrel):

    model_price_tmp = np.array(compute_VG_prices(list_model_params, S0, curve, mkt_data[['maturity','strike','type']]))
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