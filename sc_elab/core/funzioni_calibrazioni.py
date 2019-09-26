import tkMessageBox
import numpy as np
import pandas as pd
from scipy.special import ive
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
    out.columns = ['TIME','MKT']


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


def compute_zc_vsck_rate(list_model_params, T):
    r0 = list_model_params[0]
    kappa = list_model_params[1]
    theta = list_model_params[2]
    sigma = list_model_params[3]

    B0 = (1.0 / kappa) * (1.0 - np.exp(-kappa * T))
    g0 = (sigma * sigma) / (4.0 * kappa)
    G0 = (theta - sigma * sigma / (2.0 * kappa * kappa))
    A0 = np.exp(G0 * (B0 - T) - g0 * B0 * B0)

    model_rate = -(np.log(A0) - B0 * r0) / T

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


def computeCHI2(mkt, mdl):
    tmp = np.power(mkt - mdl,2)/np.absolute(mkt)

    return tmp.sum()
