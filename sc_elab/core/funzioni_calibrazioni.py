from sc_elab.core.SwpCurve import BootstrappedCurve
from Tkinter import *
import tkMessageBox
import numpy as np
import datetime


# NOT USED!
def readCurveFromDataFrame(df):
    cc          = BootstrappedCurve()
    if type(cc) == BootstrappedCurve: cc.code     = (df.loc[0,0].split("_"))[1]

    cc.curr          = df.loc[df.loc[:, 0] == 'Currency',1].values[0]
    cc.type          = df.loc[df.loc[:, 0] == 'CurveType',1].values[0]

    cc.ref_date = df.loc[df.loc[:, 0] == 'Date Ref',1].values[0]

    cc.description   = df.loc[df.loc[:, 0] == 'Description',1].values[0]
    cc.download_type = df.loc[df.loc[:, 0] == 'Download Type',1].values[0]
    cc.quotation     = df.loc[df.loc[:, 0] == 'Quotation',1].values[0]
    cc.source        = df.loc[df.loc[:, 0] == 'Source',1].values[0]

    #imposto il mercato
    if cc.curr == "EUR":
        cc.cal = "de.eurex"
    elif cc.curr == "USD":
        cc.cal = 'us'
    elif cc.curr == 'GBP':
        cc.cal = 'uk'
    elif cc.curr == 'CAD':
        cc.cal = 'ca'
    else:
        cc.cal = 'us'


    row_sep = df.loc[df.loc[:, 0] == 'Date',].index[0]

    df_new = df.loc[(row_sep+1):df.shape[0],:]
    idx = df_new.loc[:,2] == 'Y'

    cc.boot_dates = df_new.loc[idx, 0].__array__()
    cc.boot_rates = df_new.loc[idx, 1].__array__()
    cc.fit_usage  = df_new.loc[idx, 2].__array__()

    tmp = df.loc[1,1].split(",")

    if tmp.__len__() !=4:
        tkMessageBox.showinfo("Error", "Non e' presente la capitalizzazione")
    else:
        #(0) Simple", "(1) Compounded", "(2) Continuous"
        if tmp[3] == '0':
            cc.capitalization = "SMP"
        elif tmp[3] == '1':
            cc.capitalization = "CMP"
        elif tmp[3] == '2':
            cc.capitalization = "CNT"
        else:
            tkMessageBox.showinfo("Error", "Non e' presente la capitalizzazione")

    return cc


def preProcessignCurve(df):

    # verifico la capitalizzazione dei tassi
    tmp_cap = df.loc[1,1].split(",")

    # leggo la data di riferimento
    ref_date = df.loc[df.loc[:, 0] == 'Date Ref',1].values[0]

    # elimino le intestazioni
    row_sep = df.loc[df.loc[:, 0] == 'Date',].index[0]
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

    if tmp_cap.__len__() != 4:
        tkMessageBox.showinfo("Error", "Non e' presente la capitalizzazione")
        capitalization = ""

    else:
        #(0) Simple", "(1) Compounded", "(2) Continuous"
        if tmp_cap[3] == '0':
            capitalization = "SMP"
            timeTmp = out.iloc[:, 0]
            out.iloc[:, 1] = -1.0/timeTmp*np.log(1.0/(1.0 + out.iloc[:, 1]*timeTmp))

        elif tmp_cap[3] == '1':
            capitalization = "CMP"
            out.iloc[:, 1] = np.log(1.0 + out.iloc[:, 1])

        elif tmp_cap[3] == '2':
            capitalization = "CNT"

        else:
            tkMessageBox.showinfo("Error", "Non e' presente la capitalizzazione")
            capitalization = ""

    out_to_fit = out.copy()
    print out_origin
    print '\n'
    print out_to_fit
    return out_origin ,out_to_fit , capitalization


def fromContinuousToCompost(r):
    rate = np.exp(r) - 1.
    return rate

################################## CIR ###########################################################
def loss_zc_model_cir(list_model_params, mkt_prices,order_difference):

    model_price_tmp = compute_zc_cir_rate(list_model_params, mkt_prices["TIME"])
    diff = np.absolute(model_price_tmp - mkt_prices["VALUE"])
    diff = np.power(diff,order_difference)

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


################################## VSCK ###########################################################
def loss_zc_model_vsck(list_model_params, mkt_prices,order_difference):

    model_price_tmp = compute_zc_vsck_rate(list_model_params, mkt_prices["TIME"])
    diff = np.absolute(model_price_tmp - mkt_prices["VALUE"])
    diff = np.power(diff,order_difference)

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


def computeCHI2(mkt, mdl):
    tmp = np.power(mkt - mdl,2)/np.absolute(mkt)
    #tmp = np.absolute(mkt - mdl)/np.absolute(mkt)
    #tmp = np.power(mkt - mdl,2)

    return tmp.sum()

def set_output_calibration(ff, mkt_value, type_cap):

    # creo la lista dei risultati ottimali
    list_model_params_opt = []
    list_model_params_opt.append(ff.x[0])
    list_model_params_opt.append(ff.x[1])
    list_model_params_opt.append(ff.x[2])
    list_model_params_opt.append(ff.x[3])

    mkt_value['VALUE_OPT'] = compute_zc_cir_rate(list_model_params_opt, mkt_value["TIME"])

    # converto i risultati in composto nel caso in cui in input lo siano
    if type_cap == 'CMP':
        mkt_value['VALUE_OPT'] = fromContinuousToCompost(mkt_value['VALUE_OPT'])

    # calcolo il chi quadro
    chi2 = computeCHI2(mkt=mkt_value["VALUE"], mdl=mkt_value['VALUE_OPT'])

    return list_model_params_opt,mkt_value,chi2
