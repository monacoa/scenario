import time

import numpy as np
import pandas   as pd
from sc_elab.core.utils_g2pp_newton import Pt_MKT,found_opt,price_swaption
from sc_elab.core.funzioni_calibrazioni import forwardSwap #, fromVolaToPrice, fromPriceToVola


def curve_fromCMPtoCNT(curve):
    CMPrates = curve['VALUE'].values
    CNTrates = np.log(1.+CMPrates)
    curve['VALUE'] = CNTrates
    return curve

def load_param(inputFile):
    dt = pd.read_excel(inputFile,sheet_name='Parameters')
    return dt.Value.__array__()


if __name__ == "__main__":

    # ----------- File di input ------------------------------
    inputFile = 'input/20190930_SwaptionTestInput.xlsx'
    # ----------- File di output -----------------------------
    outputFile = 'output/20190930_SwaptionTest_PY.xlsx'
    # ---------------------------------------------------------
    OutSheetBase = 'Base'
    OutSheetUp = 'Up'
    OutSheetDown = 'Down'

    # -------- Caricamento dati ---------------------------------------------------------
    # Attenzione! Assumiamo che le curve Eiopa e la curva di inflazione siano fornite nel REGIME COMPOSTO
    # e che i valori non siano percentuali ma assoluti
    curve_base = curve_fromCMPtoCNT(pd.read_excel(io=inputFile, sheet_name='EIOPA_Base'))
    curve_up = curve_fromCMPtoCNT(pd.read_excel(io=inputFile, sheet_name='EIOPA_Up'))
    curve_down = curve_fromCMPtoCNT(pd.read_excel(io=inputFile, sheet_name='EIOPA_Down'))

    # Carico i parametri del modello
    list_model_params_opt = load_param(inputFile)
    # -----------------------------------------------------------------------------------

    # Creo le liste con le curve per i test e i fogli di output
    curves = [curve_base, curve_up, curve_down]
    sheet_names = [OutSheetBase, OutSheetUp, OutSheetDown]


    #parametri dei dati in input
    tenr = 1. # tenor delle swaption in anni
    call_flag = 1.0  # S - K
    # call_flag = -1.0 # K - S
    discretization_integral = 50.

    step = 1./12.
    tempo = np.arange(0., 50., step) + step
    time_Analytic = range(0,51)

    for j in [0,1,2]:
        # inizializzo il DataFrame di output
        dataOut = pd.DataFrame()

        rf_times  = curves[j]['TIME'].values.astype(float)
        rf_values = curves[j]['VALUE'].values.astype(float)

        s = time.time()
        for i in tempo:
            dataOut.loc[i,"Time"] = i
            t_exp = float(i)
            t_mat = 10.
            swp_atm , ForwardAnnuity = forwardSwap(tenr, rf_times, rf_values, 0., t_exp, t_mat)
            dataOut.loc[i,"ForwardAnnuityPrice(t)"] = ForwardAnnuity
            dataOut.loc[i,"ForwardSwapRate(t)"] = swp_atm
            if round(i,3) in time_Analytic :
                dataOut.loc[i,"AnalyticSwaptionPrice(t)"] = price_swaption(list_model_params_opt, t_exp, t_mat, tenr, swp_atm, rf_times, rf_values, call_flag,discretization_integral)

        t = time.time() - s
        print "il tempo impiegato: %.6f" %t

        if j==0:
            with pd.ExcelWriter(outputFile,mode='w') as writer:
                dataOut.to_excel(excel_writer=writer,sheet_name=sheet_names[j],na_rep = "#N/D")
        else:
            with pd.ExcelWriter(outputFile,mode='a') as writer:
                dataOut.to_excel(excel_writer=writer,sheet_name=sheet_names[j],na_rep = "#N/D")