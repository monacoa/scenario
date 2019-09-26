import numpy as np
import pandas as pd
import sc_elab.core.funzioni_calibrazioni as JY

def curve_fromCMPtoCNT(curve):
    CMPrates = curve['VALUE'].values
    CNTrates = np.log(1.+CMPrates)
    curve['VALUE'] = CNTrates

    return curve

def parameters_fromDFtoDict(parameters):
    dict = {}
    parameters_array = parameters['Parameter'].values.astype(str)
    for p in parameters_array:
        value = parameters.loc[parameters['Parameter']==p,'Value'].values[0]
        dict[p]=value
    return dict


if __name__ == "__main__":

    # ----------- File di input ------------------------------
    inputFile = 'input/20190628_InflationOptionTestInput.xlsx'
    # ----------- File di output -----------------------------
    outputFile = 'output/20190628_InflationOptionTest_PY.xlsx'
    #---------------------------------------------------------
    OutSheetBase = 'Base'
    OutSheetUp   = 'Up'
    OutSheetDown = 'Down'

    # Attenzione! Assumiamo che le curve Eiopa e la curva di inflazione siano fornite nel REGIME COMPOSTO
    # e che i valori non siano percentuali ma assoluti
    curve_base = curve_fromCMPtoCNT(pd.read_excel(io=inputFile,sheet_name='EIOPA_Base'))
    curve_up = curve_fromCMPtoCNT(pd.read_excel(io=inputFile, sheet_name='EIOPA_Up'))
    curve_down = curve_fromCMPtoCNT(pd.read_excel(io=inputFile, sheet_name='EIOPA_Down'))
    curve_infl = curve_fromCMPtoCNT(pd.read_excel(io=inputFile, sheet_name='Inflation'))

    # Le seguenti righe calcolano le curve reali a partire dalle curve nominali e quella di inflazione
    curve_real_base = JY.createRealCurve(curve_base, curve_infl)
    curve_real_up = JY.createRealCurve(curve_up, curve_infl)
    curve_real_down = JY.createRealCurve(curve_down, curve_infl)

    # Vengono letti i parametri del modello dal foglio Excel di input e creato il corrispondente dizionario
    model_parameters = pd.read_excel(io=inputFile, sheet_name='Parameters')
    param = parameters_fromDFtoDict(model_parameters)

    # creo le liste di curve su cui eseguire il test e la lista di fogli su cui scrivere
    nominal_curves = [curve_base,curve_up,curve_down]
    real_curves    = [curve_real_base,curve_real_up,curve_real_down]
    sheet_names    = [OutSheetBase,OutSheetUp,OutSheetDown]

    # imposto alcuni parametri per il test
    I0 = 1.
    t0 = 0.

    param['Strike'] = 1.
    param['I0']     = I0
    param['omega']  = 1.
    param['t0']     = t0


    step = 1./12.
    time = np.arange(0., 50., step) + step

    # Eseguo il test e scrivo i risultati
    for j in range(0,3):
        data = pd.DataFrame()
        for i in time:
            data.loc[i,"Time"] = i
            te = i
            param['te'] = te
            data.loc[i,"AnalyticOptionPrice(t)"] = JY.computeDfCurve(nominal_curves[j], te) * JY.analytical_option_jy(param, nominal_curves[j], real_curves[j])

        if j==0:
            with pd.ExcelWriter(outputFile,mode='w') as writer:
                data.to_excel(excel_writer=writer,sheet_name=sheet_names[j])
        else:
            with pd.ExcelWriter(outputFile,mode='a') as writer:
                data.to_excel(excel_writer=writer,sheet_name=sheet_names[j])