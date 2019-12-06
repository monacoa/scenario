from sc_elab.excel_hook.xls_utils import drawBox, drawLine, formatTestataCurva, findRigthPlaceBootCurveSeg
import datetime
from sc_elab.excel_hook.DEF_intef import nameSheetCalib,nameSheetCalibRes,FORMATT
from sc_elab.excel_hook.xls_utils import findCalibrationPos, writeResultPandas
from sc_elab.core.anagrafica_dati import MaturityFromIntToString


import pandas as pd
import numpy as np
from    win32com.client import constants as const



def intestazioneCalibration( xla, rng,  attributi, nCols = 2, title= 'Calibration'):

    nRows           = len(attributi.keys())
    topLeftRow      = rng.Row
    topLeftCol      = rng.Column
    drawBox            (xla, 3,topLeftRow, topLeftCol,topLeftRow + nRows, topLeftCol + nCols - 1, 0)
    formatTestataCurva (xla, topLeftRow, topLeftCol, nCols, title)

    kk = attributi.keys()
    kk.sort()
    i = 0
    for k in kk:
        xla.Cells(topLeftRow + 1+ i, topLeftCol).Value   = k[3:]
        xla.Cells(topLeftRow + 1+ i, topLeftCol+1).Value = attributi[k]
        if isinstance(attributi[k],datetime.datetime):
            xla.Cells(topLeftRow + 1 + i, topLeftCol + 1).NumberFormat = FORMATT
        xla.Cells(topLeftRow + 1 + i, topLeftCol + 1).HorizontalAlignment = const.xlCenter
        i+=1

    rangeStart = xla.Range(xla.Cells(topLeftRow + nRows + 1, topLeftCol),xla.Cells(topLeftRow + nRows + 1, topLeftCol))
    return rangeStart


def writeParameterCalibration( xla, rng , v_name , v_value,  dict, nCols = 5, rand = False):

    nRows           = len(v_name)
    topLeftRow      = rng.Row
    topLeftCol      = rng.Column
    drawBox(xla, 3 , topLeftRow, topLeftCol,topLeftRow + nRows, topLeftCol + nCols - 1, 0)

    xla.Cells(topLeftRow , topLeftCol + 0).Value = 'Parameter'
    xla.Cells(topLeftRow , topLeftCol + 1).Value = 'Value'
    xla.Cells(topLeftRow , topLeftCol + 2).Value = 'Initial guess'
    xla.Cells(topLeftRow , topLeftCol + 3).Value = 'Min'
    xla.Cells(topLeftRow , topLeftCol + 4).Value = 'Max'

    i = 0
    for k in v_name:
        xla.Cells(topLeftRow + 1+ i, topLeftCol+0).Value   = k
        xla.Cells(topLeftRow + 1+ i, topLeftCol+1).Value = v_value[i]
        if rand:
            xla.Cells(topLeftRow + 1 + i, topLeftCol + 2).Value = ''
        else:
            xla.Cells(topLeftRow + 1+ i, topLeftCol+2).Value = dict[k]['sv']
        xla.Cells(topLeftRow + 1+ i, topLeftCol+3).Value = dict[k]['min']
        xla.Cells(topLeftRow + 1+ i, topLeftCol+4).Value = dict[k]['max']
        i+=1

    rangeStart = xla.Range(xla.Cells(topLeftRow + nRows + 1, topLeftCol),xla.Cells(topLeftRow + nRows + 1, topLeftCol))
    return rangeStart


def writeCalibrationResOnXls(type_data, model, W_class, xla, chi2, opt_dict, res,capitalization_type = 'CNT'):

    nameSheet = nameSheetCalibRes
    try:
        s = xla.ActiveWorkbook.Sheets(nameSheet)
    except:
        s = xla.ActiveWorkbook.Sheets.Add()
        s.Name = nameSheet
    s.Activate()

    r = findCalibrationPos(xla, nameSheet)

    #---
    #mi posiziono nella prima cella utile per scrivere i risultati del fitting
    #---

    if type_data == 'MKT':

        df = W_class.CurveChosen
        ref_date = df.loc[df.loc[:, 0] == 'Date Ref',1].values[0]


        if W_class.loss_function_type.get() == 1:
            text_type_loss_function = "Euclidean absolute"

        elif  W_class.loss_function_type.get() == 2:
            text_type_loss_function = "Euclidean relative"

        elif W_class.loss_function_type.get() == 3:
            text_type_loss_function = "Manhattan absolute"

        elif W_class.loss_function_type.get() == 4:
            text_type_loss_function = "Manhattan relative"

        if model == 'CIR':

            res_feller_condition = 2 * opt_dict[1] * opt_dict[2] > (np.power(opt_dict[3],2))

            Attributi = \
                {     "0. Model"                 : model
                    , "1. Date ref"              : ref_date
                    , "2. Type Data Calibration" : W_class.set_mkt_ts.get()
                    , "3. Name Curve"            : W_class.NameCurve.get()
                    , "4. Name Option"           : W_class.NameOption.get()
                    , "5. Type Calibration"      : W_class.mkt_calibration_type.get()
                    , "6. Type Loss Function"    : text_type_loss_function
                    , "7. Chi-squared"           : chi2
                    , "8. Interest rate Type"    : capitalization_type
                    , "9. Feller condition"      : str(res_feller_condition)
                }

        else:
            Attributi = \
                {     "0. Model"                 : model
                    , "1. Date ref"              : ref_date
                    , "2. Type Data Calibration" : W_class.set_mkt_ts.get()
                    , "3. Name Curve"            : W_class.NameCurve.get()
                    , "4. Name Option"           : W_class.NameOption.get()
                    , "5. Type Calibration"      : W_class.mkt_calibration_type.get()
                    , "6. Type Loss Function"    : text_type_loss_function
                    , "7. Chi-squared"           : chi2
                    , "8. Interest rate Type"    : capitalization_type
                }

    else:
        if model == 'CIR':

            res_feller_condition = 2 * opt_dict[1] * opt_dict[2] > (np.power(opt_dict[3],2))

            Attributi = \
                {"0. Model": model
                    , "1. Name Time Series"     : W_class.NameTS.get()
                    , "2. Start date "          : W_class.TS_dateMIN.strftime("%d/%m/%Y")
                    , "3. End date"             : W_class.TS_dateMAX.strftime("%d/%m/%Y")
                    , "4. Chi-squared"          : chi2
                    , "5. Feller condition"     : str(res_feller_condition)
                }
        else:
            Attributi = \
                {     "0. Model"                 : model
                    , "1. Name Time Series"      : W_class.NameTS.get()
                    , "2. Start date "           : W_class.TS_dateMIN.strftime("%d/%m/%Y")
                    , "3. End date"              : W_class.TS_dateMAX.strftime("%d/%m/%Y")
                    , "4. Chi-squared"           : chi2
                }

    rand_flag = False
    if W_class.nTime.get() > 1 : rand_flag = True

    r = intestazioneCalibration(xla = xla, rng = r, attributi = Attributi , title = model)
    r = writeParameterCalibration(xla = xla, rng = r, v_name = W_class.params_names, v_value = opt_dict,  dict = W_class.param_dict, rand=rand_flag)
    r = writeResultPandas(xla = xla , rng = r, df = res, flagPrintColumns = True)
    s = xla.Cells.Columns.AutoFit()



def writeDividendsResOnXls(title, W_class, xla, res, dividend_type = 'CNT'):

    nameSheet = nameSheetCalibRes
    try:
        s = xla.ActiveWorkbook.Sheets(nameSheet)
    except:
        s = xla.ActiveWorkbook.Sheets.Add()
        s.Name = nameSheet
    s.Activate()

    r = findCalibrationPos(xla, nameSheet)

    #---
    #mi posiziono nella prima cella utile per scrivere i risultati del fitting
    #---

    df = W_class.CurveChosen
    ref_date = df.loc[df.loc[:, 0] == 'Date Ref',1].values[0]

    Attributi = \
        {     "1. Date ref"              : ref_date
            , "2. Type Data Calibration" : W_class.set_mkt_ts.get()
            , "3. Name Curve"            : W_class.NameCurve.get()
            , "4. Name Option"           : W_class.NameOption.get()
            , "5. Dividend rate Type"    : dividend_type
        }

    r = intestazioneCalibration(xla = xla, rng = r, attributi = Attributi , title = title)
    r = writeResultPandas(xla = xla , rng = r, df = res, flagPrintColumns = True)
    s = xla.Cells.Columns.AutoFit()


def writeVolResOnXls(title, W_class, xla, res):

    nameSheet = nameSheetCalibRes
    try:
        s = xla.ActiveWorkbook.Sheets(nameSheet)
    except:
        s = xla.ActiveWorkbook.Sheets.Add()
        s.Name = nameSheet
    s.Activate()

    r = findCalibrationPos(xla, nameSheet)

    #---
    #mi posiziono nella prima cella utile per scrivere i risultati del fitting
    #---

    df = W_class.CurveChosen
    ref_date = df.loc[df.loc[:, 0] == 'Date Ref',1].values[0]

    Attributi = \
        {     "1. Date ref"              : ref_date
            , "2. Type Volatility"       : 'Black-Scholes vol from VG'
        }

    r = intestazioneCalibration(xla = xla, rng = r, attributi = Attributi , title = title)
    r = writeResultPandas(xla=xla, rng=r, df=res, flagPrintColumns=True)
    s = xla.Cells.Columns.AutoFit()


def writeTemplateCalibration(xla, nameSheet):

    r = findCalibrationPos(xla, nameSheet)

    #---
    #mi posiziono nella prima cella utile per scrivere i risultati del fitting
    #---

    ##########################
    ## CURVE - Date
    ##########################

    mat = pd.DataFrame()

    mat['Date']  = pd.date_range(start=datetime.datetime.now().date() - datetime.timedelta(7), periods=7)
    mat['Value (x100)'] = np.zeros(7)
    mat['Usage'] = 'Y'

    Attributi = \
        {       "0. CurveType"           : 'Zero coupon RF, Inflation, Real'
            ,   "1. Interest rate Type"  : 'SMP, CMP, CNT'
            ,   "2. Date Ref"            :  datetime.datetime.now().strftime("%m/%d/%Y")
        }

    r = intestazioneCalibration(xla=xla, rng=r, attributi=Attributi, title='Template Calibration Curve - Date')
    r = writeResultPandas(xla=xla, rng=r, df=mat, flagPrintColumns=True)
    xla.Cells.ColumnWidth = 18


    ##########################
    ## CURVE - Time
    ##########################

    mat = pd.DataFrame()

    mat['Times']  = np.arange(1,4.5,step=0.5)
    mat['Value (x100)'] = np.zeros(7)
    mat['Usage'] = 'Y'

    Attributi = \
        {       "0. CurveType"           : 'Zero coupon RF, Inflation, Real'
            ,   "1. Interest rate Type"  : 'SMP, CMP, CNT'
            ,   "2. Date Ref"            :  datetime.datetime.now().strftime("%m/%d/%Y")
        }

    r = intestazioneCalibration(xla=xla, rng=r, attributi=Attributi, title='Template Calibration Curve - Times')
    r = writeResultPandas(xla=xla, rng=r, df=mat, flagPrintColumns=True)
    xla.Cells.ColumnWidth = 18

    ##########################
    ## TIME SERIES
    ##########################

    mat = pd.DataFrame()

    mat['Date'] = pd.date_range(start=datetime.datetime.now().date() - datetime.timedelta(7),periods=7)
    mat['Value'] = np.zeros(7)
    mat['Usage'] = 'Y'

    Attributi = \
        {     "0. TSType"      : 'Time Series'
        }

    r = intestazioneCalibration(xla = xla, rng = r, attributi = Attributi , title = 'Template Calibration Time Series')
    r = writeResultPandas(xla = xla , rng = r, df = mat, flagPrintColumns = True)
    xla.Cells.ColumnWidth = 18

    ##########################
    ## OPTION - Swaption
    ##########################

    mat = pd.DataFrame()

    mat['Expiry']   = np.arange(1,8,step=1) * 360.
    mat['Maturity'] = 10. * 360.
    mat['Expiry']   = mat['Expiry'].map(MaturityFromIntToString)
    mat['Maturity'] = mat['Maturity'].map(MaturityFromIntToString)
    mat['Value (x100)'] = np.zeros(7)
    mat['Shift'] = np.zeros(7)
    mat['Usage'] = 'Y'

    Attributi = {
             "0. Date Ref"    : datetime.datetime.now().strftime("%m/%d/%Y")
            ,"1. OptionType"  : 'Swaption'
            ,"2. SwaptionType": 'Payer, Receiver'
            ,"3. Type value"  : 'Volatility'
            ,"4. Type model"  : 'No Shifted'
            ,"5. Tenor swap"  : '6M'
    }

    r = intestazioneCalibration(xla=xla, rng=r, attributi=Attributi, title='Template Calibration Swaption')
    r = writeResultPandas(xla=xla, rng=r, df=mat, flagPrintColumns=True)
    xla.Cells.ColumnWidth = 18


    ##########################
    ## Cap Floor
    ##########################

    mat = pd.DataFrame()

    mat['Time']   = np.arange(1,8,step=1)
    mat['Strike (x100)'] = np.repeat(0.9,7)
    mat['Value (x100)'] = np.zeros(7)
    mat['Usage'] = 'Y'

    Attributi = {
             "0. Date Ref"           : datetime.datetime.now().strftime("%m/%d/%Y")
            ,"1. OptionType"         : 'Vol Caps, Vol Caplets, Caplets, Caps'
            # ,"2. Type value"         : 'Price, Volatility'
            ,"3. Shift"              : 0
    }

    r = intestazioneCalibration(xla=xla, rng=r, attributi=Attributi, title='Template Calibration Cap Floor')
    r = writeResultPandas(xla=xla, rng=r, df=mat, flagPrintColumns=True)
    xla.Cells.ColumnWidth = 18

    #######################
    ## OPTION - Opzione
    ##########################

    mat = pd.DataFrame()

    mat['Maturity']   = np.arange(1,8,step=1)
    mat['Strike'] = np.repeat(0.9,7)
    mat['Price'] = np.zeros(7)
    mat['Type'] = 'CALL'
    mat['Usage'] = 'Y'

    Attributi = {
             "0. Date Ref"           : datetime.datetime.now().strftime("%m/%d/%Y")
            ,"1. OptionType"         : 'PUT / CALL'
            ,"2. Initial Price"      : 100
            ,"3. Strike"             : 100
            ,"4. Mat"                : 2
    }

    r = intestazioneCalibration(xla=xla, rng=r, attributi=Attributi, title='Template Calibration Option')
    r = writeResultPandas(xla=xla, rng=r, df=mat, flagPrintColumns=True)
    xla.Cells.ColumnWidth = 18

    #######################
    ## INFLATION OPTION
    ##########################

    mat = pd.DataFrame()

    mat['Maturity'] = np.arange(1, 8, step=1)
    mat['Value'] = np.zeros(7)
    mat['Usage'] = 'Y'

    Attributi = {
        "0. Date Ref": datetime.datetime.now().strftime("%m/%d/%Y")
        , "1. OptionType": 'Cap/Floor'
        , "2. Type value": 'Price/Volatility'
        , "3. tenorOplet": 12
    }

    r = intestazioneCalibration(xla=xla, rng=r, attributi=Attributi, title='Template Calibration Inflation Option')
    r = writeResultPandas(xla=xla, rng=r, df=mat, flagPrintColumns=True)
    xla.Cells.ColumnWidth = 18

    ##########################
    # Vol from surface
    ##########################

    mat = pd.DataFrame()
    mat['Maturity'] = np.arange(1, 8, step=1)
    mat['Strike'] = 9999*np.ones(7)
    mat['Usage'] = 'Y'

    Attributi = {
        "0. Date Ref":  datetime.datetime.now().strftime("%m/%d/%Y")
        , "1. ElabType": "Vol from Surface"
    }

    r = intestazioneCalibration(xla=xla, rng=r, attributi=Attributi, title='Template Vol from Surface')
    r = writeResultPandas(xla=xla, rng=r, df=mat, flagPrintColumns=True)
    xla.Cells.ColumnWidth = 18