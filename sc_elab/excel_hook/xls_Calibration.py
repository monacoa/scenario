from pyxll import xl_func, xl_app, xl_menu, xl_macro, xlcAlert
from sc_elab.excel_hook.xls_utils import drawBox, drawLine, formatTestataCurva, findRigthPlaceBootCurveSeg
import datetime
from sc_elab.excel_hook.DEF_intef import nameSheetCalib,nameSheetCalibRes
from    win32com.client import constants as const


def findCalibrationPos (xla, nameSheet):
    rangeInitial = "B2"
    distanza = 2
    sheet = xla.ActiveWorkbook.Sheets(nameSheet)
    r = sheet.Range(rangeInitial)
    col = r.Column
    i = 0

    rStart = xla.Range(xla.Cells(r.Row + distanza, col), xla.Cells(r.Row + distanza, col))

    while r.Value != None or rStart.Value != None:
        i = 1
        r = xla.Range(xla.Cells(r.Row + 1, col), xla.Cells(r.Row + 1, col))
        rStart = xla.Range(xla.Cells(r.Row + distanza, col), xla.Cells(r.Row + distanza, col))

    if i == 0:
        rStart = sheet.Range(rangeInitial)

    return rStart


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
            xla.Cells(topLeftRow + 1 + i, topLeftCol + 1).NumberFormat = "gg/MM/aaaa"
        xla.Cells(topLeftRow + 1 + i, topLeftCol + 1).HorizontalAlignment = const.xlCenter
        i+=1

    rangeStart = xla.Range(xla.Cells(topLeftRow + nRows + 1, topLeftCol),xla.Cells(topLeftRow + nRows + 1, topLeftCol))
    return rangeStart


def writeParameterCalibration( xla, rng , v_name , v_value,  dict, nCols = 4):

    nRows           = len(v_name)
    topLeftRow      = rng.Row
    topLeftCol      = rng.Column
    drawBox(xla, 3 , topLeftRow, topLeftCol,topLeftRow + nRows, topLeftCol + nCols - 1, 0)

    xla.Cells(topLeftRow , topLeftCol + 0).Value = 'Parameter'
    xla.Cells(topLeftRow , topLeftCol + 1).Value = 'Value'
    xla.Cells(topLeftRow , topLeftCol + 2).Value = 'Min'
    xla.Cells(topLeftRow , topLeftCol + 3).Value = 'Max'

    i = 0
    for k in v_name:
        xla.Cells(topLeftRow + 1+ i, topLeftCol+0).Value   = k
        xla.Cells(topLeftRow + 1+ i, topLeftCol+1).Value = v_value[i]
        xla.Cells(topLeftRow + 1+ i, topLeftCol+2).Value = dict[k]['min']
        xla.Cells(topLeftRow + 1+ i, topLeftCol+3).Value = dict[k]['max']
        i+=1

    rangeStart = xla.Range(xla.Cells(topLeftRow + nRows + 1, topLeftCol),xla.Cells(topLeftRow + nRows + 1, topLeftCol))
    return rangeStart


def writeResultPandas( xla, rng , df):

    nRows           = df.shape[0]
    nCols           = df.shape[1]
    topLeftRow      = rng.Row
    topLeftCol      = rng.Column
    drawBox(xla, 3 , topLeftRow, topLeftCol,topLeftRow + nRows, topLeftCol + nCols - 1, 0)

    for j in xrange(0,nCols):
        xla.Cells(topLeftRow , topLeftCol + j).Value = df.columns.values[j]

    for i in xrange(0,nRows):
        for j in xrange(0,nCols):
            xla.Cells(topLeftRow + 1+ i, topLeftCol+j).Value   = df.iloc[i,j]

    rangeStart = xla.Range(xla.Cells(topLeftRow + nRows + 1, topLeftCol),xla.Cells(topLeftRow + nRows + 1, topLeftCol))
    return rangeStart


def writeCalibrationResOnXls(model, W_class, xla, chi2, opt_dict, res):

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
        {     "0. Model"                 : model
            , "1. Date ref"              : ref_date
            , "2. Type Data Calibration" : W_class.set_mkt_ts.get()
            , "3. Name Curve"            : W_class.NameCurve.get()
            , "4. Name Option"           : W_class.NameOption.get()
            , "5. Name Time Series"      : W_class.NameTS.get()
            , "6. Type Calibration"      : W_class.mkt_calibration_type.get()
            , "7. Type Loss Function"    : W_class.loss_function_type.get()
            , "8. Chi-squared"           : chi2
        }

    row = r.Row
    col = r.Column
    r = intestazioneCalibration(xla = xla, rng = r, attributi = Attributi , title = model)
    r = writeParameterCalibration(xla = xla, rng = r, v_name = W_class.params_names, v_value = opt_dict,  dict = W_class.param_dict)
    r = writeResultPandas(xla = xla , rng = r, df = res)
