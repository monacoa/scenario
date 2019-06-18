import pandas as pd
import datetime
from win32com.client import constants as const

from sc_elab.excel_hook.xls_utils import drawBox, formatTestataCurva, findCalibrationPos, writeResultPandas
from sc_elab.core.anagrafica_dati import MaturityFromIntToString
from xls_Calibration import intestazioneCalibration

from DEF_intef import nameSheetCapFloorVolatilities, nameSheetBVolBoot

# -----------------------------------------------------------------------
# Funzione di lettura degli elementi identificativi della curva di sconto

def readFeaturesDiscCurve(input_dict):
   col = ['keys', 'TypeObject', 'Name']
   element_on_sheet = pd.DataFrame(columns=col)

   for k in input_dict.keys():
      item = input_dict[k]
      if u'Discount Factors' in item.loc[:,1].values:
          element_on_sheet=element_on_sheet.append({'keys': k,
                                                    'TypeObject': 'Discount Curve',
                                                    'Name': item.loc[0,0]}, ignore_index=True)
   return element_on_sheet


# ----------------------------------------------------
# Funzioni scrittura su foglio Excel

def intestazioneVCapFloor( xla, rng,  attributi, nCols = 2, title= "Cap Floor volatilities"):

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


def writeVCapFloorResOnXls(data, xla, ref_date, currency, contributor, tipo_modello, shift):

    r = findCalibrationPos(xla, nameSheetCapFloorVolatilities)

    #---
    #mi posiziono nella prima cella utile per scrivere i risultati del fitting
    #---

    # if option_print == 'matrix':
    #     Attributi = \
    #         {     "1. Date ref"    : ref_date
    #             , "2. Tipo Dato"   : 'VolCapFloor'
    #             , "3. Valore"      : 'MID'
    #             , "4. Contributor" : contributor
    #             , "5. Currency"    : currency
    #             , "6. Tipo modello": tipo_modello
    #             , "7. Tenor Floater " : '6M' #per il momento scritto a mano, se verra cambiato sara da modificare
    #             , "8. Rows"        : 'Expiry'
    #             , "9. Columns"     : 'Maturity'
    #               }
    Attributi = \
        {     "1. Date ref"    : ref_date
            , "2. OptionType"  : 'Vol Cap Floor'
            , "3. Contributor" : contributor
            , "4. Currency"    : currency
            , "5. Tipo Modello": tipo_modello
            , "6. Shift"       : shift
            , "7. Tipo dato"   : 'ATM'
              }

    r = intestazioneVCapFloor(xla = xla, rng = r, attributi = Attributi)
    r = writeResultPandas(xla = xla , rng = r, df = data, flagPrintColumns = True)


def write_VolCapFloor(xla, res, ref_date, currency , contributor, tipo_modello, shift):
    # casto l'output della pandas come float
    res = res.astype('float')
    res2 = pd.DataFrame()
    res2['Maturity'] = res['MaturityInt'].map(MaturityFromIntToString)
    res2['Strike']=res['Strike']
    res2['Volatility (x100)'] = res['ValoreMid']
    res2['Usage'] = 'Y'
    writeVCapFloorResOnXls(res2,  xla, ref_date, currency, contributor, tipo_modello, shift)


def writeBootstrapVolOnXls(xla, res, df_vols, df_curva):

    nameSheet = nameSheetBVolBoot
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

    Attributi = \
        {     "1. Date ref"              : df_vols.loc[(df_vols.loc[:,0]== 'Date ref'),1].values[0]
            , "2. OptionType"            : 'Vol Cap Floor'
            , "3. Curve Name"            : df_curva.iloc[0,0]
            , "4. Type Bootstrap"        : df_vols.loc[(df_vols.loc[:,0]== 'Tipo dato'),1].values[0]
            , "5. Contributor"           : df_vols.loc[(df_vols.loc[:,0]== 'Contributor'),1].values[0]
            , "6. Currency"              : df_vols.loc[(df_vols.loc[:,0]== 'Currency'),1].values[0]
            , "7. Tipo Modello"          : df_vols.loc[(df_vols.loc[:,0]== 'Tipo Modello'),1].values[0]
            , "8. Shift"                 : df_vols.loc[(df_vols.loc[:,0]== 'Shift'),1].values[0]
        }

    r = intestazioneCalibration(xla = xla, rng = r, attributi = Attributi, title = 'Bootstrap Cap Floor Volatilities')
    r = writeResultPandas(xla = xla , rng = r, df = res, flagPrintColumns = True)
    s = xla.Cells.Columns.AutoFit()