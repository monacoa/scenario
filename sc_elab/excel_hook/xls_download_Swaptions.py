import numpy as np
import datetime
from win32com.client import constants as const

from sc_elab.excel_hook.xls_Calibration import findCalibrationPos,writeResultPandas
from sc_elab.excel_hook.xls_utils import drawBox, formatTestataCurva

from DEF_intef import nameSheetScaricoSwaption


def intestazioneSwaptions( xla, rng,  attributi, nCols = 2, title= 'Matrix Swaption'):

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


def writeSwaptionsResOnXls(data, xla, ref_date, option_print):

    r = findCalibrationPos(xla, nameSheetScaricoSwaption)

    #---
    #mi posiziono nella prima cella utile per scrivere i risultati del fitting
    #---

    if option_print == 'matrix':
        Attributi = \
            {     "1. Date ref"    : ref_date
                , "2. Tipo Dato"   : 'VSwaption'
                , "3. Valore"      : 'MID'
                , "4. Contributor" : ''
                , "5. Currency"    : ''
                , "6. Rows"        : 'Tenor'
                , "7. Columns"     : 'Maturity'
                  }
    else:
        Attributi = \
            {     "1. Date ref"    : ref_date
                , "2. Tipo Dato"   : 'VSwaption'
                , "3. Valore"      : 'MID'
                , "4. Contributor" : ''
                , "5. Currency"    : ''
                  }

    r = intestazioneSwaptions(xla = xla, rng = r, attributi = Attributi)
    r = writeResultPandas(xla = xla , rng = r, df = data)


def write_Swaptions(xla, res, ref_date, option_print = 'matrix'):
    # casto l'output della pandas come float
    res = res.astype('float')
    if option_print == 'matrix':

        res3 = res.pivot_table(index='Tenor', columns='MaturityInt', values='ValoreMid',fill_value = -1)

        #cambio il format della matrice
        res3.columns = res3.columns.astype(np.float)
        res3.index = res3.index.astype(np.float)
        res3.reset_index(level='Tenor', inplace=True)

        writeSwaptionsResOnXls(res3, xla, ref_date, option_print)

    else:
        res['Usage'] = 'Y'
        writeSwaptionsResOnXls(res,  xla, ref_date, option_print)
