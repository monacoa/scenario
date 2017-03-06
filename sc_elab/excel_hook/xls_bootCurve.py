from pyxll import xlcAlert
import sys
import datetime
from win32com.client import constants as const
from Tkinter import *

from sc_elab.core.SwpCurve import BootstrappedCurve
from xls_swapCurve import intestazioneSwapCurveSegmenti
from xls_utils import drawBox,drawLine, findRigthPlaceBootCurveSeg

from DEF_intef import FORMATT, nameSheetBootstrap



def writeBootstrapResOnXls(crv, xla, str_boot_opt, res, codeL, codeR):
    nameSheet = nameSheetBootstrap
    try:
        s = xla.ActiveWorkbook.Sheets(nameSheet)
    except:
        s = xla.ActiveWorkbook.Sheets.Add()
        s.Name = nameSheet
    s.Activate()
    rangeStart  = "B2"
    distCurve   = 1
    r           = s.Range(rangeStart)
    r           = findRigthPlaceBootCurveSeg(xla, r, distCurve, "v")
    str_segms   = res['CodiceCurva']
    # -----
    Attributi_1 = \
        { "Date Ref"        : crv.ref_date
        , "Description"     : crv.description
        , "Currency"        : crv.curr
        , "Download Type"   : crv.download_type
        , "Quotation"       : crv.quotation
        , "Source"          : crv.source
        , "CurveType"       : crv.type
        , "Return"          : "Zero Coupon"
        , "Node Type"       : "Spot"
        , "Boot Optons"     : str_boot_opt
        , "Segms"           : str_segms
        }

    Attributi_2 = \
        { "Date Ref"        : crv.ref_date
        , "Description"     : crv.description
        , "Currency"        : crv.curr
        , "Download Type"   : crv.download_type
        , "Quotation"       : crv.quotation
        , "Source"          : crv.source
        , "CurveType"       : crv.type
        , "Return"          : "-"
        , "Node Type"       : "Discount Factors"
        , "Boot Optons"     : str_boot_opt
        , "Segms"           : str_segms
        }

    r2 = s.Range(xla.Cells(r.Row, r.Column + 3), xla.Cells(r.Row, r.Column + 3))
    ra = intestazioneSwapCurveSegmenti(xla, s, r , Attributi_1, nCols=2, text = codeL+"ST"+codeR+"_"+res['CodiceCurva'])
    rb = intestazioneSwapCurveSegmenti(xla, s, r2, Attributi_2, nCols=2, text = codeL+"DF"+codeR+"_"+res['CodiceCurva'])
    r  = s.Range(ra)

    topLeftRow = r.Row
    topLeftCol = r.Column

    xla.Cells(topLeftRow - 1, topLeftCol).Value = "Date"
    xla.Cells(topLeftRow - 1, topLeftCol).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 1).Value = "Value"
    xla.Cells(topLeftRow - 1, topLeftCol + 1).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow - 1, topLeftCol + 3).Value = "Date"
    xla.Cells(topLeftRow - 1, topLeftCol + 3).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 4).Value = "Value"
    xla.Cells(topLeftRow - 1, topLeftCol + 4).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 5).Value = "Usage"
    xla.Cells(topLeftRow - 1, topLeftCol + 5).HorizontalAlignment = const.xlCenter

    nRows = len(res['DiscountFactors'])

    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol, topLeftRow + nRows - 1, topLeftCol + 1)
    drawLine(xla, topLeftRow - 1, topLeftCol, topLeftRow - 1, topLeftCol + 1, "o", const.xlThin)

    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol + 3, topLeftRow + nRows - 1, topLeftCol + 5)
    drawLine(xla, topLeftRow - 1, topLeftCol + 3, topLeftRow - 1, topLeftCol + 5, "o", const.xlThin)
    for i in range(nRows):
        value = res['DiscountFactors'][i]
        date = res['DateScadenza'][i]
        rate = res['TassiZC'][i]

        xla.Cells(topLeftRow + i, topLeftCol).Value = date
        xla.Cells(topLeftRow + i, topLeftCol).NumberFormat = FORMATT
        xla.Cells(topLeftRow + i, topLeftCol).HorizontalAlignment = const.xlCenter
        xla.Cells(topLeftRow + i, topLeftCol + 1).Value = rate
        xla.Cells(topLeftRow + i, topLeftCol + 1).NumberFormat = "0.00"
        xla.Cells(topLeftRow + i, topLeftCol + 1).HorizontalAlignment = const.xlCenter

        xla.Cells(topLeftRow + i, topLeftCol + 3).Value = date
        xla.Cells(topLeftRow + i, topLeftCol + 3).NumberFormat = FORMATT
        xla.Cells(topLeftRow + i, topLeftCol + 3).HorizontalAlignment = const.xlCenter
        xla.Cells(topLeftRow + i, topLeftCol + 4).Value = value
        xla.Cells(topLeftRow + i, topLeftCol + 4).NumberFormat = "0.00000"
        xla.Cells(topLeftRow + i, topLeftCol + 4).HorizontalAlignment = const.xlCenter
        xla.Cells(topLeftRow + i, topLeftCol + 5).Value = "Y"
        xla.Cells(topLeftRow + i, topLeftCol + 5).HorizontalAlignment = const.xlCenter

def readIntestazioneBootstrap(xla , r , cc):
    row              = r.Row
    col              = r.Column


    if type(cc) == BootstrappedCurve : cc.code          = ((r.Value).split("_"))[1]


    cc.curr          = xla.Range(xla.Cells(row + 2, col + 1), xla.Cells(row + 2, col + 1)).Value
    cc.type          = xla.Range(xla.Cells(row + 3, col + 1), xla.Cells(row + 3, col + 1)).Value

    dd = xla.Range(xla.Cells(row + 4, col + 1), xla.Cells(row + 4, col + 1)).Value
    cc.ref_date      = datetime.date(year = dd.year, month = dd.month, day = dd.day)

    cc.description   = xla.Range(xla.Cells(row + 5, col+1), xla.Cells(row+5, col+1)).Value
    cc.download_type = xla.Range(xla.Cells(row + 6, col + 1), xla.Cells(row + 6, col + 1)).Value
    cc.quotation     = xla.Range(xla.Cells(row + 5, col + 1), xla.Cells(row + 8, col + 1)).Value
    cc.source        = xla.Range(xla.Cells(row + 6, col + 1), xla.Cells(row + 11, col + 1)).Value


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

    return  xla.Range(xla.Cells(row + 13, col+3), xla.Cells(row + 13, col+3))


def readRatesAndDiscountFactors(xla , r , cc):
    i   = 0
    row = r.Row
    col = r.Column
    print "valore inziale di r:", r.Value, r.Address
    while (True):
        if (xla.Range(xla.Cells(row + i , col), xla.Cells(row + i, col)).Value == None) : break

        dd          = xla.Range(xla.Cells(row + i , col), xla.Cells(row + i, col)).Value
        data        = datetime.date(year=dd.year, month=dd.month, day=dd.day)
        rate        = xla.Range(xla.Cells(row + i , col-2), xla.Cells(row + i, col-2)).Value
        df          = xla.Range(xla.Cells(row + i , col+1), xla.Cells(row + i, col+1)).Value
        use         = xla.Range(xla.Cells(row + i , col+2), xla.Cells(row + i, col+2)).Value
        cc.boot_dates.append(data)
        cc.boot_df.append(df)
        cc.boot_rates.append(rate)
        cc.fit_usage.append(use)
        i += 1


def findBootstrappedCurveFromPos (xla, nameSheet, pos):
    rangeStart = "B2"
    distanza = 2
    sheet = xla.ActiveWorkbook.Sheets(nameSheet)
    r = sheet.Range(rangeStart)
    nCols = r.Columns.Count
    row = r.Row
    col = r.Column
    i = 1
    while i < pos:
        while (r.Value != None):
            nCols = r.Columns.Count
            row = r.Row
            col = r.Column
            r = xla.Range(xla.Cells(row + 1, col), xla.Cells(row + 1, col))

        r = xla.Range(xla.Cells(row + distanza, col), xla.Cells(row + distanza, col))
        if r.Value != None:
            i += 1
    return r


def readBootstrappedCurveFromXls(xla, des, pos, nameSheet):

    r           = findBootstrappedCurveFromPos (xla, nameSheet, pos)
    nomeCurva   = r.Value
    cc          = BootstrappedCurve()
    r           = readIntestazioneBootstrap(xla, r,cc)

    readRatesAndDiscountFactors(xla, r, cc)

    return cc
