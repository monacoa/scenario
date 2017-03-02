from pyxll import xlcAlert
import sys
import datetime
from win32com.client import constants as const
from Tkinter import *


from xls_utils import drawLine, drawBox
from DEF_intef import fitt_translate, FORMAT, nameSheetBootstrap, nameSheetCurve
from xls_swapCurve import findCurveFromPos
from xls_bootCurve import intestazioneSwapCurveSegmenti, findBootstrappedCurveFromPos



def writeFittingResLinear(xla, s, r, Attributi, res):
    ra = intestazioneSwapCurveSegmenti(xla, "", r, Attributi, nCols=2)
    r = s.Range(ra)
    topLeftRow = r.Row
    topLeftCol = r.Column
    xla.Cells(topLeftRow - 1, topLeftCol).Value = "Date"
    xla.Cells(topLeftRow - 1, topLeftCol).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 1).Value = "a"
    xla.Cells(topLeftRow - 1, topLeftCol + 1).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 2).Value = "b"
    xla.Cells(topLeftRow - 1, topLeftCol + 2).HorizontalAlignment = const.xlCenter
    nRows = len(res['Dates'])
    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol, topLeftRow + nRows - 1, topLeftCol + 2)
    drawLine(xla, topLeftRow - 1, topLeftCol, topLeftRow - 1, topLeftCol + 2, "o", const.xlThin)
    for i in range(nRows):
        date = res['Dates'][i]
        xla.Cells(topLeftRow + i, topLeftCol).Value = date
        xla.Cells(topLeftRow + i, topLeftCol).NumberFormat = FORMAT
        xla.Cells(topLeftRow + i, topLeftCol).HorizontalAlignment = const.xlCenter
        if (i > 0):
            a = res['a'][i-1]
            b = res['b'][i-1]
            xla.Cells(topLeftRow + i, topLeftCol + 1).Value = a
            xla.Cells(topLeftRow + i, topLeftCol + 1).NumberFormat = "0.0000"
            xla.Cells(topLeftRow + i, topLeftCol + 1).HorizontalAlignment = const.xlCenter
            xla.Cells(topLeftRow + i, topLeftCol + 2).Value = b
            xla.Cells(topLeftRow + i, topLeftCol + 2).NumberFormat = "0.0000"
            xla.Cells(topLeftRow + i, topLeftCol + 2).HorizontalAlignment = const.xlCenter


def writeFittingResAVD(xla, s, r, Attributi, res):
    ra = intestazioneSwapCurveSegmenti(xla, "", r, Attributi, nCols=2)
    r = s.Range(ra)
    topLeftRow = r.Row
    topLeftCol = r.Column
    xla.Cells(topLeftRow - 1, topLeftCol).Value = "Date"
    xla.Cells(topLeftRow - 1, topLeftCol).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 1).Value = "a"
    xla.Cells(topLeftRow - 1, topLeftCol + 1).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 2).Value = "b"
    xla.Cells(topLeftRow - 1, topLeftCol + 2).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 3).Value = "c"
    xla.Cells(topLeftRow - 1, topLeftCol + 3).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 4).Value = "d"
    xla.Cells(topLeftRow - 1, topLeftCol + 4).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 5).Value = "e"
    xla.Cells(topLeftRow - 1, topLeftCol + 5).HorizontalAlignment = const.xlCenter

    nRows = len(res['Dates'])
    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol, topLeftRow + nRows - 1, topLeftCol + 5)
    drawLine(xla, topLeftRow - 1, topLeftCol, topLeftRow - 1, topLeftCol + 5, "o", const.xlThin)
    for i in range(nRows):
        date = res['Dates'][i]
        xla.Cells(topLeftRow + i, topLeftCol).Value = date
        xla.Cells(topLeftRow + i, topLeftCol).NumberFormat = FORMAT
        xla.Cells(topLeftRow + i, topLeftCol).HorizontalAlignment = const.xlCenter
        if (i > 0):
            a = res['a'][i - 1]
            b = res['b'][i - 1]
            c = res['c'][i - 1]
            d = res['d'][i - 1]
            e = res['e'][i - 1]
            xla.Cells(topLeftRow + i, topLeftCol + 1).Value = a
            xla.Cells(topLeftRow + i, topLeftCol + 1).NumberFormat = "0.0000"
            xla.Cells(topLeftRow + i, topLeftCol + 1).HorizontalAlignment = const.xlCenter
            xla.Cells(topLeftRow + i, topLeftCol + 2).Value = b
            xla.Cells(topLeftRow + i, topLeftCol + 2).NumberFormat = "0.0000"
            xla.Cells(topLeftRow + i, topLeftCol + 2).HorizontalAlignment = const.xlCenter
            xla.Cells(topLeftRow + i, topLeftCol + 3).Value = c
            xla.Cells(topLeftRow + i, topLeftCol + 3).NumberFormat = "0.0000"
            xla.Cells(topLeftRow + i, topLeftCol + 3).HorizontalAlignment = const.xlCenter
            xla.Cells(topLeftRow + i, topLeftCol + 4).Value = d
            xla.Cells(topLeftRow + i, topLeftCol + 4).NumberFormat = "0.0000"
            xla.Cells(topLeftRow + i, topLeftCol + 4).HorizontalAlignment = const.xlCenter
            xla.Cells(topLeftRow + i, topLeftCol + 5).Value = e
            xla.Cells(topLeftRow + i, topLeftCol + 5).NumberFormat = "0.0000"
            xla.Cells(topLeftRow + i, topLeftCol + 5).HorizontalAlignment = const.xlCenter


def writeFittingResSVE(xla, s, r, Attributi, res):
    ra = intestazioneSwapCurveSegmenti(xla, "", r, Attributi, nCols=2)
    r = s.Range(ra)

    topLeftRow = r.Row
    topLeftCol = r.Column
    xla.Cells(topLeftRow - 1, topLeftCol).Value = "const1"
    xla.Cells(topLeftRow - 1, topLeftCol).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 1).Value = "const2"
    xla.Cells(topLeftRow - 1, topLeftCol + 1).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 2).Value = "beta 0"
    xla.Cells(topLeftRow - 1, topLeftCol + 2).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 3).Value = "beta 1"
    xla.Cells(topLeftRow - 1, topLeftCol + 3).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 4).Value = "beta 2"
    xla.Cells(topLeftRow - 1, topLeftCol + 4).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 5).Value = "beta 3"
    xla.Cells(topLeftRow - 1, topLeftCol + 5).HorizontalAlignment = const.xlCenter

    nRows = 1
    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol, topLeftRow + nRows - 1, topLeftCol + 5)
    drawLine(xla, topLeftRow - 1, topLeftCol, topLeftRow - 1, topLeftCol + 5, "o", const.xlThin)

    const1 = res['const1']
    const2 = res['const2']
    beta0  = res['beta0']
    beta1  = res['beta1']
    beta2  = res['beta2']
    beta3  = res['beta3']

    xla.Cells(topLeftRow , topLeftCol ).Value = const1
    xla.Cells(topLeftRow , topLeftCol ).NumberFormat = "0.0000"
    xla.Cells(topLeftRow , topLeftCol ).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow, topLeftCol+1).Value = const2
    xla.Cells(topLeftRow, topLeftCol+1).NumberFormat = "0.0000"
    xla.Cells(topLeftRow, topLeftCol+1).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow, topLeftCol+2).Value = beta0
    xla.Cells(topLeftRow, topLeftCol+2).NumberFormat = "0.0000"
    xla.Cells(topLeftRow, topLeftCol+2).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow, topLeftCol+3).Value = beta1
    xla.Cells(topLeftRow, topLeftCol+3).NumberFormat = "0.0000"
    xla.Cells(topLeftRow, topLeftCol+3).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow, topLeftCol+4).Value = beta2
    xla.Cells(topLeftRow, topLeftCol+4).NumberFormat = "0.0000"
    xla.Cells(topLeftRow, topLeftCol+4).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow, topLeftCol+5).Value = beta3
    xla.Cells(topLeftRow, topLeftCol+5).NumberFormat = "0.0000"
    xla.Cells(topLeftRow, topLeftCol+5).HorizontalAlignment = const.xlCenter



def writeFittingBootResOnXls(crv, xla, opt_dict, res, pos):
    nameSheet = nameSheetBootstrap
    try:
        s = xla.ActiveWorkbook.Sheets(nameSheet)
    except:
        s = xla.ActiveWorkbook.Sheets.Add()
        s.Name = nameSheet
    s.Activate()

    r = findBootstrappedCurveFromPos(xla, nameSheet,pos)
    print "------------------indirizzo trovato:", r.Value, r.Address
    #---
    #mi posiziono nella prima cella utile per scrivere i risultati del fitting
    #---
    Attributi = \
        {     "Date Ref"     : crv.ref_date
            , "Segms:"       : "-"
            , "Al"           : "-"
            , "Description"  : crv.description
            , "Currency"     : crv.curr
            , "Download Type": crv.download_type
            , "Quotation"    : crv.quotation
            , "CurveType"    : crv.type
            , "Return"       : "-"
            , "Node Type"    : "-"
            , "Interp. Model": fitt_translate[opt_dict['interp']]
        }

    row = r.Row
    col = r.Column
    r = xla.Range(xla.Cells(row+len(Attributi.keys())+1, col ),(xla.Cells(row+len(Attributi.keys())+1, col)))

    while (r.Value!= None):
        row = r.Row
        col = r.Column
        r = xla.Range(xla.Cells(row, col + 1), xla.Cells(row, col + 1))
        if (r.Value == None):
            row = r.Row
            col = r.Column
            r = xla.Range(xla.Cells(row, col + 1), xla.Cells(row, col + 1))
    print "------------------posizione calcolata:", r.Value, r.Address

    row = r.Row
    col = r.Column
    r = xla.Range(xla.Cells(row -( len(Attributi.keys()) +1), col), (xla.Cells(row -(len(Attributi.keys()) + 1), col)))

    print opt_dict
    print opt_dict['interp'], type(opt_dict['interp'])
    if      opt_dict['interp'] == '0':  writeFittingResLinear(xla, s, r, Attributi, res)
    elif    opt_dict['interp'] == '1':  writeFittingResAVD(xla, s, r, Attributi, res)
    else:                               writeFittingResSVE(xla, s, r, Attributi, res)


def writeFittingPyResOnXls(crv, xla, opt_dict, res, pos):

    nameSheet = nameSheetCurve
    try:
        s = xla.ActiveWorkbook.Sheets(nameSheet)
    except:
        s = xla.ActiveWorkbook.Sheets.Add()
        s.Name = nameSheet
    s.Activate()

    r = findCurveFromPos(xla, pos, nameSheet)
    print "------------------indirizzo trovato:", r.Value, r.Address
    #---
    #mi posiziono nella prima cella utile per scrivere i risultati del fitting
    #---
    Attributi = \
        {     "Date Ref"     : crv.ref_date
            , "Description"  : crv.description
            , "Currency"     : crv.curr
            , "Download Type": crv.download_type
            , "Quotation"    : crv.quotation
            , "CurveType"    : crv.type
            , "Interp. Model": fitt_translate[opt_dict['interp']]
        }

    row = r.Row
    col = r.Column
    r = xla.Range(xla.Cells(row+13, col ),(xla.Cells(row+13, col)))

    while (r.Value!= None):
        row = r.Row
        col = r.Column
        r = xla.Range(xla.Cells(row+1, col), xla.Cells(row+1, col ))
        if (r.Value == None):
            row = r.Row
            col = r.Column
            r = xla.Range(xla.Cells(row+1, col), xla.Cells(row+1,col))
            if (r.Value == None):
                row = r.Row
                col = r.Column
                r = xla.Range(xla.Cells(row + 1, col), xla.Cells(row + 1, col))

    print "------------------posizione calcolata per scrivere output fitting.:", r.Value, r.Address
    print opt_dict
    print opt_dict['interp'], type(opt_dict['interp'])
    if      opt_dict['interp'] == '0':  writeFittingResLinear(xla, s, r, Attributi, res)
    elif    opt_dict['interp'] == '1':  writeFittingResAVD(xla, s, r, Attributi, res)
    else:                               writeFittingResSVE(xla, s, r, Attributi, res)

