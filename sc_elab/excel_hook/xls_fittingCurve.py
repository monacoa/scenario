from pyxll import xlcAlert
import sys
import datetime
from win32com.client import constants as const
from Tkinter import *


from xls_utils import drawLine, drawBox
from DEF_intef import fitt_translate, FORMATT, nameSheetBootstrap, nameSheetCurve
from xls_swapCurve import findCurveFromPos
from xls_bootCurve import intestazioneSwapCurveSegmenti, findBootstrappedCurveFromPos



def writeFittingResLinear(xla, s, r, Attributi, res):

    # ---
    # calcolo il codice corretto per il metodo di IL
    # ---
    text = "PIL" + Attributi['Currency']+ Attributi['Return'] + "BLM" + str(Attributi['Date Ref'])[8:10] + str(Attributi['Date Ref'])[5:7] + str(Attributi['Date Ref'])[2:4]+ "_" +Attributi['Segms']


    dateRef = Attributi['Date Ref']
    ra = intestazioneSwapCurveSegmenti(xla, "", r, Attributi, nCols=2, text = text)
    r = s.Range(ra)
    topLeftRow = r.Row
    topLeftCol = r.Column
    xla.Cells(topLeftRow - 1, topLeftCol).Value = "Date"
    xla.Cells(topLeftRow - 1, topLeftCol).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 1).Value = "a"
    xla.Cells(topLeftRow - 1, topLeftCol + 1).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 2).Value = "b"
    xla.Cells(topLeftRow - 1, topLeftCol + 2).HorizontalAlignment = const.xlCenter
    nRows = len(res['Dates']) + 2
    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol, topLeftRow + nRows - 2, topLeftCol + 2)
    drawLine(xla, topLeftRow - 1, topLeftCol, topLeftRow - 1, topLeftCol + 2, "o", const.xlThin)
    for i in range(nRows):
        xla.Cells(topLeftRow + i, topLeftCol).NumberFormat = FORMATT
        xla.Cells(topLeftRow + i, topLeftCol).HorizontalAlignment = const.xlCenter
        if (i > 0):
            date = res['Dates'][i-1]
            xla.Cells(topLeftRow + i, topLeftCol).Value = date
            a = res['a'][i-1]
            b = res['b'][i-1]
            xla.Cells(topLeftRow + i, topLeftCol + 1).Value = a
            xla.Cells(topLeftRow + i, topLeftCol + 1).NumberFormat = "0.0000"
            xla.Cells(topLeftRow + i, topLeftCol + 1).HorizontalAlignment = const.xlCenter
            xla.Cells(topLeftRow + i, topLeftCol + 2).Value = b
            xla.Cells(topLeftRow + i, topLeftCol + 2).NumberFormat = "0.0000"
            xla.Cells(topLeftRow + i, topLeftCol + 2).HorizontalAlignment = const.xlCenter
        else:
            xla.Cells(topLeftRow + i, topLeftCol).Value = dateRef
            


def writeFittingResAVD(xla, s, r, Attributi, res):
    # ---
    # calcolo il codice corretto per il metodo di IL
    # ---
    text = "PIA" + Attributi['Currency'] + Attributi['Return'] + "BLM" + str(Attributi['Date Ref'])[8:10] + str(Attributi['Date Ref'])[5:7] + str(Attributi['Date Ref'])[2:4] + "_" + Attributi['Segms']

    ra = intestazioneSwapCurveSegmenti(xla, "", r, Attributi, nCols=2, text = text )
    r = s.Range(ra)
    dateRef = Attributi['Date Ref']
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

    nRows = len(res['Dates']) + 2
    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol, topLeftRow + nRows - 2, topLeftCol + 5)
    drawLine(xla, topLeftRow - 1, topLeftCol, topLeftRow - 1, topLeftCol + 5, "o", const.xlThin)
    for i in range(nRows):
        #date = res['Dates'][i]
        #xla.Cells(topLeftRow + i, topLeftCol).Value = date
        xla.Cells(topLeftRow + i, topLeftCol).NumberFormat = FORMATT
        xla.Cells(topLeftRow + i, topLeftCol).HorizontalAlignment = const.xlCenter
        if (i > 0):
            date = res['Dates'][i - 1]
            a = res['a'][i - 1]
            b = res['b'][i - 1]
            c = res['c'][i - 1]
            d = res['d'][i - 1]
            e = res['e'][i - 1]
            xla.Cells(topLeftRow + i, topLeftCol).Value = date
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
        else:
            xla.Cells(topLeftRow + i, topLeftCol).Value = dateRef


def writeFittingResSVE(xla, s, r, Attributi, res):
    text = "PIS" + Attributi['Currency'] + Attributi['Return'] +"BLM" + str(Attributi['Date Ref'])[8:10] + str(Attributi['Date Ref'])[5:7] + str( Attributi['Date Ref'])[2:4] + "_" + Attributi['Segms']
    ra = intestazioneSwapCurveSegmenti(xla, "", r, Attributi, nCols=2, text = text)
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

def writeFittingResNS(xla, s, r, Attributi, res):
    text = "PIN" + Attributi['Currency'] + Attributi['Return'] +"BLM" + str(Attributi['Date Ref'])[8:10] + str(Attributi['Date Ref'])[5:7] + str( Attributi['Date Ref'])[2:4] + "_" + Attributi['Segms']
    ra = intestazioneSwapCurveSegmenti(xla, "", r, Attributi, nCols=2, text = text)
    r = s.Range(ra)

    topLeftRow = r.Row
    topLeftCol = r.Column
    xla.Cells(topLeftRow - 1, topLeftCol).Value = "const1"
    xla.Cells(topLeftRow - 1, topLeftCol).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 1).Value = "beta 0"
    xla.Cells(topLeftRow - 1, topLeftCol + 1).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 2).Value = "beta 1"
    xla.Cells(topLeftRow - 1, topLeftCol + 2).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 3).Value = "beta 2"
    xla.Cells(topLeftRow - 1, topLeftCol + 3).HorizontalAlignment = const.xlCenter

    nRows = 1
    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol, topLeftRow + nRows - 1, topLeftCol + 3)
    drawLine(xla, topLeftRow - 1, topLeftCol, topLeftRow - 1, topLeftCol + 3, "o", const.xlThin)

    const1 = res['const1']
    beta0  = res['beta0']
    beta1  = res['beta1']
    beta2  = res['beta2']

    xla.Cells(topLeftRow , topLeftCol ).Value = const1
    xla.Cells(topLeftRow , topLeftCol ).NumberFormat = "0.0000"
    xla.Cells(topLeftRow , topLeftCol ).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow, topLeftCol+1).Value = beta0
    xla.Cells(topLeftRow, topLeftCol+1).NumberFormat = "0.0000"
    xla.Cells(topLeftRow, topLeftCol+1).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow, topLeftCol+2).Value = beta1
    xla.Cells(topLeftRow, topLeftCol+2).NumberFormat = "0.0000"
    xla.Cells(topLeftRow, topLeftCol+2).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow, topLeftCol+3).Value = beta2
    xla.Cells(topLeftRow, topLeftCol+3).NumberFormat = "0.0000"
    xla.Cells(topLeftRow, topLeftCol+3).HorizontalAlignment = const.xlCenter


def writeFittingResCIR(xla, s, r, Attributi, res):
    text = "PIC" + Attributi['Currency'] + Attributi['Return'] +"BLM" + str(Attributi['Date Ref'])[8:10] + str(Attributi['Date Ref'])[5:7] + str( Attributi['Date Ref'])[2:4] + "_" + Attributi['Segms']
    ra = intestazioneSwapCurveSegmenti(xla, "", r, Attributi, nCols=2, text = text)
    r = s.Range(ra)

    topLeftRow = r.Row
    topLeftCol = r.Column
    xla.Cells(topLeftRow - 1, topLeftCol).Value = "r0"
    xla.Cells(topLeftRow - 1, topLeftCol).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 1).Value = "kappa"
    xla.Cells(topLeftRow - 1, topLeftCol + 1).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 2).Value = "theta"
    xla.Cells(topLeftRow - 1, topLeftCol + 2).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 3).Value = "sigma"
    xla.Cells(topLeftRow - 1, topLeftCol + 3).HorizontalAlignment = const.xlCenter

    nRows = 1
    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol, topLeftRow + nRows - 1, topLeftCol + 3)
    drawLine(xla, topLeftRow - 1, topLeftCol, topLeftRow - 1, topLeftCol + 3, "o", const.xlThin)

    r0     = res['r0']
    kappa  = res['kappa']
    theta  = res['theta']
    sigma  = res['sigma']

    xla.Cells(topLeftRow , topLeftCol ).Value = r0
    xla.Cells(topLeftRow , topLeftCol ).NumberFormat = "0.0000"
    xla.Cells(topLeftRow , topLeftCol ).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow, topLeftCol+1).Value = kappa
    xla.Cells(topLeftRow, topLeftCol+1).NumberFormat = "0.0000"
    xla.Cells(topLeftRow, topLeftCol+1).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow, topLeftCol+2).Value = theta
    xla.Cells(topLeftRow, topLeftCol+2).NumberFormat = "0.0000"
    xla.Cells(topLeftRow, topLeftCol+2).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow, topLeftCol+3).Value = sigma
    xla.Cells(topLeftRow, topLeftCol+3).NumberFormat = "0.0000"
    xla.Cells(topLeftRow, topLeftCol+3).HorizontalAlignment = const.xlCenter





def writeFittingBootResOnXls(crv, xla, opt_dict, res, pos):
    nameSheet = nameSheetBootstrap
    try:
        s = xla.ActiveWorkbook.Sheets(nameSheet)
    except:
        s = xla.ActiveWorkbook.Sheets.Add()
        s.Name = nameSheet
    s.Activate()

    r = findBootstrappedCurveFromPos(xla, nameSheet,pos)

    #---
    #mi posiziono nella prima cella utile per scrivere i risultati del fitting
    #---

    Attributi = \
        {     "Date Ref"     : crv.ref_date
            , "Segms"        : crv.code
            , "Source"       : crv.source
            , "Description"  : crv.description
            , "Currency"     : crv.curr
            , "Download Type": crv.download_type
            , "Quotation"    : crv.quotation
            , "CurveType"    : crv.type
            , "Return"       : "ZC"
            , "Node Type"    : "nd"
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

    row = r.Row
    col = r.Column
    r = xla.Range(xla.Cells(row -( len(Attributi.keys()) +1), col), (xla.Cells(row -(len(Attributi.keys()) + 1), col)))

    if      opt_dict['interp'] == '0':  writeFittingResLinear(xla, s, r, Attributi, res)
    elif    opt_dict['interp'] == '1':  writeFittingResAVD(xla, s, r, Attributi, res)
    elif    opt_dict['interp'] == '2':  writeFittingResSVE(xla, s, r, Attributi, res)
    elif    opt_dict['interp'] == '3':  writeFittingResCIR(xla, s, r, Attributi, res)
    elif    opt_dict['interp'] == '4':   writeFittingResNS(xla, s, r, Attributi, res)

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
        {"Date Ref": crv.ref_date
            , "Segms": "ND"
            , "Source": crv.source
            , "Description": crv.description
            , "Currency": crv.curr
            , "Download Type": crv.download_type
            , "Quotation": crv.quotation
            , "CurveType": crv.type
            , "Return": "PY"
            , "Node Type": "nd"
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
    elif    opt_dict['interp'] == '2':  writeFittingResSVE(xla, s, r, Attributi, res)
    elif    opt_dict['interp'] == '3':  writeFittingResCIR(xla, s, r, Attributi, res)
    elif    opt_dict['interp'] == '4':  writeFittingResNS(xla, s, r, Attributi, res)

    else:                               writeFittingResSVE(xla, s, r, Attributi, res)




def readIntestazioneFitting(xla , r , cc):

    row              = r.Row
    col              = r.Column

    cc.code_interp   = xla.Range(xla.Cells(row , col), xla.Cells(row , col)).Value
    cc.curr          = xla.Range(xla.Cells(row + 1, col + 1), xla.Cells(row + 1, col + 1)).Value
    cc.type          = xla.Range(xla.Cells(row + 2, col + 1), xla.Cells(row + 2, col + 1)).Value

    dd = xla.Range(xla.Cells(row + 3, col + 1), xla.Cells(row + 3, col + 1)).Value
    cc.ref_date      = datetime.date(year = dd.year, month = dd.month, day = dd.day)

    cc.description   = xla.Range(xla.Cells(row + 4, col+1), xla.Cells(row+4, col+1)).Value
    cc.download_type = xla.Range(xla.Cells(row + 5, col + 1), xla.Cells(row + 5, col + 1)).Value
    cc.interpolation_model = xla.Range(xla.Cells(row + 6, col + 1), xla.Cells(row + 6, col + 1)).Value

    cc.quotation     = xla.Range(xla.Cells(row + 8, col + 1), xla.Cells(row + 8, col + 1)).Value
    cc.rendimento    = xla.Range(xla.Cells(row + 9, col + 1), xla.Cells(row + 9, col + 1)).Value
    cc.code          = xla.Range(xla.Cells(row + 10, col + 1), xla.Cells(row + 10, col + 1)).Value
    cc.source        = xla.Range(xla.Cells(row + 11, col + 1), xla.Cells(row + 11, col + 1)).Value

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


    return  xla.Range(xla.Cells(row + 12, col), xla.Cells(row + 12, col))


def readParmsNames(xla , r , cc):
    code = (cc.getCurveCode())[2]
    if ((code.upper() == 'S') or (code.upper() == 'C') or (code.upper() == 'N')) : offset = 0
    else                                                : offset = 1
    rp = xla.Range(xla.Cells(r.Row, r.Column+offset), xla.Cells(r.Row, r.Column+offset))
    print rp.Value, rp.Address
    while (rp.Value != None):
        par = rp.Value
        print "nome parametro:", par
        cc.parms_list.append(par)
        cc.interp_parms[par]= []
        rp = xla.Range(xla.Cells(rp.Row , rp.Column+1), xla.Cells(rp.Row , rp.Column+1))
    return

def readParms(xla, r, cc):
    code = (cc.getCurveCode())[2]

    rp = xla.Range(xla.Cells(r.Row+1, r.Column), xla.Cells(r.Row+1, r.Column))

    while rp.Value != None:
        i = 0
        print "CODE:", code
        #se non sono nel caso sve o cir, devo leggere le date
        if ((code.upper()) != 'S') and ((code.upper()) != 'C') and ((code.upper()) != 'N'):
            dd = rp.Value
            cc.interp_dates.append(datetime.date(year = dd.year, month = dd.month, day = dd.day))
            i = 1

        # ---
        for p in cc.parms_list:
            value = xla.Range(xla.Cells(rp.Row, rp.Column+i), xla.Cells(rp.Row , rp.Column+i)).Value
            if value == None: value = 0.0
            cc.interp_parms[p].append(value)
            i += 1
        # ---
        rp = xla.Range(xla.Cells(rp.Row + 1, rp.Column), xla.Cells(rp.Row + 1, rp.Column))

    cc.show()

def readInterpParmsFromXls (xla, range_curve, pos_parms):
    # ---
    # i parametri si trovano sempre alla colonna I
    # ---
    rp = xla.Range(xla.Cells(range_curve.Row, "I"), xla.Cells(range_curve.Row, "I"))
    p_i = 0
    nome_parms = ""
    while (rp.Value != None):
        p_i += 1
        nome_parms = rp.Value
        print "nome parametri:", nome_parms, "posizione:", p_i, "posizione che cerco:", pos_parms
        if p_i == pos_parms:  break

        if nome_parms[0:3] == "PIL":  distanzaO = 4
        elif nome_parms[0:3] == "PIN":  distanzaO = 7
        else                       :  distanzaO = 7

        rp = xla.Range(xla.Cells(rp.Row, rp.Column + distanzaO), xla.Cells(rp.Row, rp.Column+ distanzaO))

    from sc_elab.core.SwpCurve import InterpolatedCurve

    iCurve = InterpolatedCurve()
    rp = readIntestazioneFitting(xla, rp, iCurve)

    readParmsNames(xla, rp, iCurve)
    readParms(xla, rp, iCurve)

    return iCurve


import tkMessageBox

def readInterpolatedCurveFromXls (xla, sheet, pos_curve,  pos_parms):
    try:
        book = xla.ActiveWorkbook
        s = book.Sheets(sheet)
        s.Activate()
    except:

        msg = "Missing input sheet for Swap Curves in your workbook... \nNothing to do for me!"
        tkMessageBox.showinfo("Warning!", msg)

        return False

    rangeStart = "B2"

    r = findBootstrappedCurveFromPos(xla, sheet, pos_curve)
    print "ho trovato la curva boot da posizione e star in :", r.Address
    iCurve = readInterpParmsFromXls(xla, r, pos_parms)
    return iCurve