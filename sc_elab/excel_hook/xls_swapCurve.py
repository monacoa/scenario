from pyxll import xlcAlert
import  sys
import  datetime
from    win32com.client import constants as const
from Tkinter import *
from sc_elab.core.SwpCurve import dict_segm2, Segm, Curve

from xls_utils import drawBox, drawLine, formatTestataCurva, findRigthPlaceBootCurveSeg



def intestazioneSwapCurveSegmenti( xla, sheet, rng,  attributi, nCols = 2, text= None):
    print attributi
    txt = text if text!=None else attributi['Description']
    nRows           = len(attributi.keys())
    topLeftRow      = rng.Row
    topLeftCol      = rng.Column
    drawBox            (xla, 3,topLeftRow, topLeftCol,topLeftRow + nRows, topLeftCol + nCols - 1, 0)
    formatTestataCurva (xla, topLeftRow, topLeftCol, nCols, txt)

    kk = (attributi.keys())
    kk.sort()
    i = 0
    for k in kk:
        xla.Cells(topLeftRow + 1+ i, topLeftCol).Value   = k
        xla.Cells(topLeftRow + 1+ i, topLeftCol+1).Value = attributi[k]
        if (type(attributi[k]) == datetime.datetime) or (type(attributi[k]) == datetime.date):
            print "SONO QUI!"
            print "i:", i
            print "k", k, "attr[k]", attributi[k], type(attributi[k])
            print "topLeftRow", topLeftRow, "tlc", topLeftCol
            print "*"*120
            xla.Cells(topLeftRow + 1 + i, topLeftCol + 1).NumberFormat = "gg/MM/aaaa"
        print "sono uscita da if e setto center"
        xla.Cells(topLeftRow + 1 + i, topLeftCol + 1).HorizontalAlignment = const.xlCenter
        print "center fatto!"
        i+=1

    rangeStart = xla.Cells(topLeftRow + nRows + 2, topLeftCol).Address
    return rangeStart


def displayHWParamSwCurve (xla, rangeStart, attributi):
    nRows = 1
    nCols = 4
    # 2 parametri + descrizione
    topLeftRow = xla.Range(rangeStart).Row
    topLeftCol = xla.Range(rangeStart).Column
    nomeBox = "Par. Hull & White per Conv. Adj."

    # box parametri
    drawBox(xla, const.xlMedium, topLeftRow, topLeftCol, topLeftRow + nRows, topLeftCol + nCols - 1, 0)
    drawLine(xla, topLeftRow + 1,topLeftCol + 1, topLeftRow + nRows, topLeftCol + 1, "v", const.xlThin)

    # ' intestazione
    formatTestataCurva(xla, topLeftRow, topLeftCol, nCols, nomeBox)
    # Scirttura descrizione e parametri

    xla.Cells(topLeftRow + 1, topLeftCol).Value = "mean revertion speed"
    xla.Cells(topLeftRow + 1, topLeftCol+1).Value = "0.0"
    xla.Cells(topLeftRow + 1, topLeftCol+2).Value = "volatility"
    xla.Cells(topLeftRow + 1, topLeftCol+3).Value = "0.0"

    rangeStartN = xla.Cells(topLeftRow + 3, topLeftCol).Address
    return rangeStartN

def segmentoSwapCurve(xla, rangeS, code, segm):
    from DEF_intef import FORMATT
    rangeStart = rangeS
    topLeftRow = xla.Range(rangeStart).Row
    topLeftCol = xla.Range(rangeStart).Column
    nNodi = len(segm.dates)
    nRows = nNodi + 1
    nCols = 4

    #box intestazione
    drawBox(xla, const.xlMedium, topLeftRow, topLeftCol, topLeftRow + nRows, topLeftCol + nCols - 1)
    isFut = False
    #Formatto zona codice curva
    if   code == "G": nomeSegmento = "0. Short term swap"
    elif code == "D": nomeSegmento = "1. Depositi"
    elif code == "L": nomeSegmento = "2. Libor"
    elif code =="F" :
        nomeSegmento = "3. Futures"
        isFut = True
    elif code == "S": nomeSegmento = "4. Swap Rate"
    else:             nomeSegmento = code

    formatTestataCurva(xla, topLeftRow, topLeftCol, nCols, nomeSegmento)
    #Linea orizzontale di separazione
    drawLine(xla, topLeftRow + 1, topLeftCol, topLeftRow + 1, topLeftCol + nCols - 1, "o", const.xlThin)
    drawLine(xla, topLeftRow + 1, topLeftCol+nCols - 2, topLeftRow + nRows, topLeftCol + nCols - 2, "v", const.xlThin)

    xla.Cells(topLeftRow + 1, topLeftCol + 0).Value = "Node"
    xla.Cells(topLeftRow + 1, topLeftCol + 0).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow + 1, topLeftCol + 1).Value = "Maturity"
    xla.Cells(topLeftRow + 1, topLeftCol + 1).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow + 1, topLeftCol + 2).Value = "Value"
    xla.Cells(topLeftRow + 1, topLeftCol + 2).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow + 1, topLeftCol + 3).Value = "Usage"
    xla.Cells(topLeftRow + 1, topLeftCol + 3).HorizontalAlignment = const.xlCenter

    f = 1
    i = 1
    for tag,date,value in zip(segm.tags, segm.dates, segm.values):
        ll = [tag, date, value]
        for j in range(3):
            a = xla.Cells(topLeftRow + 1 + i, topLeftCol+j)
            if (j == 0) and (isFut):
                a.Value = f
                a.NumberFormat = "0"
                f += 1
            else :
                a.Value = ll[j]
                if (type(ll[j]) == datetime.date) or (type(ll[j]) == datetime.datetime): a.NumberFormat = FORMATT
                if (type(ll[j]) == float)         : a.NumberFormat = "0.00"
            a.HorizontalAlignment = const.xlCenter
            j +=1

        b = xla.Cells(topLeftRow + 1+ i, topLeftCol+3)
        b.Value = "Y"
        b.HorizontalAlignment = const.xlCenter
        i+=1

    rangeStart = xla.Cells(topLeftRow + nRows + 2, topLeftCol).Address
    return rangeStart


#=======================================================================================================================
def writeCurveOnXls(crv, nameSheet, xla):
    from DEF_intef import FORMATT
    rangeStart = "B2"
    distCurve  = 5
    #Individuo posizione in cui scrivere
    sheet = xla.ActiveWorkbook.Sheets(nameSheet)
    r = sheet.Range(rangeStart)
    rOut =  findRigthPlaceBootCurveSeg(xla, r, distCurve)
    #rangeStartNew = rOut.Address

    #Genero il blocco di intestazione della curva

    Attributi     =  { "Date Ref": crv.ref_date
                     , "Description":crv.description
                     , "Currency":crv.curr
                     , "Download Type": crv.download_type
                     , "Quotation": crv.quotation
                     , "Source": crv.source
                     , "Tenor Floater Rate":crv.floater_tenor}

    rangeStartNew = intestazioneSwapCurveSegmenti ( xla, sheet , rOut, Attributi)
    # Genero il blocco per i parametri di Hull e White
    rangeStartNew = displayHWParamSwCurve (xla, rangeStartNew, Attributi)
    # ' Genero i blocchi dei segmenti della curva
    # codes = "DLGFS"; D = dep, L = libor, F = futures, G=swap sotto 1y, S = swap >=1Y
    cd = "DLGFS"
    for j in  range (len(cd)):
        code =cd[j]
        for s in crv.segms.keys():
            if code == s[0]:
                # visualizzo segmento
                rangeStartNew = segmentoSwapCurve (xla, rangeStartNew, code, crv.segms[s])


def readIntestazione(xla , r , cc):
    row              = r.Row
    col              = r.Column
    cc.description   = xla.Range(xla.Cells(row , col   ), xla.Cells(row , col   )).Value
    cc.curr          = xla.Range(xla.Cells(row + 1, col + 1), xla.Cells(row + 1, col + 1)).Value
    dd = xla.Range(xla.Cells(row + 2, col + 1), xla.Cells(row + 2, col + 1)).Value

    cc.ref_date      = datetime.date(year = dd.year, month = dd.month, day = dd.day)
    cc.download_type = xla.Range(xla.Cells(row + 4, col + 1), xla.Cells(row + 4, col + 1)).Value
    cc.quotation     = xla.Range(xla.Cells(row + 5, col + 1), xla.Cells(row + 5, col + 1)).Value
    cc.source        = xla.Range(xla.Cells(row + 6, col + 1), xla.Cells(row + 6, col + 1)).Value
    cc.floater_tenor = xla.Range(xla.Cells(row + 7, col + 1), xla.Cells(row + 7, col + 1)).Value
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
    r =  xla.Range(xla.Cells(row + 9, col), xla.Cells(row + 9, col))
    cc.show()
    return r


def readParametriHW(xla,r,cc):
    row = r.Row
    col = r.Column
    cc.HWparms ['meanRS']= xla.Range(xla.Cells(row+1, col+1), xla.Cells(row+1, col+1)).Value
    cc.HWparms ['sigma'] = xla.Range(xla.Cells(row + 1, col + 3), xla.Cells(row + 1, col + 3)).Value
    r = xla.Range(xla.Cells(row + 3, col ), xla.Cells(row + 3, col))
    cc.show()
    return r

def readSegms(xla, r, cc):

    while (r.Value != None):
        row  = r.Row
        col  = r.Column
        name = r.Value
        if   name == "0. Short term swap": code = "G"
        elif name == "1. Depositi"       : code = "D"
        elif name == "2. Libor"          : code = "L"
        elif name == "3. Futures"        : code = "F"
        elif name == "4. Swap Rate"      : code = "S"
        else                             : code = name

        if code == 'G':
            add = cc.floater_tenor[0]
            print add
            code += add

        cc.segms [dict_segm2[code]] = Segm()
        ss = cc.segms [dict_segm2[code]]
        i = 2
        r = xla.Range(xla.Cells(row + 2, col), xla.Cells(row + 2, col))
        while (r.Value != None):

            tag   = r.Value
            date  = xla.Range(xla.Cells(row + i, col + 1), xla.Cells(row + i, col + 1)).Value
            value = xla.Range(xla.Cells(row + i, col + 2), xla.Cells(row + i, col + 2)).Value
            use   = xla.Range(xla.Cells(row + i, col + 3), xla.Cells(row + i, col + 3)).Value

            ss.tags.append(tag)
            ss.dates.append(datetime.date(year = date.year, month = date.month, day = date.day))
            ss.values.append(value)
            ss.usage.append(use)
            i += 1
            r = xla.Range(xla.Cells(row + i, col), xla.Cells(row + i, col))
            #---

        r = xla.Range(xla.Cells(row + i + 1, col), xla.Cells(row + i + 1, col))



    cc.fillAnagSegm()
    cc.show()
    return r


def readCurveFromXls(xla, des, pos, nameSheet):
    rangeStart = "B2"
    distCurve = 5
    sheet = xla.ActiveWorkbook.Sheets(nameSheet)
    r     = sheet.Range(rangeStart)

    nCols = r.Columns.Count
    row = r.Row
    col = r.Column

    r = xla.Range(xla.Cells(row, col + distCurve*(pos-1)), xla.Cells(row, col + (nCols-1) + distCurve*(pos-1)))

    cc = Curve()
    r = readIntestazione(xla, r,cc)
    r = readParametriHW(xla,r,cc)
    r = readSegms(xla,r,cc)

    return cc


def findCurveFromPos(xla, pos, nameSheet):
    rangeStart = "B2"
    distCurve = 5
    print "pos;", pos, type(pos)
    sheet = xla.ActiveWorkbook.Sheets(nameSheet)
    r     = sheet.Range(rangeStart)
    nCols = r.Columns.Count
    row   = r.Row
    col   = r.Column
    r     = xla.Range(xla.Cells(row, col + distCurve*(pos - 1)), xla.Cells(row, col + (nCols-1) + distCurve*(pos-1)))
    return r
