from pyxll import xlcAlert
import  sys
import  datetime
from    win32com.client import constants as const
from Tkinter import *
from sc_elab.core.SwpCurve import dict_segm2, Segm, Curve, CdsCurve
from map_bf_fields import field_map

from xls_utils import drawBox, drawLine, formatTestataCurva, findRigthPlaceBootCurveSeg



def intestazioneBondPortfolio( xla, sheet, rng,  attributi, nCols = 2, text= None):

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



def bondPortfolioAnag(xla, rangeS, crv):
    from DEF_intef import FORMATT
    rangeStart = rangeS
    topLeftRow = xla.Range(rangeStart).Row
    topLeftCol = xla.Range(rangeStart).Column
    nbonds = len(crv.bond_portfolio)
    nRows = nbonds
    nCols = 31

    #box intestazione
    drawBox(xla, const.xlMedium, topLeftRow, topLeftCol, topLeftRow + nRows+1, topLeftCol + nCols - 1)

    formatTestataCurva(xla, topLeftRow, topLeftCol, nCols, "Bond data")


    #Linea orizzontale di separazione
    drawLine(xla, topLeftRow + 1, topLeftCol, topLeftRow + 1, topLeftCol + nCols - 1, "o", const.xlThin)

    field_list = ['Isin', 'Descrizione', 'Seniority', 'Tipo tasso', 'Tipo quotazione', 'Tipo rimborso', 'Prezzo emissione', 'Prezzo rimborso / Inflation Ratio', 'Data emissione',  'Data scadenza', 
        'Tempo scadenza', 'Giorni di fixing',    'Tipo fixing',    'Basis',    'Adjustment',    'Periodicita cedola (mesi)',    'Tenor del tasso floater (anni)',    
        'Tasso cedolare annuo (Fisso/spread)',    'Cedola in corso',    'Prezzo-MID',    'Prezzo-BID',    'Prezzo-ASK',    'YTM/DM (MID)',    'YTM/DM (BID)',    'YTM/DM (ASK)',
        'Tasso di riferiemnto',    'Tasso repo',    'Data prezzo di mercato',    'Contributor',    'Peso',    'Indicizzazione']

    n_field = len(field_list)
    for i in range(0, n_field):
        xla.Cells(topLeftRow + 1, topLeftCol + i).Value = field_list[i]
        xla.Cells(topLeftRow + 1, topLeftCol + i).HorizontalAlignment = const.xlCenter


    n_bonds = len(crv.bond_portfolio)
    #field_list_n = crv.bond_portfolio[0].keys()

    for i in range(0, n_bonds):
        for j in range(0, n_field):
            
            a = xla.Cells(topLeftRow  + i + 2, topLeftCol+j)
            
            
            fieldNameTmp = field_list[j]
            fieldNameDBTmp = field_map[fieldNameTmp][0]
            
            
            
            if (fieldNameDBTmp != None):
                fieldValTmp  = crv.bond_portfolio[i][fieldNameDBTmp]
                
                if (fieldValTmp == None):
                    fieldValTmp  = field_map[fieldNameTmp][1]
                
                if (fieldNameDBTmp == 'FrequenzaCoupon'):
                    fieldValTmp = int(1.0/float(fieldValTmp)*12.0)

                cTypeTmp  = crv.bond_portfolio[i]['CouponType']
                if (fieldNameTmp == 'Tasso cedolare annuo (Fisso/spread)' and cTypeTmp == 'FIXED'):
                    fieldValTmp = crv.bond_portfolio[i]['CedolaInCorso']
                
            else:
                fieldValTmp  = field_map[fieldNameTmp][1]
            
            #if (i == 0):
            #    print 'fieldNameTmp: ', fieldNameTmp
            #    print 'fieldValTmp: ', fieldValTmp
            #    print '-----------------------------------------'
            

            a.Value = fieldValTmp
            if (type(fieldValTmp) == datetime.date) or (type(fieldValTmp) == datetime.datetime): a.NumberFormat = FORMATT
            a.HorizontalAlignment = const.xlCenter
            

    """
    f = 1
    i = 2
    for tag, mat, value in zip(crv.tags, crv.mats, crv.values):
        ll = [tag, mat, value]
        for j in range(3):
            a = xla.Cells(topLeftRow  + i, topLeftCol+j)
            a.Value = ll[j]
            if (type(ll[j]) == datetime.date) or (type(ll[j]) == datetime.datetime): a.NumberFormat = FORMATT
            if (type(ll[j]) == float)         : a.NumberFormat = "0.00"
            a.HorizontalAlignment = const.xlCenter
            j +=1
        i+=1
    """
    
    rangeStart = xla.Cells(topLeftRow + nRows + 2, topLeftCol).Address
    return rangeStart
#=======================================================================================================================
def writePortfoliOnXls(crv, nameSheet, xla, curve_type):
    # ---
    from DEF_intef import FORMATT
    # ---
    rangeStart = "B2"
    distCurve  = 5
    # ---
    #Individuo posizione in cui scrivere
    # ---
    sheet = xla.ActiveWorkbook.Sheets(nameSheet)
    r     = sheet.Range(rangeStart)
    #rOut  =  findRigthPlaceBootCurveSeg(xla, r, distCurve)
    rOut  = findRigthPlaceBootCurveSeg(xla, r, distCurve, "v")


    # ---
    #Genero il blocco di intestazione della curva
    # ---

    
    
    
    Attributi =   {"Date Ref": crv.ref_date
                      , "Description": crv.description
                      , "Currency": crv.curr
                      , "Issuer": crv.emittente
                      , "Sector": crv.settore
                      , "Recovery Rate (%)": crv.recoveryRate
                      , "Rating": crv.rating
                      , "Source": crv.source
                      }


    rangeStartNew = intestazioneBondPortfolio ( xla, sheet , rOut, Attributi)
    rangeStartNew = bondPortfolioAnag(xla, rangeStartNew, crv)


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

def readIntestazioneCds(xla , r , cc):
    row             = r.Row
    col             = r.Column

    cc.capitalization = xla.Range(xla.Cells(row+1, col+1), xla.Cells(row+1, col+1)).Value
    cc.curr           = xla.Range(xla.Cells(row+2, col+1), xla.Cells(row+2, col+1)).Value
    # imposto il mercato
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

    print ".................", r.Value
    cc.type = xla.Range(xla.Cells(row + 3, col + 1), xla.Cells(row + 3, col + 1)).Value
    dd = xla.Range(xla.Cells(row + 4, col + 1), xla.Cells(row + 4, col + 1)).Value
    cc.ref_date = datetime.date(year=dd.year, month=dd.month, day=dd.day)

    cc.dayAdj   = xla.Range(xla.Cells(row + 5, col + 1), xla.Cells(row + 5, col + 1)).Value
    cc.dayCount = xla.Range(xla.Cells(row + 6, col + 1), xla.Cells(row + 6, col + 1)).Value

    cc.description   = xla.Range(xla.Cells(row + 7, col + 1), xla.Cells(row + 7, col + 1)).Value
    cc.download_type = xla.Range(xla.Cells(row + 8, col + 1), xla.Cells(row + 8, col + 1)).Value
    cc.lag           = xla.Range(xla.Cells(row + 9, col + 1), xla.Cells(row + 9, col + 1)).Value
    cc.frequency     = xla.Range(xla.Cells(row + 10, col + 1), xla.Cells(row + 10, col + 1)).Value
    cc.emittente     = xla.Range(xla.Cells(row + 11, col + 1), xla.Cells(row + 11, col + 1)).Value
    cc.node_type     = xla.Range(xla.Cells(row + 12, col + 1), xla.Cells(row + 12, col + 1)).Value
    cc.quotation     = xla.Range(xla.Cells(row + 13, col + 1), xla.Cells(row + 13, col + 1)).Value
    cc.rating        = xla.Range(xla.Cells(row + 14, col + 1), xla.Cells(row + 14, col + 1)).Value
    cc.recovery      = xla.Range(xla.Cells(row + 15, col + 1), xla.Cells(row + 15, col + 1)).Value
    cc.rendimento    = xla.Range(xla.Cells(row + 16, col + 1), xla.Cells(row + 16, col + 1)).Value
    cc.sector        = xla.Range(xla.Cells(row + 17, col + 1), xla.Cells(row + 17, col + 1)).Value
    cc.seniority     = xla.Range(xla.Cells(row + 18, col + 1), xla.Cells(row + 18, col + 1)).Value
    cc.source        = xla.Range(xla.Cells(row + 19, col + 1), xla.Cells(row + 19, col + 1)).Value

    r = xla.Range(xla.Cells(row + 21, col), xla.Cells(row + 21, col))
    cc.show()
    return r




def readSegm(xla, r, cc):
    row  = r.Row
    col  = r.Column
    name = r.Value
    print "SONO QUIIIIIIIIIII", name
    i = 2

    while r.Value != None:
        r = xla.Range(xla.Cells(row + 2, col), xla.Cells(row + 2, col))
        tag   = r.Value
        print "tag", tag
        time  = xla.Range(xla.Cells(row + i, col + 1), xla.Cells(row + i, col + 1)).Value
        value = xla.Range(xla.Cells(row + i, col + 2), xla.Cells(row + i, col + 2)).Value

        cc.tags.append(tag)
        cc.mats.append(time)
        cc.values.append(value)
        i += 1
        r = xla.Range(xla.Cells(row + i, col), xla.Cells(row + i, col))

    cc.show()
    return r


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



