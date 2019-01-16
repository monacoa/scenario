from pyxll import xlcAlert
import  sys
import  datetime
from    win32com.client import constants as const
from Tkinter import *
from sc_elab.core.SwpCurve import dict_segm2, Segm, Curve, CdsCurve
from sc_elab.core.BondPortfolio import BondPortfolio
from xls_swapCurve import intestazioneSwapCurveSegmenti

from DEF_intef import FORMATT, nameSheetBootstrap, nameSheetBondSurv, nameSheetBondSpread, nameSheetBondResults, nameSheetBondOptPrms


from map_bf_fields import field_map

from xls_utils import drawBox, drawLine, formatTestataCurva, findRigthPlaceBootCurveSeg



def writeIntestazionePortfolio( xla, sheet, rng,  attributi, nCols = 2, text= None):

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
            #print "SONO QUI!"
            #print "i:", i
            #print "k", k, "attr[k]", attributi[k], type(attributi[k])
            #print "topLeftRow", topLeftRow, "tlc", topLeftCol
            #print "*"*120
            xla.Cells(topLeftRow + 1 + i, topLeftCol + 1).NumberFormat = "gg/MM/aaaa"
        #print "sono uscita da if e setto center"
        xla.Cells(topLeftRow + 1 + i, topLeftCol + 1).HorizontalAlignment = const.xlCenter
        #print "center fatto!"
        i+=1

    rangeStart = xla.Cells(topLeftRow + nRows + 2, topLeftCol).Address
    return rangeStart



def writePortfolioAnag(xla, rangeS, crv):
    from DEF_intef import FORMATT
    rangeStart = rangeS
    topLeftRow = xla.Range(rangeStart).Row - 1
    topLeftCol = xla.Range(rangeStart).Column
    nbonds = len(crv.bond_portfolio)
    nRows = nbonds
    nCols = 31

    #box intestazione
    drawBox(xla, const.xlMedium, topLeftRow, topLeftCol, topLeftRow + nRows, topLeftCol + nCols - 1)

    #formatTestataCurva(xla, topLeftRow, topLeftCol, nCols, "Bond data")


    #Linea orizzontale di separazione
    drawLine(xla, topLeftRow, topLeftCol, topLeftRow, topLeftCol + nCols - 1, "o", const.xlThin)

    field_list = ['Isin', 'Descrizione', 'Seniority', 'Tipo tasso', 'Tipo quotazione', 'Tipo rimborso', 'Prezzo emissione', 'Prezzo rimborso / Inflation Ratio', 'Data emissione',  'Data scadenza', 
        'Tempo scadenza', 'Giorni di fixing',    'Tipo fixing',    'Basis',    'Adjustment',    'Periodicita cedola (mesi)',    'Tenor del tasso floater (mesi)',
        'Tasso cedolare annuo (Fisso/spread)',    'Cedola in corso',    'Prezzo-MID',    'Prezzo-BID',    'Prezzo-ASK',    'YTM/DM (MID)',    'YTM/DM (BID)',    'YTM/DM (ASK)',
        'Tasso di riferiemnto',    'Tasso repo',    'Data prezzo di mercato',    'Contributor',    'Peso',    'Indicizzazione']

    n_field = len(field_list)
    for i in range(0, n_field):
        xla.Cells(topLeftRow , topLeftCol + i).Value = field_list[i]
        xla.Cells(topLeftRow , topLeftCol + i).HorizontalAlignment = const.xlCenter


    n_bonds = len(crv.bond_portfolio)
    #field_list_n = crv.bond_portfolio[0].keys()

    for i in range(0, n_bonds):
        for j in range(0, n_field):
            
            
            a = xla.Cells(topLeftRow  + i + 1, topLeftCol+j)
            
            
            fieldNameTmp = field_list[j]
            fieldNameDBTmp = field_map[fieldNameTmp][0]
            
            #print 'fieldNameDBTmp: ', fieldNameDBTmp
            
            
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

def writeBondFittingRes4OnXls(crv, xla, str_boot_opt, res, codice_curva):


    nameSheet = nameSheetBondOptPrms
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


    Attributi_1 = \
        { "Data riferimento" : crv.ref_date
        , "Descrizione"   : crv.description
        , "Valuta"        : crv.curr
        , "Tipo nodo"     : "HR model prms"
        , "Emittente"     : crv.emittente
        , "Recovery rate" : crv.recoveryRate
        , "Rating" : crv.rating
        , "Modello hazard rate" : crv.HRateModel
        , "Modello di valutazione" : crv.BondModel
        }

    

    ra = intestazioneSwapCurveSegmenti(xla, s, r , Attributi_1, nCols=2, text = codice_curva)

    r  = s.Range(ra)

    topLeftRow = r.Row
    topLeftCol = r.Column

    xla.Cells(topLeftRow - 1, topLeftCol).Value = "Name params"
    xla.Cells(topLeftRow - 1, topLeftCol).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow - 1, topLeftCol + 1).Value = "Value params"
    xla.Cells(topLeftRow - 1, topLeftCol + 1).HorizontalAlignment = const.xlCenter


    
    
    name_prms = res['dict_opt_prms'].keys()
    nRows     = len(name_prms)
    

    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol, topLeftRow + nRows - 1, topLeftCol + 1)
    drawLine(xla, topLeftRow - 1, topLeftCol, topLeftRow - 1, topLeftCol + 1, "o", const.xlThin)



    
    for i in range(nRows):
        


        name_prms_tmp    = name_prms[i]
        value_prms_tmp   = res['dict_opt_prms'][name_prms_tmp]


        xla.Cells(topLeftRow + i, topLeftCol).Value = name_prms_tmp
        #xla.Cells(topLeftRow + i, topLeftCol).NumberFormat = "0.00"
        xla.Cells(topLeftRow + i, topLeftCol).HorizontalAlignment = const.xlCenter

        xla.Cells(topLeftRow + i, topLeftCol + 1).Value = value_prms_tmp
        xla.Cells(topLeftRow + i, topLeftCol + 1).NumberFormat = "0.00"
        xla.Cells(topLeftRow + i, topLeftCol + 1).HorizontalAlignment = const.xlCenter
        

#=======================================================================================================================

def writeBondFittingRes3OnXls(crv, xla, str_boot_opt, res, codice_curva):


    nameSheet = nameSheetBondResults
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


    Attributi_1 = \
        { "Data riferimento" : crv.ref_date
        , "Descrizione"   : crv.description
        , "Valuta"        : crv.curr
        , "Tipo nodo"     : "TQ Price"
        , "Emittente"     : crv.emittente
        , "Recovery rate" : crv.recoveryRate
        , "Rating" : crv.rating
        , "Modello hazard rate" : crv.HRateModel
        , "Modello di valutazione" : crv.BondModel
        }

    
    #r2 = s.Range(xla.Cells(r.Row, r.Column + 3), xla.Cells(r.Row, r.Column + 3))
    #r3 = s.Range(xla.Cells(r.Row, r.Column + 6), xla.Cells(r.Row, r.Column + 7))

    ra = intestazioneSwapCurveSegmenti(xla, s, r , Attributi_1, nCols=3, text = codice_curva)
    #rb = intestazioneSwapCurveSegmenti(xla, s, r2, Attributi_2, nCols=2, text = codice_curva)
    #rc = intestazioneSwapCurveSegmenti(xla, s, r3, Attributi_3, nCols=2, text = codice_curva)

    r  = s.Range(ra)

    topLeftRow = r.Row
    topLeftCol = r.Column

    xla.Cells(topLeftRow - 1, topLeftCol).Value = "Maturity"
    xla.Cells(topLeftRow - 1, topLeftCol).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow - 1, topLeftCol + 1).Value = "MKT Price"
    xla.Cells(topLeftRow - 1, topLeftCol + 1).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow - 1, topLeftCol + 2).Value = "Model Price"
    xla.Cells(topLeftRow - 1, topLeftCol + 2).HorizontalAlignment = const.xlCenter

    

    nRows = len(res['bondTimes'])

    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol, topLeftRow + nRows - 1, topLeftCol + 2)
    drawLine(xla, topLeftRow - 1, topLeftCol, topLeftRow - 1, topLeftCol + 2, "o", const.xlThin)

    #drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol, topLeftRow + nRows - 1, topLeftCol + 2)
    #drawLine(xla, topLeftRow - 1, topLeftCol+3, topLeftRow - 1, topLeftCol + 2, "o", const.xlThin)


    
    for i in range(nRows):
        


        bond_times         = res['bondTimes'][i]
        opt_clean_prices   = res['opt_clean_prices'][i]
        mkt_clean_prices   = res['mkt_clean_prices'][i]


        xla.Cells(topLeftRow + i, topLeftCol).Value = bond_times
        xla.Cells(topLeftRow + i, topLeftCol).NumberFormat = "0.00"
        xla.Cells(topLeftRow + i, topLeftCol).HorizontalAlignment = const.xlCenter

        xla.Cells(topLeftRow + i, topLeftCol + 1).Value = mkt_clean_prices
        xla.Cells(topLeftRow + i, topLeftCol + 1).NumberFormat = "0.00"
        xla.Cells(topLeftRow + i, topLeftCol + 1).HorizontalAlignment = const.xlCenter
        
        xla.Cells(topLeftRow + i, topLeftCol + 2).Value = opt_clean_prices
        xla.Cells(topLeftRow + i, topLeftCol + 2).NumberFormat = "0.00"
        xla.Cells(topLeftRow + i, topLeftCol + 2).HorizontalAlignment = const.xlCenter


        



        
#=======================================================================================================================

def writeBondFittingRes2OnXls(crv, xla, str_boot_opt, res, codice_curva):


    nameSheet = nameSheetBondSpread
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


    Attributi_1 = \
        { "Data riferimento"        : crv.ref_date
        , "Descrizione"     : crv.description
        , "Valuta"        : crv.curr
        , "Tipo nodo"     : "ZSpread"
        , "Emittente"     : crv.emittente
        , "Recovery rate" : crv.recoveryRate
        , "Rating" : crv.rating
        , "Modello hazard rate" : crv.HRateModel
        , "Modello di valutazione" : crv.BondModel
        }

    Attributi_2 = \
        { "Data riferimento"        : crv.ref_date
        , "Descrizione"     : crv.description
        , "Valuta"        : crv.curr
        , "Tipo nodo"     : "PY Spread"
        , "Emittente"     : crv.emittente
        , "Recovery rate" : crv.recoveryRate
        , "Rating" : crv.rating
        , "Modello hazard rate" : crv.HRateModel
        , "Modello di valutazione" : crv.BondModel

        }

    Attributi_3 = \
        { "Data riferimento"        : crv.ref_date
        , "Descrizione"     : crv.description
        , "Valuta"        : crv.curr
        , "Tipo nodo"     : "PY Risk free"
        , "Emittente"     : crv.emittente
        , "Recovery rate" : crv.recoveryRate
        , "Rating" : crv.rating
        , "Modello risk free" : crv.HRateModel
        , "Modello di valutazione" : crv.BondModel

        }


    r2 = s.Range(xla.Cells(r.Row, r.Column + 3), xla.Cells(r.Row, r.Column + 3))
    r3 = s.Range(xla.Cells(r.Row, r.Column + 6), xla.Cells(r.Row, r.Column + 7))
    r4 = s.Range(xla.Cells(r.Row, r.Column + 9), xla.Cells(r.Row, r.Column + 11))

    ra = intestazioneSwapCurveSegmenti(xla, s, r , Attributi_1, nCols=2, text = codice_curva)
    rb = intestazioneSwapCurveSegmenti(xla, s, r2, Attributi_2, nCols=2, text = codice_curva)
    rc = intestazioneSwapCurveSegmenti(xla, s, r3, Attributi_3, nCols=2, text = codice_curva)

    r  = s.Range(ra)

    topLeftRow = r.Row
    topLeftCol = r.Column

    xla.Cells(topLeftRow - 1, topLeftCol).Value = "Date"
    xla.Cells(topLeftRow - 1, topLeftCol).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 1).Value = "Value (x 100)"
    xla.Cells(topLeftRow - 1, topLeftCol + 1).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow - 1, topLeftCol + 3).Value = "Date"
    xla.Cells(topLeftRow - 1, topLeftCol + 3).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 4).Value = "Value (x 100)"
    xla.Cells(topLeftRow - 1, topLeftCol + 4).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow - 1, topLeftCol + 6).Value = "Date"
    xla.Cells(topLeftRow - 1, topLeftCol + 6).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 7).Value = "Value (x 100)"
    xla.Cells(topLeftRow - 1, topLeftCol + 7).HorizontalAlignment = const.xlCenter

    
    #xla.Cells(topLeftRow - 1, topLeftCol + 6).Value = "Date"
    #xla.Cells(topLeftRow - 1, topLeftCol + 6).HorizontalAlignment = const.xlCenter
    #xla.Cells(topLeftRow - 1, topLeftCol + 7).Value = "Value"
    #xla.Cells(topLeftRow - 1, topLeftCol + 7).HorizontalAlignment = const.xlCenter

    

    nRows = len(res['outputDates'])

    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol, topLeftRow + nRows - 1, topLeftCol + 1)
    drawLine(xla, topLeftRow - 1, topLeftCol, topLeftRow - 1, topLeftCol + 1, "o", const.xlThin)

    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol + 3, topLeftRow + nRows - 1, topLeftCol + 4)
    drawLine(xla, topLeftRow - 1, topLeftCol+3, topLeftRow - 1, topLeftCol + 4, "o", const.xlThin)

    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol + 6, topLeftRow + nRows - 1, topLeftCol + 7)
    drawLine(xla, topLeftRow - 1, topLeftCol+6, topLeftRow - 1, topLeftCol + 7, "o", const.xlThin)

    #drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol + 6, topLeftRow + nRows - 1, topLeftCol + 7)
    #drawLine(xla, topLeftRow - 1, topLeftCol +6, topLeftRow - 1, topLeftCol + 7, "o", const.xlThin)

    
    for i in range(nRows):
        


        date   = res['outputDates'][i]
        zspread   = float(res['zcSpread'][i]*100.0)
        pypread = float(res['pySpread'][i]*100.0)
        pyRiskFree = float(res['pyRiskFree'][i]*100.0)


        xla.Cells(topLeftRow + i, topLeftCol).Value = date
        xla.Cells(topLeftRow + i, topLeftCol).NumberFormat = FORMATT
        xla.Cells(topLeftRow + i, topLeftCol).HorizontalAlignment = const.xlCenter
        xla.Cells(topLeftRow + i, topLeftCol + 1).Value = zspread
        xla.Cells(topLeftRow + i, topLeftCol + 1).NumberFormat = "0.00000"
        xla.Cells(topLeftRow + i, topLeftCol + 1).HorizontalAlignment = const.xlCenter

        xla.Cells(topLeftRow + i, topLeftCol + 3).Value = date
        xla.Cells(topLeftRow + i, topLeftCol + 3).NumberFormat = FORMATT
        xla.Cells(topLeftRow + i, topLeftCol + 3).HorizontalAlignment = const.xlCenter
        xla.Cells(topLeftRow + i, topLeftCol + 4).Value = pypread
        xla.Cells(topLeftRow + i, topLeftCol + 4).NumberFormat = "0.00000"
        xla.Cells(topLeftRow + i, topLeftCol + 4).HorizontalAlignment = const.xlCenter

        xla.Cells(topLeftRow + i, topLeftCol + 6).Value = date
        xla.Cells(topLeftRow + i, topLeftCol + 6).NumberFormat = FORMATT
        xla.Cells(topLeftRow + i, topLeftCol + 6).HorizontalAlignment = const.xlCenter
        xla.Cells(topLeftRow + i, topLeftCol + 7).Value = pyRiskFree
        xla.Cells(topLeftRow + i, topLeftCol + 7).NumberFormat = "0.00000"
        xla.Cells(topLeftRow + i, topLeftCol + 7).HorizontalAlignment = const.xlCenter
        



#=======================================================================================================================

def writeBondFittingRes1OnXls(crv, xla, str_boot_opt, res, codice_curva):


    nameSheet = nameSheetBondSurv
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


    Attributi_1 = \
        { "Data riferimento"        : crv.ref_date
        , "Descrizione"     : crv.description
        , "Valuta"        : crv.curr
        , "Tipo nodo"     : "Survival probability"
        , "Emittente"     : crv.emittente
        , "Recovery rate" : crv.recoveryRate
        , "Rating" : crv.rating
        , "Modello hazard rate" : crv.HRateModel
        , "Modello di valutazione" : crv.BondModel
        }

    Attributi_2 = \
        { "Data riferimento"        : crv.ref_date
        , "Descrizione"     : crv.description
        , "Valuta"        : crv.curr
        , "Tipo nodo"     : "Hazard rate"
        , "Emittente"     : crv.emittente
        , "Recovery rate" : crv.recoveryRate
        , "Rating" : crv.rating
        , "Modello hazard rate" : crv.HRateModel
        , "Modello di valutazione" : crv.BondModel

        }

    r2 = s.Range(xla.Cells(r.Row, r.Column + 3), xla.Cells(r.Row, r.Column + 3))
    r3 = s.Range(xla.Cells(r.Row, r.Column + 6), xla.Cells(r.Row, r.Column + 7))

    ra = intestazioneSwapCurveSegmenti(xla, s, r , Attributi_1, nCols=2, text = codice_curva)
    rb = intestazioneSwapCurveSegmenti(xla, s, r2, Attributi_2, nCols=2, text = codice_curva)
    #rc = intestazioneSwapCurveSegmenti(xla, s, r3, Attributi_3, nCols=2, text = codice_curva)

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
    
    #xla.Cells(topLeftRow - 1, topLeftCol + 6).Value = "Date"
    #xla.Cells(topLeftRow - 1, topLeftCol + 6).HorizontalAlignment = const.xlCenter
    #xla.Cells(topLeftRow - 1, topLeftCol + 7).Value = "Value"
    #xla.Cells(topLeftRow - 1, topLeftCol + 7).HorizontalAlignment = const.xlCenter

    

    nRows = len(res['outputDates'])

    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol, topLeftRow + nRows - 1, topLeftCol + 1)
    drawLine(xla, topLeftRow - 1, topLeftCol, topLeftRow - 1, topLeftCol + 1, "o", const.xlThin)

    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol + 3, topLeftRow + nRows - 1, topLeftCol + 4)
    drawLine(xla, topLeftRow - 1, topLeftCol+3, topLeftRow - 1, topLeftCol + 4, "o", const.xlThin)

    #drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol + 6, topLeftRow + nRows - 1, topLeftCol + 7)
    #drawLine(xla, topLeftRow - 1, topLeftCol +6, topLeftRow - 1, topLeftCol + 7, "o", const.xlThin)

    
    for i in range(nRows):

        date   = res['outputDates'][i]
        surv   = res['survProbCum'][i]
        #marg   = res['marginalDefault'][i]
        h_rate = res['hazardRate'][i]


        xla.Cells(topLeftRow + i, topLeftCol).Value = date
        xla.Cells(topLeftRow + i, topLeftCol).NumberFormat = FORMATT
        xla.Cells(topLeftRow + i, topLeftCol).HorizontalAlignment = const.xlCenter
        xla.Cells(topLeftRow + i, topLeftCol + 1).Value = surv
        xla.Cells(topLeftRow + i, topLeftCol + 1).NumberFormat = "0.00000"
        xla.Cells(topLeftRow + i, topLeftCol + 1).HorizontalAlignment = const.xlCenter

        xla.Cells(topLeftRow + i, topLeftCol + 3).Value = date
        xla.Cells(topLeftRow + i, topLeftCol + 3).NumberFormat = FORMATT
        xla.Cells(topLeftRow + i, topLeftCol + 3).HorizontalAlignment = const.xlCenter
        xla.Cells(topLeftRow + i, topLeftCol + 4).Value = h_rate
        xla.Cells(topLeftRow + i, topLeftCol + 4).NumberFormat = "0.00000"
        xla.Cells(topLeftRow + i, topLeftCol + 4).HorizontalAlignment = const.xlCenter

        """
        xla.Cells(topLeftRow + i, topLeftCol + 6).Value = date
        xla.Cells(topLeftRow + i, topLeftCol + 6).NumberFormat = FORMATT
        xla.Cells(topLeftRow + i, topLeftCol + 6).HorizontalAlignment = const.xlCenter
        xla.Cells(topLeftRow + i, topLeftCol + 7).Value = h_rate
        xla.Cells(topLeftRow + i, topLeftCol + 7).NumberFormat = "0.00000"
        xla.Cells(topLeftRow + i, topLeftCol + 7).HorizontalAlignment = const.xlCenter
        """

"""
def writeBondFittingRes2OnXls(crv, xla, str_boot_opt, res, codice_curva):

    
    nameSheet = nameSheetBondSpread
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

    boot_type_0   = int(crv.cds_boot_method)
    interp_type_0 = crv.rf_interp_type
    
    boot_type = crv.mapBootCDS(boot_type_0)
    interp_type = crv.mapCodeModelInv(interp_type_0)
  
    # -----
    Attributi_1 = \
        { "Data riferimento"        : crv.ref_date
        , "Descrizione"     : crv.description
        , "Valuta"        : crv.curr
        , "Tipo curva"    : crv.type
        , "Tipo nodo"     : "Par Yield spread"
        , "Emittente"     : crv.emittente
        , "Recovery rate" : crv.recovery
        , "Rating" : crv.rating
        , "Seniority" : crv.seniority
        , "Modello interpolante" : interp_type
        , "Tipo Bootstrap Hazard Rate" : boot_type
        }

    Attributi_2 = \
      { "Data riferimento"        : crv.ref_date
        , "Descrizione"     : crv.description
        , "Valuta"        : crv.curr
        , "Tipo curva"    : crv.type
        , "Tipo nodo"     : "Zero Coupon spread"
        , "Emittente"     : crv.emittente
        , "Recovery rate" : crv.recovery
        , "Rating" : crv.rating
        , "Seniority" : crv.seniority
        , "Modello interpolante" :interp_type
        , "Tipo Bootstrap Hazard Rate" : boot_type
        }


    r2 = s.Range(xla.Cells(r.Row, r.Column + 3), xla.Cells(r.Row, r.Column + 3))
    r3 = s.Range(xla.Cells(r.Row, r.Column + 6), xla.Cells(r.Row, r.Column + 7))
    
    ra = intestazioneSwapCurveSegmenti(xla, s, r , Attributi_1, nCols=2, text = codice_curva)
    rb = intestazioneSwapCurveSegmenti(xla, s, r2, Attributi_2, nCols=2, text = codice_curva)


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

    nRows = len(res['outputDates'])

    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol, topLeftRow + nRows - 1, topLeftCol + 1)
    drawLine(xla, topLeftRow - 1, topLeftCol, topLeftRow - 1, topLeftCol + 1, "o", const.xlThin)

    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol + 3, topLeftRow + nRows - 1, topLeftCol + 4)
    drawLine(xla, topLeftRow - 1, topLeftCol+3, topLeftRow - 1, topLeftCol + 4, "o", const.xlThin)
    
    for i in range(nRows):

        date   = res['outputDates'][i]
        pySpread   = res['pySpread'][i]
        zcSpread   = res['zcSpread'][i]

        xla.Cells(topLeftRow + i, topLeftCol).Value = date
        xla.Cells(topLeftRow + i, topLeftCol).NumberFormat = FORMATT
        xla.Cells(topLeftRow + i, topLeftCol).HorizontalAlignment = const.xlCenter
        xla.Cells(topLeftRow + i, topLeftCol + 1).Value = pySpread
        xla.Cells(topLeftRow + i, topLeftCol + 1).NumberFormat = "0.00"
        xla.Cells(topLeftRow + i, topLeftCol + 1).HorizontalAlignment = const.xlCenter

        xla.Cells(topLeftRow + i, topLeftCol + 3).Value = date
        xla.Cells(topLeftRow + i, topLeftCol + 3).NumberFormat = FORMATT
        xla.Cells(topLeftRow + i, topLeftCol + 3).HorizontalAlignment = const.xlCenter
        xla.Cells(topLeftRow + i, topLeftCol + 4).Value = zcSpread
        xla.Cells(topLeftRow + i, topLeftCol + 4).NumberFormat = "0.00000"
        xla.Cells(topLeftRow + i, topLeftCol + 4).HorizontalAlignment = const.xlCenter

"""












def writePortfoliOnXls(crv, nameSheet, xla, curve_type):
    # ---
    from DEF_intef import FORMATT
    # ---
    rangeStart = "B2"
    distCurve  = 1
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


    rangeStartNew = writeIntestazionePortfolio(xla, sheet , rOut, Attributi)
    rangeStartNew = writePortfolioAnag(xla, rangeStartNew, crv)


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


def readIntestazionePortfolio(xla , r , cc):
    
    row             = r.Row
    col             = r.Column

    cc.curr           = xla.Range(xla.Cells(row+1, col+1), xla.Cells(row+1, col+1)).Value
    dd = xla.Range(xla.Cells(row + 2, col + 1), xla.Cells(row + 2, col + 1)).Value
    cc.ref_date = datetime.date(year=dd.year, month=dd.month, day=dd.day)
    cc.description   = xla.Range(xla.Cells(row + 3, col + 1), xla.Cells(row + 3, col + 1)).Value
    cc.emittente     = xla.Range(xla.Cells(row + 4, col + 1), xla.Cells(row + 4, col + 1)).Value
    cc.rating        = xla.Range(xla.Cells(row + 5, col + 1), xla.Cells(row + 5, col + 1)).Value
    cc.recoveryRate  = xla.Range(xla.Cells(row + 6, col + 1), xla.Cells(row + 6, col + 1)).Value
    cc.source        = xla.Range(xla.Cells(row + 7, col + 1), xla.Cells(row + 7, col + 1)).Value

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
    
    """
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
    cc.node_type     = xla.Range(xla.Cells(row + 12, col + 1), xla.Cells(row + 12, col + 1)).Value
    cc.quotation     = xla.Range(xla.Cells(row + 13, col + 1), xla.Cells(row + 13, col + 1)).Value
    cc.rating        = xla.Range(xla.Cells(row + 14, col + 1), xla.Cells(row + 14, col + 1)).Value
    cc.recovery      = xla.Range(xla.Cells(row + 15, col + 1), xla.Cells(row + 15, col + 1)).Value
    cc.rendimento    = xla.Range(xla.Cells(row + 16, col + 1), xla.Cells(row + 16, col + 1)).Value
    cc.sector        = xla.Range(xla.Cells(row + 17, col + 1), xla.Cells(row + 17, col + 1)).Value
    cc.seniority     = xla.Range(xla.Cells(row + 18, col + 1), xla.Cells(row + 18, col + 1)).Value
    cc.source        = xla.Range(xla.Cells(row + 19, col + 1), xla.Cells(row + 19, col + 1)).Value

    """
    r = xla.Range(xla.Cells(row + 9, col), xla.Cells(row + 9, col))
    #cc.show()
    
    return r




def readPortfolioFromXls(xla, des, pos, nameSheet):
    rangeStart = "B2"
    sheet = xla.ActiveWorkbook.Sheets(nameSheet)
    r     = sheet.Range(rangeStart)
    

    row = r.Row
    col = r.Column
    i = 0
    k = 1
    
    for i in range(0, 1000):

        if (r.Value == des) and (k == pos):
            break
        
        else:
            
            r = xla.Range(xla.Cells(row + i, col), xla.Cells(row + i, col))
            i = i  + 1 
            if (r.Value == None):
                k = k + 1

    
    cc = BondPortfolio()
    r = readIntestazionePortfolio(xla, r, cc)
    r = readPortfolioAnag(xla, r, cc)
    
    return cc

def readPortfolioAnag(xla, r, cc):
    row  = r.Row + 1
    col  = r.Column
    name = r.Value
    
    
    field_list = ['Isin', 'Descrizione', 'Seniority', 'Tipo tasso', 'Tipo quotazione', 'Tipo rimborso', 'Prezzo emissione', 'Prezzo rimborso / Inflation Ratio', 'Data emissione',  'Data scadenza', 
        'Tempo scadenza', 'Giorni di fixing',    'Tipo fixing',    'Basis',    'Adjustment',    'Periodicita cedola (mesi)',    'Tenor del tasso floater (mesi)',
        'Tasso cedolare annuo (Fisso/spread)',    'Cedola in corso',    'Prezzo-MID',    'Prezzo-BID',    'Prezzo-ASK',    'YTM/DM (MID)',    'YTM/DM (BID)',    'YTM/DM (ASK)',
        'Tasso di riferiemnto',    'Tasso repo',    'Data prezzo di mercato',    'Contributor',    'Peso',    'Indicizzazione']

    field_date_list = ['Data emissione', 'Data scadenza', 'Data prezzo di mercato']

    n_fields = 31
    
    portfolio_anag = {}

    i = 0
    r_pesi = ''

    
    #while r.Value != None and r_pesi != None:
    while r_pesi != None:

        portfolio_anag[i] = {}
        
        for j in range(0, n_fields):
            
        
            fieldNameTmp = field_list[j]
            r = xla.Range(xla.Cells(row + i, col + j), xla.Cells(row + i, col + j))
            
            
            
            
            if (fieldNameTmp)in field_date_list:

                v0 = r.Value

                if (v0 == None):
                    fieldValueTmp = None
                else:
                    fieldValueTmp = datetime.datetime(year = v0.year, month = v0.month, day = v0.day)
            
            else:
                
                fieldValueTmp = r.Value

            portfolio_anag[i][fieldNameTmp] = fieldValueTmp
        r_pesi = portfolio_anag[i]['Peso']

        
        
        i = i + 1
        


    cc.portfolio_anag = portfolio_anag
    


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



