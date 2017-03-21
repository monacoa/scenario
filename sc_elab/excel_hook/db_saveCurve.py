import connection
from Tkinter import *
import tkMessageBox
from xls_bootCurve import readBootstrappedCurveFromXls


def qrySaveAnag (curve,code,label_r, label_tn, cursor):
    curve.show()
    qry =  '''
	        INSERT INTO MKT_Curve
		    (codice_rating, codice_emittente, valuta, data_rif, fonte_dato, descrizione, rendimento, tipo_curva,
		    codice_settore, quotazione, tipo_scarico, codice_seniority, tipo_nodo, codice_curva)
		    VALUES
		    ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s',
		   '%s', '%s', '%s', '%s', '%s', '%s')
		  ''' %(curve.rating
                , curve.emittente
                , curve.curr
                , curve.ref_date
                , curve.source
                , curve.description
                , label_r
                , curve.type
                , curve.settore
                , curve.quotation
                , int(curve.download_type)
                , curve.seniority
                , label_tn
                , code)
    print qry
    cursor.execute(qry)

def qrySaveData(curve, code , ds, vs, cursor):

    dd      = "%02d" % curve.ref_date.day
    mm      = "%02d" % curve.ref_date.month
    yy      = str(curve.ref_date.year)
    settl   = yy + "-" + mm + "-" + dd

    sep = ""
    sr = ""
    for d,v, in zip(ds, vs):
        dd      = "%02d" % d.day
        m       = "%02d" % d.month
        y       = str(d.year)
        term    = y + "-" + m + "-" + dd
        tup     = "( '" + code + "', '" + settl + "', '" + term + "' , '%.8f')"%v
        sr     += sep + tup
        sep     = ","
    qry ='''
        INSERT INTO MKT_Curve_D ( CODICE_CURVA, SETTLEMENT, TERM, VALORE)
        VALUES %s
        '''%sr
    print qry
    cursor.execute(qry)



def saveCurve(curve, ZC, DF, SP, cur):
    if ZC:
        codeL, codeR = curve.getCurveCode()
        code         = codeL + "ZC" + codeR + "_" + curve.code
        #print code
        label_r      = "Zero Coupon"
        label_tn     = "Spot"
        #print label_r
        #print label_tn
        qrySaveAnag(curve, code, label_r, label_tn, cur)
        dates        = curve.boot_dates
        values       = curve.boot_rates
        value2       = [rate/100. for rate in values]
        qrySaveData(curve, code, dates, value2, cur)
    if DF:
        codeL, codeR = curve.getCurveCode()
        code         = codeL + "DF" + codeR + "_" + curve.code
        label_r      = "Zero Coupon"
        label_tn     = "Discount Factor"
        qrySaveAnag(curve, code, label_r, label_tn, cur)
        dates        = curve.boot_dates
        values       = curve.boot_df
        qrySaveData(curve, code, dates, values, cur)
    if SP:
        xxxx


def checkCurveExistenceOnDb(curve, cur, ZC = False, DF = False, SP = False):
    codeL, codeR = curve.getCurveCode()
    print "codeL, codeR", codeL, codeR
    str = ""
    sep = ""
    # ---

    codes = []
    if (ZC == True):
        codeZC = "'"+codeL+"ZC"+codeR+"_"+ curve.code+"'"
        codes.append(codeZC)
        str += codeZC
        sep = ","


    if (DF == True):
        codeDF = "'"+codeL + "DF" + codeR + "_" + curve.code + "'"
        codes.append(codeDF)
        str += sep + codeDF
    print str

    if (SP == True) : zzzzzzzzzzzzzzzzz

    qry = '''
            SELECT * from MKT_Curve where codice_curva in (%s)
          '''%str
    print qry
    cur.execute(qry)
    res = cur.fetchall()
    print "---->res:", res
    if (len(res)>0):
        return True, str
    return False, None



def deletingDbCurves(codes):
    con = connection.Connection()
    cur = con.db_save()

    print "codici da cancellare:", codes
    qry = ''' delete FROM MKT_Curve_D where codice_curva in (%s)'''%codes
    print qry
    cur.execute(qry)
    qry = ''' delete FROM MKT_Curve where codice_curva in (%s)''' % codes
    print qry
    cur.execute(qry)
    con.commit()
    con.close()

def saveZcDfOnDB(xla, nameSheet, code, pos, DF, ZC):
    try:
        book = xla.ActiveWorkbook
        s = book.Sheets(nameSheet)
        s.Activate()
    except:
        #root = Tk()
        msg = "Missing input sheet for Swap Curves in your workbook... \nNothing to do for me!"
        tkMessageBox.showinfo("Warning!", msg)
        #root.destroy()
        return False

    #1) leggo i dati da salvare
    Bootcurve = readBootstrappedCurveFromXls(xla, "", pos, nameSheet)

    con = connection.Connection()
    cur = con.db_save()
    exist, codes = checkCurveExistenceOnDb(Bootcurve, cur, ZC, DF, False)
    if exist :
        return False, codes

    else:
        saveCurve(Bootcurve, ZC, DF, False, cur)
        con.commit()
        con.close()
        return True, None




