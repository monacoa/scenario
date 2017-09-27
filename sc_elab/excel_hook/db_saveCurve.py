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
        label_r      = "Zero Coupon"
        label_tn     = "Spot"
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
        msg = "Missing input sheet for Swap Curves in your workbook... \nNothing to do for me!"
        tkMessageBox.showinfo("Warning!", msg)
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


def checkParmsExistenceOnDb(curve, c):
    code = curve.getCurveCode()
    print "code", code
    qry = '''
            SELECT * from MKT_ParametriInterp where codice_modello = '%s'
          '''%code
    print qry
    c.execute(qry)
    res = c.fetchall()
    print "---->res:", res
    if (len(res)>0):
        return True
    return False


def deletingDbParms(code, conn):
    cursor = conn.db_anag()
    qry = '''
            delete from MKT_ParametriInterp_D where codice_modello = '%s'
          '''%code
    print qry
    cursor.execute(qry)
    qry = '''
                delete from MKT_ParametriInterp where codice_modello = '%s'
              ''' % code
    print qry
    cursor.execute(qry)
    conn.commit()

def saveParmsOnDB(curve, cursor):

    cod_strip = curve.getCurveCode()[1:3]
    qry = '''
                SELECT TIPO_MODELLO FROM MKT_Modello where COD = '%s'
            '''%cod_strip
    print qry
    cursor.execute(qry)

    res = cursor.fetchall()
    model_code = res[0][0]

    qry = '''
            INSERT INTO MKT_ParametriInterp
            (CODICE_MODELLO, DATA_RIF, DESCRIZIONE, VALUTA, QUOTAZIONE, TIPO_SCARICO, FONTE_DATO, TIPO_CURVA, TIPO_MODELLO, CODICE_EMITTENTE, CODICE_SETTORE, CODICE_RATING, CODICE_SENIORITY)
            VALUES
            ('%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s', '%s')
            '''%(curve.getCurveCode(),str(curve.ref_date), curve.description, curve.curr, curve.quotation, str(int(curve.download_type)), curve.source,
                                       curve.type, model_code, curve.emittente, curve.settore, curve.rating, curve.seniority)
    print qry
    cursor.execute(qry)

    print curve.parms_list
    print curve.interp_parms

    qry2 = '''INSERT INTO MKT_ParametriInterp_D (codice_modello, codice_parametro, valore) VALUES '''

    counter = 1
    sep = ""
    for par in curve.parms_list:
        for value in (curve.interp_parms[par]):
            qry2 += sep + "('%s', '%d', '%.8f')"%(curve.getCurveCode(), counter, value )
            counter += 1
            sep = ","
    

    print qry2
    cursor.execute(qry2)
    return True

# ---
# cc = classe curva, c = cursore
# ---
from connection import Connection
def saveInterpolationParmsOnDb(cc):
    connection = Connection()
    c = connection.db_anag()

    root2 = Tk()
    root2.withdraw()
    exst = checkParmsExistenceOnDb(cc, c)
    if exst:
        ans= tkMessageBox.askquestion("Unable to save Bootstrap results because they're already on DB.",
                                 "DELETING... Are You Sure?", icon='warning')
        if ans == 'yes':
            print "ho risposto yes: entro in deleting"
            deletingDbParms(cc.getCurveCode(), connection)
            # ---
        else:
            print "sono entrata nel caso no"
            msg = "Unable to save Bootstrap results because they're already on DB... Please delete IT before!!"
            tkMessageBox.showinfo("x@!#!", msg)
            connection.close()
            root2.destroy()
            return False

    print "ORA CERCO DI SALVARE I PARAMETRI"
    res = saveParmsOnDB(cc, c)
    connection.commit()
    if not res:
        msg = "Something's wrong saving parms.....!"
        tkMessageBox.showinfo("x@!#!", msg)
    else:
        msg = "Fitting results are on DB!"
        tkMessageBox.showinfo("YES WE CAN!", msg)
    connection.close()
    root2.destroy()
