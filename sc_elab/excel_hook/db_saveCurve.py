import connection
from Tkinter import *
import tkMessageBox
from xls_bootCurve import readBootstrappedCurveFromXls


def checkCurveExistenceOnDb(Bootcurve, ZC, DF):
    con = connection.Connection()
    cur = con.db_save()

    codeL, codeR = Bootcurve.getCurveCode()
    str = ""
    sep = ""
    if (ZC == True):
        codeZC = codeL+"ST"+codeR+"_"+ Bootcurve.code
        str += codeZC
        sep = ","
    if (DF == True):
        codeDF = codeL + "DF" + codeR + "_" + Bootcurve.code
        str += sep + codeDF

    qry = '''
            SELECT * from MKT_Curve where codice_curva in (%s)
          '''%str

    print qry


def saveZcDfOnDB(xla, nameSheet, code, pos, DF, ZC):
    try:
        book = xla.ActiveWorkbook
        s = book.Sheets(nameSheet)
        s.Activate()
    except:
        root = Tk()
        msg = "Missing input sheet for Swap Curves in your workbook... \nNothing to do for me!"
        tkMessageBox.showinfo("Warning!", msg)
        root.destroy()
        return

    #1) leggo i dati da salvare
    Bootcurve = readBootstrappedCurveFromXls(xla, "", pos, nameSheet)
    print  "MOSTRO LA CURVA!!"
    Bootcurve.show()
    print "##################################"
    checkCurveExistenceOnDb(Bootcurve, ZC, DF)
