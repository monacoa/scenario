import connection
from Tkinter import *
import tkMessageBox
from xls_bootCurve import readBootstrappedCurveFromXls

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

    #creo la connessione
    con = connection.Connection()
    cur = con.db_save()

    #1) leggo i dati da salvare
    Bootcurve = readBootstrappedCurveFromXls(xla, "", pos, nameSheet)

