
from pyxll import xl_func, xl_app, xl_menu, xl_macro, xlcAlert
import win32api
import win32con

from Tkinter import *
import tkMessageBox
import Tkinter
import pyodbc
import datetime
import time
from connection import *




def donothing():
    tkMessageBox.showinfo("Nothing To do", "bye bye")
    # filewin = Toplevel(root)
    # button = Button(filewin, text="Do nothing button")
    # button.pack()
    # filewin.mainloop()



def getDatesListFromDb():
    cn = Connection()
    cur = cn.db_data()
    qry = '''
                     select distinct Data from alm.DProTS_master where BloombergTicker in
                     ( select distinct BloombergTicker from alm.DProCurve where TipoDato in
                         ('CDepositi', 'CSwap', 'CLibor', 'CFuture')
                      )
                     order by Data desc
                    '''
    cur.execute(qry)
    res = cur.fetchall()
    date_list = []
    for record in res:
        date_list.append(record[0])
    return date_list



class W_curveType (Frame):
    def __init__(self, master = None):

        Frame.__init__(self, master)
        self.master = master
        self.menubar = Menu(self)
        filemenu = Menu(self.menubar, tearoff=0)
        filemenu.add_command(label="swap - RF curve", command=self.load_swp_RF_curve)
        filemenu.add_command(label="swap - EUR multitenor curve", command=donothing)

        filemenu.add_command(label="govt curve", command=donothing)
        filemenu.add_command(label="sector curve", command=donothing)

        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=donothing)

        self.menubar.add_cascade(label="select Curve Type", menu=filemenu)
        self.master.config(menu=self.menubar)
        self.canvas = Canvas(self, bg="grey", width=200, height=200,
                             bd=0, highlightthickness=0)
        self.canvas.pack()

    def load_swp_RF_curve(self):
        self.new_window = W_curveDate(self)

    def close_window(self):
        self.destroy()



class W_curveDate(LabelFrame):

    def __init__(self, parent = None):
        c_date = None
        if parent:
            self.master = parent.master
            parent.close_window()
        # create labelframe
        #self.master.geometry("360x480")
        LabelFrame.__init__(self, self.master)
        self.config(text = "Select reference date:")#, bg="white", width=400, height=400 )
        self.pack(fill="both", expand="yes")

        # create scrollbar
        self.bar = Scrollbar(self)
        # create mylist
        self.mylist = Listbox(self, yscrollcommand=self.bar.set)
        dl = getDatesListFromDb()
        for d in dl:
           date = d.date()
           self.mylist.insert(END, str(date))

        self.mylist.pack(side=LEFT, fill='y')
        self.bar.pack(side=LEFT, fill='y')

        self.mylist.config(yscrollcommand=self.bar.set)
        self.bar.config(command=self.mylist.yview)

        # -----------
        # cretae button
        self.btn2 = Button(self, text="Cancel",  command=self.close_window)
        self.btn2.pack(side=BOTTOM, fill = 'x')
        # cretae button
        self.btn1 = Button(self, text="Select", command=self.selected_date)
        self.btn1.pack(side=BOTTOM, fill='x')



    def close_window(self):
        self.destroy()


    def selected_date(self):
        #recupero la data selezionata
        c_date = str((self.mylist.get(ACTIVE)))
        print "====value", c_date
        self.date = c_date
        self.new_window = W_curveSelection (parent = self, curve_date=self.date)



def getCurvesListFromDb(curve_date):
    c_date_qry = str(curve_date).replace("-", "")
    con = Connection()
    qry = '''
          select distinct BloombergTicker FROM alm.DProTS_master where BloombergTicker in
          ( select distinct BloombergTicker from alm.DProCurve
            where TipoDato in ('CDepositi', 'CSwap', 'CLibor', 'CFuture')
          )and Data = '%s'
          ''' % c_date_qry
    print qry
    c_data = con.db_data()
    c_data.execute(qry)
    tickers = c_data.fetchall()
    tkrs = ""
    sep  = ""
    for tic in tickers:
        tkrs += sep + "'" + tic[0] + "'"
        sep = ","
    #================
    c_anag = con.db_anag()
    qry = '''
          SELECT DISTINCT description, ID_object FROM AnagMktObject WHERE ID_object in
          ( SELECT DISTINCT ID_object FROM AnagMktObject_D WHERE ticker IN
            (%s
            )
          )ORDER BY description
          ''' % tkrs
    print "qry:", qry
    c_anag.execute(qry)
    res = c_anag.fetchall()
    print res
    ll = []
    for c in res:
        print "c", c
        curva = str(c[0])
        id    = int(c[1])
        a =[curva,id]
        ll.append(a)
        print "ll", ll
    con.close()
    return ll
    # -----------


class W_curveSelection (LabelFrame):
    def __init__(self, parent = None, curve_date =None):
        c_date = None
        if parent:
            self.master = parent.master
            parent.close_window()
        LabelFrame.__init__(self, self.master)
        self.config(text="Available curves for the selected date:")
        self.pack(fill="both", expand="yes")
        # create scrollbar
        self.bar = Scrollbar(self)
        # create mylist
        self.mylist = Listbox(self, yscrollcommand=self.bar.set)
        dl = getCurvesListFromDb(curve_date)
        for l in dl:
            crv = l[0]
            self.mylist.insert(END, crv)

        self.mylist.pack(side=LEFT, fill='y')
        self.bar.pack(side=LEFT, fill='y')
        self.mylist.config(yscrollcommand=self.bar.set)
        self.bar.config(command=self.mylist.yview)
        # -----------
        # cretae button
        self.btn2 = Button(self, text="Cancel", command=self.close_window)
        self.btn2.pack(side=BOTTOM, fill='x')
        # cretae button
        self.btn1 = Button(self, text="Select", command=self.selected_curve)
        self.btn1.pack(side=BOTTOM, fill='x')

    def close_window(self):
        self.destroy()

    def selected_curve(self):
        #recupero la data selezionata
        curve_des = str((self.mylist.get(ACTIVE)))
        self.curve = curve_des
        #self.close_window()



from sc_elab.core.SwpCurve import *


def loadCurvedataFromDB(crv):

    con = Connection()
    c_a = con.db_anag()
    qry = "SELECT DISTINCT  ID_object, Tipo, rating FROM AnagMktObject WHERE description = '%s' "%(crv.description)
    print qry
    c_a.execute(qry)
    res = c_a.fetchall()
    print qry
    if len(res)!= 1:
        aaaaaaaaaaaaaaaaaaaaaaa
    crv.rating = res[0][2]
    crv.type = res[0][1]
    #----

    qry = '''
            SELECT ticker FROM AnagMktObject_D where ID_object ='%s'
            '''%res[0][0]
    print qry

    c_a.execute(qry)

    res = c_a.fetchall()
    blm_tckrs_lst_str  = ""
    sep = ""

    for record in res:
        blm_tckrs_lst_str+= sep+ "'" + record[0] +"'"
        sep = ","

    #---------
    # RECUPERO CURRENCY E NUMERO/TIPOLOGIA DI SEGMENTI
    #---------

    c_d = con.db_data()

    qry = '''
         SELECT DISTINCT currency, TipoDato FROM alm.DProCurve where BloombergTicker IN
         (
            SELECT DISTINCT BLOOMBERGTICKER FROM alm.DProTS_master WHERE
            BLOOMBERGTICKER IN (%s) AND
             DATA = '%s'
            )
        '''%(blm_tckrs_lst_str, str(crv.ref_date).replace("-", "") )

    print qry

    c_d.execute(qry)
    res = c_d.fetchall()
    print "********************************************************"
    print qry
    print res
    print "********************************************************"

    crv.curr = record[0][0]
    if   crv.curr == "EUR" : crv.floater_tenor = "6M"
    elif crv.curr == "USD" : crv.floater_tenor = "3M"

    if   crv.quotation == "MID": quotazione = "valoremid"
    elif crv.quotation == "ASK": quotazione = "valoreask"
    elif crv.quotation == "BID": quotazione = "valorebid"
    else: mmmmmmmmmmmmmmmm

    for record in res:
        segm = record[1]
        qry = ""
        if ((segm == "CDepositi") or (segm == "CSwap") or (segm == 'CLibor')):

            qry = '''
                    SELECT      alm.DProCurve.maturityInt, alm.DProTS_master.%s
                    FROM        alm.DProCurve
                    INNER JOIN  alm.DProTS_master
                    ON          alm.DProCurve.BloombergTicker = alm.DProTS_master.BloombergTicker
                    WHERE       alm.DProTS_master.data='%s'
                    AND         alm.DProCurve.TipoDato = '%s'
                    AND         alm.DProCurve.BloombergTicker in (%s)
                    ORDER BY    alm.DProCurve.maturityInt
                    '''%(quotazione,str(crv.ref_date).replace("-", ""), segm, blm_tckrs_lst_str)
        elif segm == "CFuture":
            qry = '''
                    SELECT      alm.DProTS_master.ScadenzaFuture, alm.DProTS_master.%s
                    FROM        alm.DProCurve
                    INNER JOIN  alm.DProTS_master
                    ON          alm.DProCurve.BloombergTicker = alm.DProTS_master.BloombergTicker
                    WHERE       alm.DProTS_master.data= '%s'
                    AND         alm.DProCurve.TipoDato ='%s'
                    AND         alm.DProCurve.BloombergTicker in (%s)
                    ORDER BY    alm.DProTS_master.ScadenzaFuture
                    '''%(quotazione,str(crv.ref_date).replace("-", ""), segm, blm_tckrs_lst_str)

        else: xxxxxxxxxxxxxxxxxxxxxxxx

        print qry
        c_d.execute(qry)
        res = c_d.fetchall()
        print "segm:", segm, "RES:", res
        for record in res:
            l = [record[0], record[1]]
            (crv.raw_data[dict_segm[segm]]).append(l)

        print "DIZIONARIO:", crv.raw_data

    con.close()



def findRigthPlaceBootCurveSeg(xla, r, distCurve, dir="O"):
    rOut = None
    if dir == "V" :
        while (r.Value != None):
            righe = r.Offset(1, 0).Rows.count()
            r = r.Offset(righe + distCurve, 0)
        rOut = r
    else:
        while (r.Value != None):
                col= r.Offset(1, 0).Columns.count()
                r = r.Offset(-1, col + distCurve)
        rOut = r
    #-----
    if (rOut == None):
        msg = "Unable to compute the output range for your curve"
        popup_messagebox(msg)
        sys.exit()

    return rOut



import win32com.client
from win32com.client import constants as const

def drawBox(xla, r , spessore , rTopLeft = 0, cTopLeft = 0,rBottomRight=0,cBottomRight=0, Colore=0):
    print "----------------------------------------------"
    if (rTopLeft <= 0) or(cTopLeft <= 0) or (rBottomRight <= 0)or (cBottomRight <= 0):
        msg = "Le coordinate del box devono essere maggiori di zero"
        popup_messagebox(msg)
        sys.exit()
    print ":::::::::::::::::::::::::::::::::::::::::::"

    RR = xla.Range(xla.Cells(rTopLeft, cTopLeft), xla.Cells(rBottomRight, cBottomRight))
    RR.Borders(const.xlDiagonalDown).LineStyle = const.xlNone
    RR.Borders(const.xlDiagonalUp).LineStyle = const.xlNone
    RR.Borders(const.xlEdgeLeft).LineStyle = const.xlContinuous
    RR.Borders(const.xlEdgeLeft).Weight = spessore
    RR.Borders(const.xlEdgeLeft).ColorIndex = Colore


    RR.Borders(const.xlEdgeTop).LineStyle = const.xlContinuous
    RR.Borders(const.xlEdgeTop).Weight = spessore
    RR.Borders(const.xlEdgeTop).ColorIndex = Colore
    print ";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;", "FQ!"

    RR.Borders(const.xlEdgeBottom).LineStyle = const.xlContinuous
    RR.Borders(const.xlEdgeBottom).Weight = spessore
    RR.Borders(const.xlEdgeBottom).ColorIndex = Colore

    RR.Borders(const.xlEdgeRight).LineStyle = const.xlContinuous
    RR.Borders(const.xlEdgeRight).Weight = spessore
    RR.Borders(const.xlEdgeRight).ColorIndex = Colore



def formatTestataCurva(xla,nRiga,nColonna,nLarghezza,testo):

    (xla.Range(xla.Cells(nRiga, nColonna), xla.Cells(nRiga, nColonna + nLarghezza - 1))).Select()

    xla.Selection.HorizontalAlignment = const.xlCenter
    xla.Selection.VerticalAlignment = const.xlBottom
    xla.Selection.WrapText = False
    xla.Selection.Orientation = 0
    xla.Selection.AddIndent = False
    xla.Selection.IndentLevel = 0
    xla.Selection.ShrinkToFit = False
    xla.Selection.ReadingOrder = const.xlContext
    xla.Selection.MergeCells = False

    xla.Selection.Merge()
    xla.Selection.Font.ColorIndex = 2
    xla.Selection.Font.Bold = True

    xla.Selection.Interior.ColorIndex = 55
    xla.Selection.Interior.Pattern = const.xlSolid

    xla.Selection.Value = testo



def intestazioneSwapCurveSegmenti( xla, sheet, rng,  attributi,nCols = 2):

    nRows = len(attributi.keys())
    topLeftRow = rng.Row
    topLeftCol = rng.Column
    #disegno il box con intestazione
    #addS = rng.Address
    #addE = sheet.Cells(topLeftRow + nRows, topLeftCol+nCols).Address
    #print "addS, addE:", addS,addE

    #rng_intest = sheet.Range(rng.Address+":")
    #exit()

    drawBox            (xla, rng,3,topLeftRow, topLeftCol,topLeftRow + nRows, topLeftCol + nCols - 1, 0)
    formatTestataCurva (xla, topLeftRow, topLeftCol, nCols, attributi["Description"])
    # '
    # ' Scrivo contenuto Intestazione
    # '
    # With
    print topLeftRow, topLeftCol
    #rtemp = xla.Range(xla.Cells(topLeftRow + 1, topLeftCol), xla.Cells(topLeftRow + 1, topLeftCol))


    kk = (attributi.keys())
    kk.sort()


    i = 0
    for k in kk:
        print "i:", i
        a= xla.Cells(topLeftRow + 1+ i, topLeftCol)
        print "a:", a.Address
        a.Value = k
        b = xla.Cells(topLeftRow + 1+ i, topLeftCol+1)
        print "b:", b.Address
        b.Value = attributi[k]
        b.HorizontalAlignment = const.xlCenter
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
    drawBox(xla, xla.Range(rangeStart), const.xlMedium, topLeftRow, topLeftCol, topLeftRow + nRows, topLeftCol + nCols - 1, 0)
    #     #DrawLine
    #     rTopLeft:=topLeftRow + 1, _
    #     cTopLeft:=topLeftCol + 1, _
    #     rBottomRight:=topLeftRow + nRows, _
    # cBottomRight:=topLeftCol + 1, _
    # Vertical:=True
    # '
    # ' intestazione
    formatTestataCurva(xla, topLeftRow, topLeftCol, nCols, nomeBox)
    # Scirttura descrizione e parametri

    xla.Cells(topLeftRow + 1, topLeftCol).Value = "mean revertion speed"
    xla.Cells(topLeftRow + 1, topLeftCol+1).Value = "0.0"
    xla.Cells(topLeftRow + 1, topLeftCol+2).Value = "volatility"
    xla.Cells(topLeftRow + 1, topLeftCol+3).Value = "0.0"

    rangeStartN = xla.Cells(topLeftRow + 2, topLeftCol).Address
    print "rangeStartN", rangeStartN
    return rangeStartN


'''
def segmentoSwapCurve(xla, code,rangeSN, griglia, dateQuot, dateDefault nodiQuot  = Empty):
    rangeStart = rangeSN

    topLeftRow = xla.Range(rangeStart).Row
    topLeftCol = Range(rangeStart).Column
    nNodi = len (dateQuot)
    nRows = nNodi + 1
    nCols = 4
    '
    If IsEmpty(nodiQuot) Then ReDim nodiQuot(0 To nNodi - 1) As String
    '
    ' box intestazione
    '
    DrawBox rTopLeft:=topLeftRow, _
            cTopLeft:=topLeftCol, _
            rBottomRight:=topLeftRow + nRows, _
            cBottomRight:=topLeftCol + nCols - 1, _
            spessore:=xlMedium
    '
    ' Formatto zona codice curva
    '
    Select Case code
        Case "G"
            nomeSegmento = "0. Short term swap"
        Case "D"
            nomeSegmento = "1. Depositi"
        Case "L"
            nomeSegmento = "2. Libor"
        Case "F"
            nomeSegmento = "3. Futures"
        Case "S"
            nomeSegmento = "4. Swap Rate"
        Case Else
            nomeSegmento = code
    End Select
    '
    FormatTestataCurva topLeftRow, topLeftCol, nCols, nomeSegmento
    SegmentoSwapCurve = False
    '
    ' Linea orizzontale di separazione
    '
    DrawLine rTopLeft:=topLeftRow + 1, _
            cTopLeft:=topLeftCol, _
            rBottomRight:=topLeftRow + 1, _
            cBottomRight:=topLeftCol + nCols - 1
    '
    ' Linea verticale di separazione
    '
    DrawLine rTopLeft:=topLeftRow + 1, _
            cTopLeft:=topLeftCol + nCols - 2, _
            rBottomRight:=topLeftRow + nRows, _
            cBottomRight:=topLeftCol + nCols - 2, _
            Vertical:=True
    '
    ' scrittura dati
    '
    With Range(rangeStart).Offset(1, 0)
        '
        ' intestazione campi
        '
        .value = "Nodo"
        .Offset(0, 1).value = "Data scadenza"
        .Offset(0, 2).value = "Valori"
        .Offset(0, 3).value = "Usa Nodo"
        '
        ' scrittura nodi
        '
        For i = 1 To nNodi
            .Offset(i, 0).value = griglia(i - 1)
            .Offset(i, 1).value = CDate(dateQuot(i - 1)): .Offset(i, 1).NumberFormat = "dd-mm-yyyy"
            .Offset(i, 2).value = nodiQuot(i - 1):
            .Offset(i, 2).NumberFormat = "0.00%": If code = "F" Then .Offset(i, 2).NumberFormat = "0.00"
            .Offset(i, 3).value = dateDefault(i - 1)
        Next
        '
        .CurrentRegion.HorizontalAlignment = xlCenter
        '
    End With
    '
    ' aggiorno rangeStart
    '
    rangeStartNew = Range(rangeStartNew).Offset(nRows + 2, 0).Address
    SegmentoSwapCurve = True
    Exit Function

'''
#=======================================================================================================================
def writeCurveOnXls(crv, nameSheet, xla):
    print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ECCOMI!!!!!!!!!!"
    rangeStart = "B2"
    distCurve  = 3
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

    print "::::::::::::::::::::::::::::::::::::::::"


    rangeStartNew = intestazioneSwapCurveSegmenti ( xla, sheet , rOut, Attributi)

    # Genero il blocco per i parametri di Hull e White

    displayHWParamSwCurve (xla, rangeStartNew, Attributi)

    # ' Genero i blocchi dei segmenti della curva
    # codes = "DLGFS"

    cd = "DLGFS"
    for j in  range (len(cd)):
        code =cd[j]
        for s in crv.raw_data.keys():
            if code == s[0]:
                # visualizzo segmento
                print "segmento:", s
                print "code:", code
                N = len(crv.raw_data[s])
                listaYes = ["Y" for n in range(N)]
                print "listaY", listaYes
                #segmentoSwapCurve (code, rangeStartNew, s("MATURITY"), s("DATE"), vettoreYes, nodiQuot:=s("NODI"))

    # '
    # rOut.CurrentRegion.Columns.AutoFit
    # DisplaySwapCurve = True
    # Exit Function

#=======================================================================================================================
#=======================================================================================================================




@xl_func
def load_swap_curve_from_db(control):
    nameSheet = "Curvette"
    xla = xl_app()
    book = xla.ActiveWorkbook
    #--- creo foglio Curvette se non esiste
    try:
        s = book.Sheets(nameSheet)
    except:
        s = book.Sheets.Add()
        s.Name = nameSheet
   #------------------
    root = Tk()


    app = W_curveType(root)
    root.mainloop()


    curve_des = app.new_window.new_window.curve
    curve_date= app.new_window.date


    print "------------FQ ----------------------------------------:", curve_des, curve_date

    cc = Curve()
    cc.ref_date = curve_date
    cc.description= curve_des
    loadCurvedataFromDB(cc)
    print "CURVA CARICATA!"
    print "calcolo le date!"

    print "ora richiamo la scrittura!!!!!!!!!!!!!!!!!!!!!!!!"
    writeCurveOnXls(cc, nameSheet, xla)





@xl_func
def popup_messagebox(msg):
    xlcAlert(msg)