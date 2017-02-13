
from pyxll import xl_func, xl_app, xl_menu, xl_macro, xlcAlert
from connection import *
from Tkinter import *
import tkMessageBox
from db_qrys import getCurvesListFromDb
from sc_elab.core.SwpCurve import *
from win32com.client import constants as const

@xl_func
def popup_messagebox(msg):
    xlcAlert(msg)


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


#-------------------
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

        self.menubar.add_cascade(label="Click & Select Curve Type", menu=filemenu)
        self.master.config(menu=self.menubar)
        self.canvas = Canvas(self, bg="grey", width=200, height=200,
                             bd=0, highlightthickness=0)
        self.canvas.pack()

    def load_swp_RF_curve(self):
        self.new_window = W_curveDate(self)

    def close_window(self):
        self.destroy()


#-------------
class W_curveDate(LabelFrame):

    def __init__(self, parent = None):
        c_date = None
        if parent:
            self.master = parent.master
            parent.close_window()
        # create labelframe
        LabelFrame.__init__(self, self.master)
        self.config(text = "Select reference date:")
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

        self.date = c_date
        self.new_window = W_curveSelection (parent = self, curve_date=self.date)



#----------------
class W_curveSelection (LabelFrame):
    def __init__(self, parent = None, curve_date = None):
        c_date = None
        if parent:
            self.master = parent.master
            parent.close_window()
        LabelFrame.__init__(self, self.master)
        #self.geometry("800x600")
        #self.master.geometry("400x500")
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
        self.mylist.config(width=50)
        self.mylist.pack(side=LEFT, fill='y')

        self.bar.pack   (side=LEFT, fill='y')
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
        self.master.destroy()

    def selected_curve(self):
        #recupero la data selezionata
        curve_des = str((self.mylist.get(ACTIVE)))
        self.curve = curve_des
        self.close_window()




def findRigthPlaceBootCurveSeg(xla, r, distCurve, dir="O"):
    rOut = None
    #if dir == "V" :
    #     while (r.Value != None):
    #         righe = r.Offset(1, 0).Rows.count()
    #         r = r.Offset(righe + distCurve, 0)
    #     rOut = r
    # else:

    while (r.Value != None):
        print "r.Value----------", r.Value
        nCols = r.Columns.Count
        row =r.Row
        col =r.Column
        print nCols, row, col, distCurve
        r = xla.Range(xla.Cells(row, col + distCurve), xla.Cells(row, col + (nCols-1) + distCurve))


    rOut = r
    #-----
    if (rOut == None):
        msg = "Unable to compute the output range for your curve"
        popup_messagebox(msg)
        sys.exit()

    return rOut





def drawBox(xla, spessore , rTopLeft = 0, cTopLeft = 0,rBottomRight=0,cBottomRight=0, Colore=0):

    if (rTopLeft <= 0) or(cTopLeft <= 0) or (rBottomRight <= 0)or (cBottomRight <= 0):
        msg = "Le coordinate del box devono essere maggiori di zero"
        popup_messagebox(msg)
        sys.exit()


    RR = xla.Range(xla.Cells(rTopLeft, cTopLeft), xla.Cells(rBottomRight, cBottomRight))
    RR.Borders(const.xlDiagonalDown).LineStyle = const.xlNone
    RR.Borders(const.xlDiagonalUp).LineStyle = const.xlNone
    RR.Borders(const.xlEdgeLeft).LineStyle = const.xlContinuous
    RR.Borders(const.xlEdgeLeft).Weight = spessore
    RR.Borders(const.xlEdgeLeft).ColorIndex = Colore


    RR.Borders(const.xlEdgeTop).LineStyle = const.xlContinuous
    RR.Borders(const.xlEdgeTop).Weight = spessore
    RR.Borders(const.xlEdgeTop).ColorIndex = Colore


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

    drawBox            (xla, 3,topLeftRow, topLeftCol,topLeftRow + nRows, topLeftCol + nCols - 1, 0)
    formatTestataCurva (xla, topLeftRow, topLeftCol, nCols, attributi["Description"])

    kk = (attributi.keys())
    kk.sort()
    i = 0
    for k in kk:
        a= xla.Cells(topLeftRow + 1+ i, topLeftCol)
        a.Value = k

        b = xla.Cells(topLeftRow + 1+ i, topLeftCol+1)

        b.Value = attributi[k]
        print a, b.Value, type(attributi[k])
        if (type(attributi[k]) == datetime.datetime) or (type(attributi[k]) == datetime.date): b.NumberFormat = "dd-mm-yyyy"
        b.HorizontalAlignment = const.xlCenter
        i+=1
    rangeStart = xla.Cells(topLeftRow + nRows + 2, topLeftCol).Address
    return rangeStart


def drawLine (xla, rTopLeft, cTopLeft, rBottomRight, cBottomRight, hor, spessore):

    (xla.Range(xla.Cells(rTopLeft, cTopLeft), xla.Cells(rBottomRight, cBottomRight))).Select()

    if (hor =="o"):
        xla.Selection.Borders(const.xlEdgeBottom).LineStyle = const.xlContinuous
        xla.Selection.Borders(const.xlEdgeBottom).Weight = spessore
        xla.Selection.Borders(const.xlEdgeBottom).ColorIndex = const.xlAutomatic

    else:
        xla.Selection.Borders(const.xlEdgeRight).LineStyle = const.xlContinuous
        xla.Selection.Borders(const.xlEdgeRight).Weight = spessore
        xla.Selection.Borders(const.xlEdgeRight).ColorIndex = const.xlAutomatic


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
    if  code == "G":
            nomeSegmento = "0. Short term swap"
    elif code == "D":
            nomeSegmento = "1. Depositi"
    elif code == "L":
            nomeSegmento = "2. Libor"
    elif code =="F":
            nomeSegmento = "3. Futures"
            isFut = True
    elif code == "S":
            nomeSegmento = "4. Swap Rate"
    else:
            nomeSegmento = code

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

    #recupero linguaggio xls per formato date
    #lngCode = xla.LanguageSettings.LanguageID(const.msoLanguageIDExeMode)
    #print lngCode,lngCode,lngCode,lngCode


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
                if (type(ll[j]) == datetime.date) or (type(ll[j]) == datetime.datetime): a.NumberFormat = "dd-mm-yyyy"
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
                print "*" * 120
                print "segmento:", s
                print "code:", code
                print "chiamo la funzione dei segmenti!"
                rangeStartNew = segmentoSwapCurve (xla, rangeStartNew, code, crv.segms[s])

    # '
    # rOut.CurrentRegion.Columns.AutoFit
    # DisplaySwapCurve = True
    # Exit Function

def readCurvesNames(xla, s, rangeStart, direzione, distanza):
    r = xla.Range(rangeStart)
    print "sono qui!"
    if direzione.lower() == "o":
        curveL = []
        while (r.Value != None):
            print r
            print "r.Value----------", r.Value
            nomeCurva = r.Value
            curveL.append(nomeCurva)

            nCols = r.Columns.Count
            row   = r.Row
            col   = r.Column
            r = xla.Range(xla.Cells(row, col + distanza), xla.Cells(row, col + (nCols - 1) + distanza))
    return curveL



class W_bootstrapSelection (LabelFrame):

    def __init__(self, master = None, curveL = []):
        c_date = None
        LabelFrame.__init__(self, master)
        self.master = master
        #self.geometry("800x600")
        #self.master.geometry("400x500")
        self.config(text="Available curves for bootstrapping:")
        self.pack(fill="both", expand="yes")
        # create scrollbar
        self.bar = Scrollbar(self)

        # create mylist
        self.mylist = Listbox(self, yscrollcommand=self.bar.set)

        for l in curveL:
            crv = l
            self.mylist.insert(END, crv)
        self.mylist.config(width=50)
        self.mylist.pack(side=LEFT, fill='y')

        self.bar.pack   (side=LEFT, fill='y')
        self.mylist.config(yscrollcommand=self.bar.set)
        self.bar.config(command=self.mylist.yview)
        # -----------
        # cretae button
        self.btn2 = Button(self, text="Cancel", command=self.close_window)
        self.btn2.pack(side=BOTTOM, fill='x')
        # cretae button
        self.btn1 = Button(self, text="Select", command=self.selected_curve)
        self.btn1.pack(side=BOTTOM, fill='x')
        print "SONO QUIIIIIII"

    def close_window(self):
        self.destroy()
        #self.master.destroy()

    def selected_curve(self):
        #recupero la data selezionata
        curve_des = str((self.mylist.get(ACTIVE)))
        self.curve = curve_des
        self.new_window = W_boot_opt(parent=self)

class  W_boot_opt(LabelFrame):
    def sel(self):
        self.close_window()

    def close_window(self):
        self.destroy()
        self.master.destroy()

    def __init__ (self,  parent = None):
        if parent:
            self.master = parent.master
            parent.close_window()

        LabelFrame.__init__(self, self.master)

        self.config(text="Set Bootstrap Options:", width=400, height=200)

        #--- Boot swaps rates
        #self.master.grid_rowconfigure(0, weight=1)
        #self.master.grid_columnconfigure(0, weight=1)
        T1 = Label(self,height=1, width=30, text = "Swaps:").grid(row=0, sticky = "e")
        self.variable1 = StringVar(self)
        self.variable1.set("Costant Fwd Rate")  # default value
        w1 = OptionMenu(self, self.variable1, "Costant Fwd Swap Rate", "Linear Swap Rate")
        w1.grid(row=0, column=1)
        w1.config(width=30)

        T2 = Label(self, height=1, width=30, text="Futures' Gap:").grid(row=1, sticky="e")
        self.variable2 = StringVar(self)
        self.variable2.set("Previous Spot Rate")
        w2 = OptionMenu(self, self.variable2, "Previous Spot Rate", "Next Forward Rate")
        w2.grid(row=1, column=1)
        w2.config(width=30)

        T3 = Label(self, height=1, width=30, text="Futures' Convexity Adj.:").grid(row=2, sticky="e")
        self.variable3 = StringVar(self)
        self.variable3.set("No Adj.")
        w3 = OptionMenu(self, self.variable3, "No Adj", "H&W/HoLee Model")
        w3.grid(row=2, column=1)
        w3.config(width=30)

        T4 = Label(self, height=1, width=30, text="Outpur IR Capitalization:").grid(row=3, sticky="e")
        self.variable4 = StringVar(self)
        self.variable4.set("Continuous")
        w4 = OptionMenu(self, self.variable4, "Continuous", "Compounded", "Simple")
        w4.grid(row=3, column=1)
        w4.config(width=30)

        #self.pack(fill="both", expand="yes")

        B1 = Button(self, text="Select",  width=20, command=self.sel).grid(row=4, column=1, sticky ='e')
        B2 = Button(self, text="Cancel",  width=20, command=self.close_window).grid(row=5, column=1, sticky ='e')

        self.pack(fill="both", expand="yes")





#=======================================================================================================================
# punto di ingresso per load curve
# =======================================================================================================================

@xl_func
def load_swap_curve_from_db(control):
    nameSheet = "Curvette"
    xla = xl_app()
    book = xla.ActiveWorkbook
    #-----
    #creo foglio Curvette se non esiste
    try:
        s = book.Sheets(nameSheet)
    except:
        s = book.Sheets.Add()
        s.Name = nameSheet
   #------------------
    root = Tk()
    app  = W_curveType(root)
    root.mainloop()

    curve_des = app.new_window.new_window.curve
    curve_date= app.new_window.date

    cc = Curve()
    cc.ref_date = datetime.date(day=int(curve_date[-2:]), month=int(curve_date[5:7]), year=int(curve_date[:4]))
    cc.description= curve_des
    cc.loadDataFromDB()
    cc.init_finalize()

    writeCurveOnXls(cc, nameSheet, xla)


import ctypes

@xl_func
def bootstrap_from_xls(control):
    nameSheet = "Curvette"
    xla = xl_app()
    book = xla.ActiveWorkbook


    try:
        s = book.Sheets(nameSheet)
    except:
        root = Tk()
        msg = "Missing input sheet 'Curvette' in your workbook... \nNothing to do for me!"
        tkMessageBox.showinfo("Warning!", msg)
        root.destroy()
        return

    rangeStart = "B2"
    distance = 5

    curveL = readCurvesNames(xla,s,rangeStart,"o", distance)
    root = Tk()
    #root.wm_withdraw()
    W = W_bootstrapSelection(root, curveL)
    root.mainloop()

    curveDes = W.curve

    print "descrizione curva selezionata:", curveDes
    print "opzione1:", str(W.new_window.variable1.get())
    print "opzione2:",str(W.new_window.variable2.get())
    print "opzione3:", str(W.new_window.variable3.get())
    print "opzione4:", str(W.new_window.variable4.get())

    #---leggo la curva dal foglio xls







    root.destroy()

