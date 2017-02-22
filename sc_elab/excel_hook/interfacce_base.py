
from pyxll import xl_func, xl_app, xl_menu, xl_macro, xlcAlert
from connection import *
from Tkinter import *
import tkMessageBox
from db_qrys import getCurvesListFromDb
from sc_elab.core.SwpCurve import *
from win32com.client import constants as const

FORMAT = "dd/mm/yyyy"

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

    if dir == "v" :
        if (r.Value == None): return r
        nCols = r.Columns.Count
        row   = r.Row
        col   = r.Column
        j = 0
        while (r.Value != None):
            j += 1
            r = xla.Range(xla.Cells(row + j, col), xla.Cells(row + j, col))
            #porto avanti ancora per controllare
            if (r.Value == None):  r = xla.Range(xla.Cells(row + j + distCurve , col), xla.Cells(row + j+distCurve, col))
        r = xla.Range(xla.Cells(row + j + distCurve, col), xla.Cells(row+j+distCurve, col))
    else:
        while (r.Value != None):

            nCols = r.Columns.Count
            row =r.Row
            col =r.Column
            r = xla.Range(xla.Cells(row, col + distCurve), xla.Cells(row, col + (nCols-1) + distCurve))

    rOut = r
    #-----
    if (rOut == None):
        msg = "Unable to compute the output range for your curve"
        popup_messagebox(msg)
        sys.exit()

    return rOut





def drawBox(xla, spessore , rTopLeft = 0, cTopLeft = 0,rBottomRight=0,cBottomRight=0, Colore=0):

    """

    :rtype: object
    """
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

        if (type(attributi[k]) == datetime.datetime) or (type(attributi[k]) == datetime.date):
            b.NumberFormat = FORMAT
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
                if (type(ll[j]) == datetime.date) or (type(ll[j]) == datetime.datetime): a.NumberFormat = FORMAT
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
                rangeStartNew = segmentoSwapCurve (xla, rangeStartNew, code, crv.segms[s])



def writeBootstrapResOnXls (crv,xla, str_boot_opt, res):
    nameSheet = "BootstrapSwapCurve"

    try:
        s = xla.ActiveWorkbook.Sheets(nameSheet)
    except:
        s = xla.ActiveWorkbook.Sheets.Add()
        s.Name = nameSheet

    s.Activate()
    rangeStart = "B2"
    distCurve  = 1
    r    = s.Range(rangeStart)
    r    = findRigthPlaceBootCurveSeg( xla, r, distCurve, "v")
    str_segms = crv.getStrSegms()
    #-----
    Attributi_1 = { "Date Ref"     : crv.ref_date
                  , "Description"  : crv.description
                  , "Currency"     : crv.curr
                  , "Download Type": crv.download_type
                  , "Quotation"    : crv.quotation
                  , "Source"       : crv.source
                  , "CurveType"    : crv.type
                  , "Return"       : "Zero Coupon"
                  , "Node Type"    : "Spot"
                  , "Boot Optons"  : str_boot_opt
                  , "Segms"        : str_segms
                    }

    Attributi_2 =   { "Date Ref": crv.ref_date
                    , "Description": crv.description
                    , "Currency": crv.curr
                    , "Download Type": crv.download_type
                    , "Quotation": crv.quotation
                    , "Source": crv.source
                    , "CurveType": crv.type
                    , "Return": "-"
                    , "Node Type": "Discount Factors"
                    , "Boot Optons": str_boot_opt
                    , "Segms": str_segms
                    }


    r2 = s.Range(xla.Cells(r.Row, r.Column+3),xla.Cells(r.Row, r.Column+3))

    ra = intestazioneSwapCurveSegmenti(xla, s, r, Attributi_1, nCols = 2 )
    rb = intestazioneSwapCurveSegmenti(xla, s, r2, Attributi_2, nCols = 2 )
    r=s.Range(ra)

    topLeftRow = r.Row
    topLeftCol = r.Column

    xla.Cells(topLeftRow-1, topLeftCol  ).Value               = "Date"
    xla.Cells(topLeftRow-1, topLeftCol  ).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow-1, topLeftCol+1).Value               = "Value"
    xla.Cells(topLeftRow-1, topLeftCol+1).HorizontalAlignment = const.xlCenter

    xla.Cells(topLeftRow - 1, topLeftCol+3).Value = "Date"
    xla.Cells(topLeftRow - 1, topLeftCol+3).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 4).Value = "Value"
    xla.Cells(topLeftRow - 1, topLeftCol + 4).HorizontalAlignment = const.xlCenter
    xla.Cells(topLeftRow - 1, topLeftCol + 5).Value = "Usage"
    xla.Cells(topLeftRow - 1, topLeftCol + 5).HorizontalAlignment = const.xlCenter

    nRows = len(res['DiscountFactors'])

    drawBox (xla, const.xlMedium, topLeftRow-1, topLeftCol  , topLeftRow+nRows-1, topLeftCol+1)
    drawLine(xla, topLeftRow-1  , topLeftCol  , topLeftRow-1, topLeftCol + 1    , "o", const.xlThin)

    drawBox(xla, const.xlMedium, topLeftRow - 1, topLeftCol+3, topLeftRow + nRows - 1, topLeftCol + 5)
    drawLine(xla, topLeftRow - 1, topLeftCol+3, topLeftRow - 1, topLeftCol + 5, "o", const.xlThin)
    for i in range(nRows):
        date  = res['DateScadenza'][i]
        value = res['DiscountFactors'][i]
        rate  = res ['Nodi']
        xla.Cells(topLeftRow + i, topLeftCol).Value = date
        xla.Cells(topLeftRow + i, topLeftCol).NumberFormat = FORMAT
        xla.Cells(topLeftRow + i, topLeftCol).HorizontalAlignment = const.xlCenter
        xla.Cells(topLeftRow + i, topLeftCol + 1).Value = value
        xla.Cells(topLeftRow + i, topLeftCol + 1).NumberFormat = "0.000000"
        xla.Cells(topLeftRow + i, topLeftCol + 1).HorizontalAlignment = const.xlCenter

        xla.Cells(topLeftRow + i, topLeftCol+3).Value = date
        xla.Cells(topLeftRow + i, topLeftCol+3).NumberFormat = FORMAT
        xla.Cells(topLeftRow + i, topLeftCol+3).HorizontalAlignment = const.xlCenter
        xla.Cells(topLeftRow + i, topLeftCol + 4).Value = rate
        xla.Cells(topLeftRow + i, topLeftCol + 4).NumberFormat = "0.00"
        xla.Cells(topLeftRow + i, topLeftCol + 4).HorizontalAlignment = const.xlCenter
        xla.Cells(topLeftRow + i, topLeftCol + 5).Value = "Y"
        xla.Cells(topLeftRow + i, topLeftCol + 5).HorizontalAlignment = const.xlCenter


def readIntestazione(xla , r , cc):
    row     = r.Row
    col     = r.Column
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
        elif name == "1. Depositi"      : code = "D"
        elif name == "2. Libor"         : code = "L"
        elif name == "3. Futures"       : code = "F"
        elif name == "4. Swap Rate"     : code = "S"
        else                            : code = name

        if code == 'G':
            add = cc.floater_tenor[0]
            print add
            code += add
            print "newCODE:::::::::::::::::::", code

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



def readCurvesNames(xla, s, rangeStart, direzione, distanza):
    r = xla.Range(rangeStart)

    if direzione.lower() == "o":
        curveL = []
        i = 0
        while (r.Value != None):
            i += 1

            nomeCurva = r.Value
            curveL.append((nomeCurva,i))

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
            crv = l[0]
            pos = l[1]
            crv = "("+str(pos)+") "+ crv
            self.mylist.insert(END, crv)
        self.mylist.config(width=50)
        self.mylist.pack(side=LEFT, fill='y')
        # ----------.
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
        #self.master.destroy()

    def selected_curve(self):
        #recupero la data selezionata
        curve = str((self.mylist.get(ACTIVE)))
        curve  = curve.replace("(", "")
        tmp = curve.split(") ")
        self.curve = tmp[1]

        self.pos = int(tmp[0])
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

        T1 = Label(self,height=1, width=30, text = "Filling Swaps:").grid(row=0, sticky = "e")
        self.variable1 = StringVar(self)
        self.variable1.set("(0) Costant Fwd Rate")  # default value
        w1 = OptionMenu(self, self.variable1, "(0) Costant Fwd Swap Rate", "(1) Linear Swap Rate")
        w1.grid(row=0, column=1)
        w1.config(width=30)

        T2 = Label(self, height=1, width=30, text="Futures' Gap:").grid(row=1, sticky="e")
        self.variable2 = StringVar(self)
        self.variable2.set("(1) Previous Spot Rate")
        w2 = OptionMenu(self, self.variable2, "(1) Previous Spot Rate", "(0) Next Forward Rate")
        w2.grid(row=1, column=1)
        w2.config(width=30)

        T3 = Label(self, height=1, width=30, text="Futures' Convexity Adj.:").grid(row=2, sticky="e")
        self.variable3 = StringVar(self)
        self.variable3.set("(1) H&W/HoLee Model")
        w3 = OptionMenu(self, self.variable3, "(0) No Adj", "(1) H&W/HoLee Model")
        w3.grid(row=2, column=1)
        w3.config(width=30)

        T4 = Label(self, height=1, width=30, text="Output IR Capitalization:").grid(row=3, sticky="e")
        self.variable4 = StringVar(self)
        self.variable4.set("(2) Continuos")
        w4 = OptionMenu(self, self.variable4, "(0) Simple", "(1) Compounded", "(2) Continuous")
        w4.grid(row=3, column=1)
        w4.config(width=30)

        T5 = Label(self, height=1, width=30, text="Graphs location (save):").grid(row=4, sticky="e")
        self.variable5 = StringVar(self)
        self.variable5.set("C://")
        w5 = Entry(self, textvariable=self.variable5)
        w5.grid(row=4, column=1)
        w5.config(width=30)

        self.pack(fill="both", expand="yes")

        #questa istruzione va messa dopo il pack altrimenti non vienne intercettato il valore dal .get() successivo
        self.variable5.set("C://")

        B1 = Button(self, text="Select",  width=20, command=self.sel).grid(row=5, column=1, sticky ='e')
        B2 = Button(self, text="Cancel",  width=20, command=self.close_window).grid(row=6, column=1, sticky ='e')

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
        s.Activate()
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
    curvePos = W.pos

    #opt
    opt_swaps     = (str(W.new_window.variable1.get()).strip(""))[1]
    opt_fut_gap   = (str(W.new_window.variable2.get()).strip(""))[1]
    opt_conv_adj  = (str(W.new_window.variable3.get()).strip(""))[1]
    opt_out_cap   = (str(W.new_window.variable4.get()).strip(""))[1]
    opt_path_graph=  W.new_window.variable5.get()

    str_boot_opt = opt_swaps+","+opt_fut_gap + "," + opt_conv_adj + "," + opt_out_cap

    data_opt                    = {}

    data_opt['GapFutures']      = opt_fut_gap
    data_opt['SwapGapMethod']   = opt_swaps
    data_opt['RegimeOutput']    = opt_out_cap
    data_opt['Convexity']       = opt_conv_adj
    data_opt['MakeGraph']       = True
    data_opt['SaveGraph']       = True
    data_opt['FutureTenor']     = 90
    data_opt['Path']            = opt_path_graph

    curve    = readCurveFromXls(xla, curveDes, curvePos, nameSheet)
    boot_out = curve.bootstrap(data_opt)

    writeBootstrapResOnXls(curve, xla, str_boot_opt,boot_out)

@xl_func
def fitting_from_xls(control):
    nameSheet = "FittingSwapCurve"
    xla = xl_app()
    book = xla.ActiveWorkbook
    try:
        s = book.Sheets(nameSheet)
        s.Activate()
    except:
        root = Tk()
        msg = "Missing sheet  'FittingSwapCurve' in your workbook... \nNothing to do for me!"
        tkMessageBox.showinfo("Warning!", msg)
        root.destroy()
        return

