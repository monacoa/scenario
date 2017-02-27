from pyxll import xl_func, xl_app, xl_menu, xl_macro, xlcAlert
from xls_utils import *
from Tkinter import *
import tkMessageBox
from db_qrys import getCurvesListFromDb, getDatesListFromDb
from sc_elab.core.SwpCurve import *
from win32com.client import constants as const


#------
def donothing():
    tkMessageBox.showinfo("Nothing To do", "bye bye")
    # filewin = Toplevel(root)
    # button = Button(filewin, text="Do nothing button")
    # button.pack()
    # filewin.mainloop()

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
        s.Activate()
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

