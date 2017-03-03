from pyxll import xl_func, xl_app, xl_menu, xl_macro, xlcAlert
from xls_utils import *
from Tkinter import *
import tkMessageBox
from db_qrys import getCurvesListFromDb, getDatesListFromDb
from sc_elab.core.SwpCurve import *
from win32com.client import constants as const

from DEF_intef import nameSheetBootstrap,nameSheetCurve

#------


class W_fittingType (Frame):
    def donothing(self):
        tkMessageBox.showinfo("Nothing To do!", "bye bye")
        self.master.destroy()
        return

    def __init__(self, master = None):
        Frame.__init__(self, master)
        self.master = master
        self.menubar = Menu(self)
        filemenu = Menu(self.menubar, tearoff=0)
        filemenu.add_command(label="Fit Bootstrapped Curve", command=self.fit_boot_curve)
        filemenu.add_command(label="Fit par-yield Swap Curve", command=self.fit_py_curve)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command = self.donothing)

        self.menubar.add_cascade(label="Click & Select Fitting Type", menu=filemenu)
        self.master.config(menu=self.menubar)
        self.canvas = Canvas(self, bg="grey", width=200, height=200, bd=0, highlightthickness=0)
        self.canvas.pack()

    def fit_boot_curve(self):
        self.fit_type = "boot"
        nameSheet = nameSheetBootstrap

        xla = xl_app()
        book = xla.ActiveWorkbook
        try:
            s = book.Sheets(nameSheet)
            s.Activate()
        except:
            msg ="Missing input sheet for Bootstrapped Curves in your workbook... \nNothing to do for me!"
            tkMessageBox.showinfo("Warning!", msg)
            self.master.destroy()
            return None

        curveL = readCurvesNames(xla, s, "B2", "v", 2)
        self.new_window = W_fittingSelection(self, curveL)


    def fit_py_curve(self):
        self.fit_type = "py"
        nameSheet = nameSheetCurve
        xla = xl_app()
        book = xla.ActiveWorkbook
        try:
            s = book.Sheets(nameSheet)
            s.Activate()
        except:
            msg = "Missing input sheet for Swap py Curves in your workbook... \nNothing to do for me!"
            tkMessageBox.showinfo("Warning!", msg)
            self.master.destroy()
            return None

        curveL = readCurvesNames(xla, s, "B2", "o", 5)
        self.new_window = W_fittingSelection(self, curveL)






    def close_window(self):
        self.destroy()


class W_fittingSelection(LabelFrame):

    def donothing(self):
        tkMessageBox.showinfo("Nothing To do!", "bye bye")
        self.master.destroy()
        return

    def __init__(self, master=None, curveL=[]):
        c_date = None
        if master:
            self.master = master.master
            master.close_window()

        LabelFrame.__init__(self, self.master)
        #self.master = master
        # self.geometry("800x600")
        # self.master.geometry("400x500")
        self.config(text="Available bootstraped curves for fitting:")
        self.pack(fill="both", expand="yes")
        # create scrollbar
        self.bar = Scrollbar(self)

        # create mylist
        self.mylist = Listbox(self, yscrollcommand=self.bar.set)

        for l in curveL:
            crv = l[0]
            pos = l[1]
            crv = "(" + str(pos) + ") " + crv
            self.mylist.insert(END, crv)

        self.mylist.config(width=50)
        self.mylist.pack(side=LEFT, fill='y')
        # ----------.
        self.bar.pack(side=LEFT, fill='y')
        self.mylist.config(yscrollcommand=self.bar.set)
        self.bar.config(command=self.mylist.yview)
        # -----------
        # cretae button
        self.btn2 = Button(self, text="Cancel", command=self.donothing)
        self.btn2.pack(side=BOTTOM, fill='x')
        # cretae button
        self.btn1 = Button(self, text="Select", command=self.selected_fit_curve)
        self.btn1.pack(side=BOTTOM, fill='x')

    def close_window(self):
        self.destroy()

    def selected_fit_curve(self):
        # recupero la data selezionata
        curve = str((self.mylist.get(ACTIVE)))
        curve = curve.replace("(", "")
        tmp = curve.split(") ")
        self.curve = tmp[1]

        self.pos = int(tmp[0])
        self.new_window = W_fit_opt(parent=self)

class W_fit_opt(LabelFrame):
    def sel(self):
        self.close_window()

    def close_window(self):
        self.destroy()
        self.master.destroy()


    def donothing(self):
        self.destroy()
        self.master.destroy()

    def __init__(self, parent=None):
        if parent:
            self.master = parent.master
            parent.close_window()

        LabelFrame.__init__(self, self.master)

        self.config(text="Set Fitting Options:", width=400, height=200)

        # --- Interpolation Type

        T1 = Label(self, height=1, width=30, text="Interpolation Type:").grid(row=0, sticky="e")
        self.variable1 = StringVar(self)
        self.variable1.set("(0) Linear")  # default value
        w1 = OptionMenu(self, self.variable1, "(0) Linear", "(1) Adam Van Deventer", "(2) Svensson")
        w1.grid(row=0, column=1)
        w1.config(width=30)

        T2 = Label(self, height=1, width=30, text="Tenor Fwd Rates (Graphs):").grid(row=1, sticky="e")
        self.variable2 = StringVar(self)
        self.variable2.set("(1) 1M")
        w2 = OptionMenu(self, self.variable2, "1M", "3M", "6M", "12M")
        w2.grid(row=1, column=1)
        w2.config(width=30)

        T3 = Label(self, height=1, width=30, text="Graphs location (save):").grid(row=2, sticky="e")
        self.variable5 = StringVar(self)
        self.variable5.set("C://")
        w5 = Entry(self, textvariable=self.variable5)
        w5.grid(row=2, column=1)
        w5.config(width=30)

        self.pack(fill="both", expand="yes")

        # questa istruzione va messa dopo il pack altrimenti non vienne intercettato il valore dal .get() successivo
        self.variable5.set("C://")

        B1 = Button(self, text="Submit", width=20, command=self.sel).grid(row=3, column=1, sticky='e')


        self.pack(fill="both", expand="yes")