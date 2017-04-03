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


class W_bootstrapSelection (LabelFrame):
    def __init__(self, master = None, curveL = [], type = "SWP"):
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
        self.btn2 = Button(self, text="Cancel", command=self.donothing)
        self.btn2.pack(side=BOTTOM, fill='x')
        # cretae button
        self.btn1 = Button(self, text="Select", command=lambda: self.selected_curve(type))
        self.btn1.pack(side=BOTTOM, fill='x')

    def close_window(self):
        self.destroy()
        #self.master.destroy()

    def donothing(self):
        self.destroy()
        self.master.destroy()
    def selected_curve(self, type):
        #recupero la curva selezionata
        curve = str((self.mylist.get(ACTIVE)))
        curve  = curve.replace("(", "")
        tmp = curve.split(") ")
        self.curve = tmp[1]
        self.pos = int(tmp[0])
        self.new_window = W_boot_opt(parent=self, type = type)


class  W_boot_opt(LabelFrame):

    def sel(self):
        self.donothing()


    def donothing(self):
        self.destroy()
        self.master.destroy()
        return

    def __init__ (self,  parent = None, type = "SWP"):
        if parent:
            self.master = parent.master
            parent.close_window()

        LabelFrame.__init__(self, self.master)

        self.config(text="Set Bootstrap Options:", width=400, height=200)
        #--- Boot swaps rates
        T1 = Label(self,height=1, width=30, text = "Bootstrap method (filling):").grid(row=0, sticky = "e")
        self.variable1 = StringVar(self)
        self.variable1.set("(0) Costant Fwd Rate")  # default value
        w1 = OptionMenu(self, self.variable1, "(0) Costant Fwd py Rate", "(1) Linear py Rate")
        w1.grid(row=0, column=1)
        w1.config(width=30)

        if type == "SWP":
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

        elif type == "CDS":
            T6 = Label(self, height=1, width=30, text="Interpolation Type:").grid(row=1, sticky="e")
            self.variable6 = StringVar(self)
            self.variable6.set("(0) Linear")  # default value
            w6 = OptionMenu(self, self.variable6, "(0) Linear", "(1) SVE")
            w6.grid(row=1, column=1)
            w6.config(width=30)
            # ---
            T7 = Label(self, height=1, width=30, text="Bootstrap Method for HR:").grid(row=2, sticky="e")
            self.variable7 = StringVar(self)
            self.variable7.set("(0) Constant Hazard Rate")  # default value
            w7 = OptionMenu(self, self.variable7, "(0) Constant Hazard Rate", "(1) Linear Credit spread")
            w7.grid(row=2, column=1)
            w7.config(width=30)
            # ---

        self.pack(fill="both", expand="yes")

        #questa istruzione va messa dopo il pack altrimenti non vienne intercettato il valore dal .get() successivo
        if type == "SWP": self.variable5.set("C://")

        B1 = Button(self, text="Submit",  width=20, command=self.sel).grid(row=5, column=1, sticky ='e')

        self.pack(fill="both", expand="yes")