from pyxll import xl_func, xl_app, xl_menu, xl_macro, xlcAlert
from Tkinter import *
import tkMessageBox
from win32com.client import constants as const
from DEF_intef import nameSheetBootstrap

from xls_utils import readCurvesNames

class W_saveType (Frame):
    def donothing(self):
        tkMessageBox.showinfo("Nothing To do!", "bye bye")
        self.master.destroy()
        return

    def __init__(self, master = None):
        Frame.__init__(self, master)
        self.master = master
        self.menubar = Menu(self)
        filemenu = Menu(self.menubar, tearoff=0)
        filemenu.add_command(label="Bootstrapped Data (ZC &| DF)", command=self.saveBoot)
        filemenu.add_command(label="Fitting Parms (bootstrapped curves)", command=self.saveFittFromBoot)
        filemenu.add_command(label="Fitting Parms (par-yield curves)", command=self.saveFittFromPy)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command = self.donothing)

        self.menubar.add_cascade(label="Click & Select DataType to be saved on DB", menu=filemenu)
        self.master.config(menu=self.menubar)
        self.canvas = Canvas(self, bg="grey", width=200, height=200, bd=0, highlightthickness=0)
        self.canvas.pack()

    def saveBoot(self):
        self.saveType = "Boot"
        #apro il foglio contenente i risultati del bootstrap
        nameSheet = nameSheetBootstrap
        try:
            xla  = xl_app()
            book = xla.ActiveWorkbook
            s    = book.Sheets(nameSheet)
            s.Activate()
        except:
            msg = "Missing input sheet for Swap py Curves in your workbook... \nNothing to do for me!"
            tkMessageBox.showinfo("Warning!", msg)
            self.master.destroy()
            return None
        curveL = readCurvesNames(xla, s, "B2", "v", 2,  5)
        self.new_window = W_saveBootSelection(self, curveL)

    def saveFittFromBoot(self):
        self.saveType = "FitFromBoot"
        aaaaaaaa
        # self.fit_type = "py"
        # nameSheet = nameSheetCurve
        # xla = xl_app()
        # book = xla.ActiveWorkbook
        # try:
        #     s = book.Sheets(nameSheet)
        #     s.Activate()
        # except:
        #     msg = "Missing input sheet for Swap py Curves in your workbook... \nNothing to do for me!"
        #     tkMessageBox.showinfo("Warning!", msg)
        #     self.master.destroy()
        #     return None
        #
        # curveL = readCurvesNames(xla, s, "B2", "o", 5)
        # self.new_window = W_fittingSelection(self, curveL)

    def saveFittFromPy(self):
        self.saveType = "FitFromPy"
        zzzzzzzzzzz

    def close_window(self):
        self.destroy()




class W_saveBootSelection(LabelFrame):

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
        self.config(text="Available curves for save:")
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
        self.btn1 = Button(self, text="Select", command=self.selected_data)
        self.btn1.pack(side=BOTTOM, fill='x')


    def close_window(self):
        self.destroy()

    def selected_data(self):
        # recupero la data selezionata
        curve = str((self.mylist.get(ACTIVE)))
        curve = curve.replace("(", "")
        tmp = curve.split(") ")
        self.curve = tmp[1]

        self.pos = int(tmp[0])
        self.new_window = W_saveBoot_opt(parent=self)



class W_saveBoot_opt(LabelFrame):
    def selected(self):
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
        self.config(text="Select Data to save:", width=400, height=200)

        self.var1 = IntVar()
        c1 = Checkbutton( self, text="Discount Factors", variable=self.var1) .grid(row=0, sticky="w")
        #c1.pack()

        self.var2 = IntVar()
        c2 = Checkbutton(self, text="ZC rates", variable=self.var2).grid(row=1, sticky="w")
        #c2.pack()
        #c2.pack(side=LEFT, anchor=W, expand=YES)

        B1 = Button(self, text="Submit", width=20, command=self.selected).grid(row=2, column=1,sticky="e")
        #B1.pack()

        self.pack(fill="both", expand="yes")
        #self.pack(fill="both", expand="yes")