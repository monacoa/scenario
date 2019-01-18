from pyxll import xl_func, xl_app, xl_menu, xl_macro, xlcAlert
from Tkinter import *
import tkMessageBox
from win32com.client import constants as const
from DEF_intef import nameSheetBootstrap

from xls_utils import readCurvesNames, readCurvesParmsNames

class W_saveType (Frame):
    def donothing(self):
        tkMessageBox.showinfo("Nothing To do!", "bye bye")
        self.master.destroy()
        return

    def __init__(self, master = None):
        Frame.__init__(self, master)
        self.master = master

        self.dataToSave = StringVar()
        self.dataToSave.set('DAT_BOOT')

        Label(self,text="""Click & Select DataType to be saved on DB :""",justify=LEFT,padx=100).pack()

        self.rb1 = Radiobutton(self,text="Bootstrapped Data (ZC &| DF)"       ,justify='left',variable=self.dataToSave,value='DAT_BOOT').pack(anchor=W)
        self.rb2 = Radiobutton(self,text="Fitting Parms (bootstrapped curves)",justify='left',variable=self.dataToSave,value='FIT_BOOT').pack(anchor=W)
        self.rb3 = Radiobutton(self,text="Fitting Parms (par-yield curves)"   ,justify='left',variable=self.dataToSave,value='FIT_PAR').pack(anchor=W)

        # create button
        self.btn2 = Button(self, text="Cancel",  command=lambda:self.close_window())
        self.btn2.pack(side=BOTTOM, fill = 'x')
        # create button
        self.btn1 = Button(self, text="Select", command= lambda:self.saveMain())
        self.btn1.pack(side=BOTTOM, fill='x')
        self.pack()

    def saveMain(self):
        if    self.dataToSave.get() == 'DAT_BOOT': self.saveBoot()
        elif  self.dataToSave.get() == 'FIT_BOOT': self.saveFittFromBoot()
        elif  self.dataToSave.get() == 'FIT_PAR' : self.saveFittFromPy()
        else: self.donothing()

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
            msg = "Missing input sheet for Bootstrapped Curves in your workbook... \nNothing to do for me!"
            tkMessageBox.showinfo("Warning!", msg)
            self.master.destroy()
            return None
        curveL = readCurvesNames(xla, s, "B2", "v", 2,  5)
        self.new_window = W_saveBootSelection(parent = self, curveL = curveL)

    def saveFittFromBoot(self):

        self.saveType = "FitFromBoot"
        nameSheet = nameSheetBootstrap
        try:
            xla = xl_app()
            book = xla.ActiveWorkbook
            s = book.Sheets(nameSheet)
            s.Activate()
        except:
            msg = "Missing input sheet for Fitted Curves in your workbook... \nNothing to do for me!"
            tkMessageBox.showinfo("Warning!", msg)
            self.master.destroy()
            return None

        curveParmsDict = readCurvesParmsNames (xla, s, "B2")

        self.new_window = W_saveFitSelection(self, curveParmsDict)

    def saveFittFromPy(self):
        self.saveType = "FitFromPy"
        self.close_window()

    def close_window(self):
        self.master.destroy()




class W_saveBootSelection(LabelFrame):

    def donothing(self):
        tkMessageBox.showinfo("Nothing To do!", "bye bye")
        self.master.destroy()
        return

    def __init__(self, parent=None, curveL=[]):
        c_date = None
        if parent:
            self.master = parent.master
            parent.destroy()

        LabelFrame.__init__(self, self.master)
        #self.master = master
        # self.geometry("800x600")
        # self.master.geometry("400x500")
        self.config(text="Available curves to save:")
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



# -----

class W_saveFitSelection(LabelFrame):

    def donothing(self):
        tkMessageBox.showinfo("Nothing To do!", "bye bye")
        self.master.destroy()
        return

    def __init__(self, parent=None, curveDict={}):
        self.curveDict = curveDict
        c_date = None
        if parent:
            self.master = parent.master
            parent.destroy()

        LabelFrame.__init__(self, self.master)
        #self.master = master
        # self.geometry("800x600")
        # self.master.geometry("400x500")
        self.config(text="Available curves to save:")
        self.pack(fill="both", expand="yes")
        # create scrollbar
        self.bar = Scrollbar(self)
        # ---
        # create curve list
        # ---
        self.mylist = Listbox(self, yscrollcommand=self.bar.set)


        for l in curveDict.keys():
            l_list = l.split("-")
            crv = l_list[0]
            pos = l_list[1]
            crv = "(" + str(pos) + ") " + crv
            self.mylist.insert(int(pos)-1, crv)

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
        self.btn1 = Button(self, text="Select", command=self.selected_curve)
        self.btn1.pack(side=BOTTOM, fill='x')


    def close_window(self):
        self.destroy()

    def selected_curve(self):
        # ---
        # recupero la data selezionata
        # ---
        curve = str((self.mylist.get(ACTIVE)))
        curve = curve.replace("(", "")
        tmp = curve.split(") ")
        self.curve = tmp[1]
        self.pos_curve = int(tmp[0])
        self.new_window = W_saveFitParms(parent=self, curve = self.curve, pos = self.pos_curve, cd = self.curveDict)


class W_saveFitParms(LabelFrame):

    def close_window(self):
        self.destroy()
        self.master.destroy()

    def donothing(self):
        tkMessageBox.showinfo("Nothing To do!", "bye bye")
        self.master.destroy()
        return

    def __init__(self, parent=None, curve = "", pos = "", cd ={}):
        if parent:
            self.master = parent.master
            parent.destroy()
        LabelFrame.__init__(self, self.master)
        self.config(text="Available parms to save:")
        self.pack(fill="both", expand="yes")
        # create scrollbar
        self.bar = Scrollbar(self)
        # ---
        # create curve list
        # ---
        self.mylist = Listbox(self, yscrollcommand=self.bar.set)
        curve_key = curve+"-"+str(pos)
        for l in cd[curve_key]:
            par     = l[0]
            pos_par = l[1]
            par = "(" + str(pos_par) + ") " + par
            self.mylist.insert(END, par)

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
        self.btn1 = Button(self, text="Select", command=self.selected)
        self.btn1.pack(side=BOTTOM, fill='x')

    def selected(self):
        stringa = str((self.mylist.get(ACTIVE)))
        stringa = stringa.replace("(", "")
        tmp = stringa.split(") ")
        self.parms = tmp[1]
        self.pos_parms = int(tmp[0])
        self.close_window()

