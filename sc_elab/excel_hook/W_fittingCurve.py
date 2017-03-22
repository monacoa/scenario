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

import ttk
class W_fit_opt(LabelFrame):
#class W_fit_opt(ttk.Notebook):

    def sel(self):
        self.close_window()

    def close_window(self):
        self.destroy()
        self.master.destroy()


    def donothing(self):
        self.destroy()
        self.master.destroy()

    # def __init__(self, parent=None):
    #     if parent:
    #         self.master = parent.master
    #         parent.close_window()
    #
    #     LabelFrame.__init__(self, self.master)
    #
    #     self.config(text="Set Fitting Options:", width=400, height=200)
    #
    #     # --- Interpolation Type
    #
    #     T1 = Label(self, height=1, width=30, text="Interpolation Type:").grid(row=0, sticky="e")
    #     self.variable1 = StringVar(self)
    #     self.variable1.set("(0) Linear")  # default value
    #     w1 = OptionMenu(self, self.variable1, "(0) Linear", "(1) Adam Van Deventer", "(2) Svensson")
    #     w1.grid(row=0, column=1)
    #     w1.config(width=30)
    #
    #     T2 = Label(self, height=1, width=30, text="Tenor Fwd Rates (Graphs):").grid(row=1, sticky="e")
    #     self.variable2 = StringVar(self)
    #     self.variable2.set("(1) 1M")
    #     w2 = OptionMenu(self, self.variable2, "1M", "3M", "6M", "12M")
    #     w2.grid(row=1, column=1)
    #     w2.config(width=30)
    #
    #     T3 = Label(self, height=1, width=30, text="Graphs location (save):").grid(row=2, sticky="e")
    #     self.variable5 = StringVar(self)
    #     self.variable5.set("C://")
    #     w5 = Entry(self, textvariable=self.variable5)
    #     w5.grid(row=2, column=1)
    #     w5.config(width=30)
    #
    #     self.pack(fill="both", expand="yes")
    #
    #     # questa istruzione va messa dopo il pack altrimenti non vienne intercettato il valore dal .get() successivo
    #     self.variable5.set("C://")
    #
    #     B1 = Button(self, text="Submit", width=20, command=self.sel).grid(row=3, column=1, sticky='e')
    #
    #
    #     self.pack(fill="both", expand="yes")



    '''
    def __init__(self, parent=None):
        if parent:
            self.master = parent.master
            parent.close_window()


        ttk.Notebook.__init__(self, self.master)
        #self.pack(fill=BOTH, padx=2, pady=3)

        master_frame = Frame(self, name='master_frame')


        self.add(master_frame, text="OPTIONS")

        T1 = Label(master_frame, height=1, width=30, text="Interpolation Type:").grid(row=0, sticky="e")
        self.variable1 = StringVar(self)
        self.variable1.set("(0) Linear")  # default value
        w1 = OptionMenu(master_frame, self.variable1, "(0) Linear", "(1) Adam Van Deventer", "(2) Svensson")
        w1.grid(row=0, column=1)
        w1.config(width=30)

        T2 = Label(master_frame, height=1, width=30, text="Tenor Fwd Rates (Graphs):").grid(row=1, sticky="e")
        self.variable2 = StringVar(self)
        self.variable2.set("(1) 1M")
        w2 = OptionMenu(master_frame, self.variable2, "1M", "3M", "6M", "12M")
        w2.grid(row=1, column=1)
        w2.config(width=30)

        T3 = Label(master_frame, height=1, width=30, text="Graphs location (save):").grid(row=2, sticky="e")
        self.variable5 = StringVar(self)
        #self.variable5.set("C://")
        w5 = Entry(master_frame, textvariable=self.variable5)
        w5.grid(row=2, column=1)
        w5.config(width=30)

        # print "impacchetto"

        param_frame = Frame(self, name='param_frame')
        self.add(param_frame, text="PARAMS")

        def createInsertRange(var_min, var_max, row_nr, frame, par_name):
            p1 = Entry(frame, textvariable=var_min)
            p1.grid(row=row_nr, column=0)
            p1.config(width=30)
            # -
            p2 = Label(frame, height=1, width=30, text=par_name).grid(row=row_nr, column=1)
            # -
            p3 = Entry(frame, textvariable=var_max)
            p3.grid(row=row_nr, column=2)
            p3.config(width=30)

        #------- CIR SECTION
        tl1 = Label(param_frame, height=2, width=30, text="Modello CIR").grid(row=0, column=0)

        i1 = Label(param_frame, height=1, width=30, text="min").grid(row=1, column=0)
        i2 = Label(param_frame, height=1, width=30, text="param").grid(row=1, column=1)
        i3 = Label(param_frame, height=1, width=30, text="max").grid(row=1, column=2)
        #---
        self.CIR_r0min = StringVar(self)
        self.CIR_r0max = StringVar(self)
        createInsertRange(self.CIR_r0min, self.CIR_r0max, 2, param_frame, "r(0)")
        self.CIR_kmin = StringVar(self)
        self.CIR_kmax = StringVar(self)
        createInsertRange(self.CIR_kmin, self.CIR_kmax, 3, param_frame, "kappa")
        self.CIR_thetamin = StringVar(self)
        self.CIR_thetamax = StringVar(self)
        createInsertRange(self.CIR_thetamin, self.CIR_thetamax, 4, param_frame, "theta")
        self.CIR_sigmamin = StringVar(self)
        self.CIR_sigmamax = StringVar(self)
        createInsertRange(self.CIR_sigmamin, self.CIR_sigmamax, 5, param_frame, "sigma")


        # ------- SVD SECTION
        tl1 = Label(param_frame, height=2, width=30, text="Modello SVE").grid(row=7, column=0)

        i1 = Label(param_frame, height=1, width=30, text="min").grid(row=9, column=0)
        i2 = Label(param_frame, height=1, width=30, text="param").grid(row=9, column=1)
        i3 = Label(param_frame, height=1, width=30, text="max").grid(row=9, column=2)
        # ---
        self.SVE_b0min = StringVar(self)
        self.SVE_b0max = StringVar(self)
        createInsertRange(self.SVE_b0min, self.SVE_b0max, 10, param_frame, "beta0")
        self.SVE_b1min = StringVar(self)
        self.SVE_b1max = StringVar(self)
        createInsertRange(self.SVE_b1min, self.SVE_b1max, 11, param_frame, "beta1")
        self.SVE_b2min = StringVar(self)
        self.SVE_b2max = StringVar(self)
        createInsertRange(self.SVE_b2min, self.SVE_b2max, 12, param_frame, "beta2")
        self.SVE_b3min = StringVar(self)
        self.SVE_b3max = StringVar(self)
        createInsertRange(self.SVE_b3min, self.SVE_b3max, 13, param_frame, "beta3")
        self.SVE_c1min = StringVar(self)
        self.SVE_c1max = StringVar(self)
        createInsertRange(self.SVE_c1min, self.SVE_c1max, 14, param_frame, "const1")
        self.SVE_c2min = StringVar(self)
        self.SVE_c2max = StringVar(self)
        createInsertRange(self.SVE_c2min, self.SVE_c2max, 15, param_frame, "const2")

        self.pack(fill=BOTH, padx=2, pady=3)
        #-- set default values
        self.SVE_b0min.set("-0.05")
        self.SVE_b0max.set("10.0")
        self.SVE_b1min.set("-0.05")
        self.SVE_b1max.set("10.0")
        self.SVE_b2min.set("-0.05")
        self.SVE_b2max.set("10.0")
        self.SVE_b3min.set("-0.05")
        self.SVE_b3max.set("10.0")
        self.SVE_c1min.set("0.001")
        self.SVE_c1max.set("100.0")
        self.SVE_c2min.set("0.001")
        self.SVE_c2max.set("100.0")
        #--
        self.CIR_r0min.set("0.0")
        self.CIR_r0max.set("0.02")
        self.CIR_kmin.set("0.001")
        self.CIR_kmax.set("10.0")
        self.CIR_sigmamin.set("0.001")
        self.CIR_sigmamax.set("0.50")
        self.CIR_thetamin.set("0.005")
        self.CIR_thetamax.set("1.0")
        # --- Interpolation Type
    '''

    def __init__(self, parent=None):
        if parent:
            self.master = parent.master
            parent.close_window()

        LabelFrame.__init__(self, self.master)
        self.config(text="Set Fitting Options:", width=400, height=200)
        self.pack(fill=BOTH)

        self.nb = ttk.Notebook(self)
        master_frame = Frame(self.nb, name='master_frame')
        self.nb.add(master_frame, text="OPTIONS")

        T1 = Label(master_frame, height=1, width=30, text="Interpolation Type:").grid(row=0, sticky="e")
        self.variable1 = StringVar(self)
        self.variable1.set("(0) Linear")  # default value
        w1 = OptionMenu(master_frame, self.variable1, "(0) Linear", "(1) Adam Van Deventer", "(2) Svensson", "(3) CIR")
        w1.grid(row=0, column=1)
        w1.config(width=30)

        T2 = Label(master_frame, height=1, width=30, text="Tenor Fwd Rates (Graphs):").grid(row=1, sticky="e")
        self.variable2 = StringVar(self)
        self.variable2.set("(1) 1M")
        w2 = OptionMenu(master_frame, self.variable2, "1M", "3M", "6M", "12M")
        w2.grid(row=1, column=1)
        w2.config(width=30)

        T3 = Label(master_frame, height=1, width=30, text="Graphs location (save):").grid(row=2, sticky="e")
        self.variable5 = StringVar(self)
        #
        w5 = Entry(master_frame, textvariable=self.variable5)
        w5.grid(row=2, column=1)
        w5.config(width=30)

        # print "impacchetto"

        param_frame = Frame(self, name='param_frame')
        self.nb.add(param_frame, text="PARAMS")

        def createInsertRange(var_min, var_max, row_nr, frame, par_name):
            p1 = Entry(frame, textvariable=var_min)
            p1.grid(row=row_nr, column=0)
            p1.config(width=30)
            # -
            p2 = Label(frame, height=1, width=30, text=par_name).grid(row=row_nr, column=1)
            # -
            p3 = Entry(frame, textvariable=var_max)
            p3.grid(row=row_nr, column=2)
            p3.config(width=30)

        # ------- CIR SECTION
        tl1 = Label(param_frame, height=2, width=30, text="Modello CIR").grid(row=0, column=0)

        i1 = Label(param_frame, height=1, width=30, text="min").grid(row=1, column=0)
        i2 = Label(param_frame, height=1, width=30, text="param").grid(row=1, column=1)
        i3 = Label(param_frame, height=1, width=30, text="max").grid(row=1, column=2)
        # ---
        self.CIR_r0min = StringVar(self)
        self.CIR_r0max = StringVar(self)
        createInsertRange(self.CIR_r0min, self.CIR_r0max, 2, param_frame, "r(0)")
        self.CIR_kmin = StringVar(self)
        self.CIR_kmax = StringVar(self)
        createInsertRange(self.CIR_kmin, self.CIR_kmax, 3, param_frame, "kappa")
        self.CIR_thetamin = StringVar(self)
        self.CIR_thetamax = StringVar(self)
        createInsertRange(self.CIR_thetamin, self.CIR_thetamax, 4, param_frame, "theta")
        self.CIR_sigmamin = StringVar(self)
        self.CIR_sigmamax = StringVar(self)
        createInsertRange(self.CIR_sigmamin, self.CIR_sigmamax, 5, param_frame, "sigma")

        # ------- SVD SECTION
        tl1 = Label(param_frame, height=2, width=30, text="Modello SVE").grid(row=7, column=0)

        i1 = Label(param_frame, height=1, width=30, text="min").grid(row=9, column=0)
        i2 = Label(param_frame, height=1, width=30, text="param").grid(row=9, column=1)
        i3 = Label(param_frame, height=1, width=30, text="max").grid(row=9, column=2)
        # ---
        self.SVE_b0min = StringVar(self)
        self.SVE_b0max = StringVar(self)
        createInsertRange(self.SVE_b0min, self.SVE_b0max, 10, param_frame, "beta0")
        self.SVE_b1min = StringVar(self)
        self.SVE_b1max = StringVar(self)
        createInsertRange(self.SVE_b1min, self.SVE_b1max, 11, param_frame, "beta1")
        self.SVE_b2min = StringVar(self)
        self.SVE_b2max = StringVar(self)
        createInsertRange(self.SVE_b2min, self.SVE_b2max, 12, param_frame, "beta2")
        self.SVE_b3min = StringVar(self)
        self.SVE_b3max = StringVar(self)
        createInsertRange(self.SVE_b3min, self.SVE_b3max, 13, param_frame, "beta3")
        self.SVE_c1min = StringVar(self)
        self.SVE_c1max = StringVar(self)
        createInsertRange(self.SVE_c1min, self.SVE_c1max, 14, param_frame, "const1")
        self.SVE_c2min = StringVar(self)
        self.SVE_c2max = StringVar(self)
        createInsertRange(self.SVE_c2min, self.SVE_c2max, 15, param_frame, "const2")

        self.nb.pack(fill=BOTH, padx=2, pady=3)

        # -- set default values

        self.variable5.set("C://")
        self.SVE_b0min.set("-0.05")
        self.SVE_b0max.set("10.0")
        self.SVE_b1min.set("-0.05")
        self.SVE_b1max.set("10.0")
        self.SVE_b2min.set("-0.05")
        self.SVE_b2max.set("10.0")
        self.SVE_b3min.set("-0.05")
        self.SVE_b3max.set("10.0")
        self.SVE_c1min.set("0.001")
        self.SVE_c1max.set("100.0")
        self.SVE_c2min.set("0.001")
        self.SVE_c2max.set("100.0")
        # --
        self.CIR_r0min.set("0.0")
        self.CIR_r0max.set("0.02")
        self.CIR_kmin.set("0.001")
        self.CIR_kmax.set("10.0")
        self.CIR_sigmamin.set("0.001")
        self.CIR_sigmamax.set("0.50")
        self.CIR_thetamin.set("0.005")
        self.CIR_thetamax.set("1.0")
        # --- Final Button
        button = Button(self, text='Submit',command=self.sel)
        button.pack()
        # ---