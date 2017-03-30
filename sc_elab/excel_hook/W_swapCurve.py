from pyxll import xl_func, xl_app, xl_menu, xl_macro, xlcAlert
from xls_utils import *
from Tkinter import *
import tkMessageBox
from db_qrys import getCurvesListFromDb, getDatesListFromDb, getProvidersFromDb
from sc_elab.core.SwpCurve import *
from win32com.client import constants as const


#------
def donothing():
    tkMessageBox.showinfo("Nothing To do", "bye bye")

#-------------------
class W_curveType (Frame):
    def __init__(self, master = None):
        Frame.__init__(self, master)
        self.master = master
        self.menubar = Menu(self)
        filemenu = Menu(self.menubar, tearoff=0)

        filemenu.add_command(label="swap - RF curve", command= lambda : self.load_swp_RF_curve("SWP"))
        filemenu.add_command(label="govt curve", command=donothing)
        filemenu.add_command(label="sector curve", command=donothing)

        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=donothing)

        self.menubar.add_cascade(label="Click & Select Curve Type", menu=filemenu)
        self.master.config(menu=self.menubar)
        self.canvas = Canvas(self, bg="grey", width=200, height=200,
                             bd=0, highlightthickness=0)
        self.canvas.pack()

    def load_swp_RF_curve(self, type):
        print "type", type
        self.new_window = W_curveDate(master = None, parent = self, type = type)

    def close_window(self):
        self.destroy()


#-------------
class W_curveDate(LabelFrame):
    def __init__(self, master = None, parent = None, type = "SWP"):
        c_date = None
        if master: self.master = master
        elif parent:
            self.master = parent.master
            parent.close_window()
        # create labelframe
        LabelFrame.__init__(self, master)

        self.config(text = "Select reference date:")
        self.pack(fill="both", expand="yes")

        # create scrollbar
        self.bar = Scrollbar(self)
        # create mylist
        self.mylist = Listbox(self, yscrollcommand=self.bar.set)
        dl = getDatesListFromDb(type)
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
        self.btn1 = Button(self, text="Select", command= lambda:self.selected_date(type = type))
        self.btn1.pack(side=BOTTOM, fill='x')

    def close_window(self):
        self.destroy()

    def selected_date(self, type):
        #recupero la data selezionata
        c_date = str((self.mylist.get(ACTIVE)))

        self.date = c_date
        self.type = type
        self.new_window = W_curveSelection (parent = self, curve_date=self.date, type = type)

#----------------
class W_curveSelection (LabelFrame):
    def __init__(self, parent = None, curve_date = None, type = ""):
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
        dl = getCurvesListFromDb(curve_date, type)
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
        self.btn1 = Button(self, text="Select", command=lambda:self.selected_curve(type))
        self.btn1.pack(side=BOTTOM, fill='x')

    def close_window(self, type):
        self.destroy()
        if (type == "SWP") : self.master.destroy()

    def selected_curve(self, type):
        #recupero la data selezionata
        curve_des = str((self.mylist.get(ACTIVE)))
        self.curve = curve_des
        if type == "SWP":
            self.close_window(type)
        elif type == "CDS":
            print "XXXXXXXXXXXXXX caso CDS!!!!!!!!!"
            self.new_window = W_settoreRatingSelection(parent=self, type = type)



class W_settoreRatingSelection(LabelFrame):
    def __init__(self, parent = None, type = "CDS"):
        if parent:
            self.master = parent.master
            parent.close_window(type)
        LabelFrame.__init__(self, self.master)
        self.config(text="Provider selection:")
        # ---
        # recupero info su settore e rating
        # ---
        sectors = getProvidersFromDb("MKT_EmittenteSettore")
        ratings  = getProvidersFromDb("MKT_EmittenteRating")

        T1 = Label(self, height=1, width=30, text="Sector Provider:").grid(row=0, sticky="e")
        self.sector = StringVar(self)
        self.sector.set(sectors[0])  # default value
        w1 = OptionMenu(self, self.sector, *sectors)
        w1.grid(row=0, column=1)
        w1.config(width=30)
        #
        T2 = Label(self, height=1, width=30, text="Rating Provider:").grid(row=1, sticky="e")
        self.rating = StringVar(self)
        self.rating.set(ratings[0])
        w2 = OptionMenu(self, self.rating, *ratings)
        w2.grid(row=1, column=1)
        w2.config(width=30)

        B1 = Button(self, text="Submit", width=20, command=self.sel).grid(row=2, column=1, sticky='e')
        self.pack(fill="both", expand="yes")

        # create scrollbar
    def close_window(self):
        self.destroy()
        self.master.destroy()

    def sel(self):
        self.close_window()