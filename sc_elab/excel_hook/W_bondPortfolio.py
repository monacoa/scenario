from pyxll import xl_func, xl_app, xl_menu, xl_macro, xlcAlert
from xls_utils import *
from Tkinter import *
import tkMessageBox
from db_qrys import getCurvesListFromDb, getDatesListFromDb, getProvidersFromDb, getBondListFromDb
from sc_elab.core.SwpCurve import *
import ttk

from win32com.client import constants as const


#------
def donothing():
    tkMessageBox.showinfo("Nothing To do", "bye bye")

#-------------------
class W_bondType (Frame):
    def __init__(self, master = None):
        Frame.__init__(self, master)
        self.master = master
        self.menubar = Menu(self)
        filemenu = Menu(self.menubar, tearoff=0)

        filemenu.add_command(label="swap - RF curve", command= lambda : self.load_swp_RF_curve("BOND"))
        filemenu.add_command(label="govt curve", command=donothing)
        filemenu.add_command(label="sector curve", command=donothing)

        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=donothing)

        self.menubar.add_cascade(label="Click & Select Curve Type", menu=filemenu)
        self.master.config(menu=self.menubar)
        self.canvas = Canvas(self, bg="grey", width=200, height=200,
                             bd=0, highlightthickness=0)
        self.canvas.pack()

    def load_bondData(self, type):
        print "type", type
        self.new_window = W_bondDate(master = None, parent = self, type = type)

    def close_window(self):
        self.destroy()


class W_bondDate(LabelFrame):
    def __init__(self, master = None, parent = None, type = "BOND"):
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
        self.new_window = W_bondSelection (parent = self, bond_date=self.date, type = type)


class W_bondSelection (LabelFrame):
    def __init__(self, parent = None, bond_date = None, type = ""):
        c_date = None
        if parent:
            self.master = parent.master
            parent.close_window()
        LabelFrame.__init__(self, self.master)
        #self.geometry("800x600")
        #self.master.geometry("400x500")
        self.config(text="Available bonds for the selected date:")
        self.pack(fill="both", expand="yes")
        # create scrollbar
        self.bar = Scrollbar(self)

        # create mylist
        self.mylist = Listbox(self, yscrollcommand=self.bar.set)
        dl = getBondListFromDb(bond_date, type)
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
        self.btn1 = Button(self, text="Select", command=lambda:self.selected_bond(type))
        self.btn1.pack(side=BOTTOM, fill='x')

    def close_window(self, type):
        self.destroy()
        if (type == "BOND") : self.master.destroy()

    def selected_bond(self, type):
        #recupero la data selezionata
        curve_des = str((self.mylist.get(ACTIVE)))
        self.bond = curve_des
        if type == "BOND":
            self.close_window(type)
        elif type == "CDS":
            print "XXXXXXXXXXXXXX caso CDS!!!!!!!!!"
            self.new_window = W_settoreRatingSelection(parent=self, type = type)



class W_setParameters(LabelFrame):
    '''
    classdocs
    '''  
    """
    def __init__ (self,  parent = None):
        if parent:
            self.master = parent.master
            parent.close_window()
            print 'aaaaaaaaaaaaaaaaaaaaaaa'

        #LabelFrame.__init__(self, self.master)
        LabelFrame.__init__(self, parent)
        self.parent=parent

        print 'bbbbbbbbbbbbbbbbbbbbbbbbbbb'

        self.initialize_user_interface()
    """


    
    def __init__(self, parent, hr_model):
        '''
        Constructor
        '''
        LabelFrame.__init__(self, parent)
        self.hr_model = hr_model
        self.parent=parent
        self.initialize_user_interface()
    


    def initialize_user_interface(self):
        """Draw a user interface allowing the user to type
        items and insert them into the treeview
        """
        self.parent.title("HR parameters")       
        self.parent.grid_rowconfigure(0,weight=1)
        self.parent.grid_columnconfigure(0,weight=1)
        #self.parent.config(background="lavender")


        # Define the different GUI widgets

        self.min_label = Label(self.parent, text = "min")
        self.min_label.grid(row = 0, column = 1, pady=10, padx=50, sticky = W)

        self.mean_label = Label(self.parent, text = "x0")
        self.mean_label.grid(row = 0, column = 2, pady=10, padx=50, sticky = W)
        
        self.max_label = Label(self.parent, text = "max")
        self.max_label.grid(row = 0, column = 3, pady=10, padx=50, sticky = W)
        
        #----------------------------------------------------------
        
        if (self.hr_model == 'SVE'):
        
            self.tau1_label = Label(self.parent, text = "Tau1:")
            self.tau1_label.grid(row = 1, column = 0, sticky = W)
    
            self.tau1_min = Entry(self.parent)
            self.tau1_min.grid(row = 1, column = 1)
            self.tau1_min.insert(0, 0.0001)
    
            self.tau1_x0 = Entry(self.parent)
            self.tau1_x0.grid(row = 1, column = 2)
            self.tau1_x0.insert(0, 1.0)
    
            self.tau1_max = Entry(self.parent)
            self.tau1_max.grid(row = 1, column = 3)
            self.tau1_max.insert(0, 10.0)
            
            #----------------------------------------------------------
    
            self.tau2_label = Label(self.parent, text = "Tau2:")
            self.tau2_label.grid(row = 2, column = 0, sticky = W)
    
            self.tau2_min = Entry(self.parent)
            self.tau2_min.grid(row = 2, column = 1)
            self.tau2_min.insert(0, 0.0001)
    
            self.tau2_x0 = Entry(self.parent)
            self.tau2_x0.grid(row = 2, column = 2)
            self.tau2_x0.insert(0, 10.0)
    
            self.tau2_max = Entry(self.parent)
            self.tau2_max.grid(row = 2, column = 3)
            self.tau2_max.insert(0, 50.0)
    
            #----------------------------------------------------------
    
            self.b0_label = Label(self.parent, text = "b0:")
            self.b0_label.grid(row = 3, column = 0, sticky = W)
    
            self.b0_min = Entry(self.parent)
            self.b0_min.grid(row = 3, column = 1)
            self.b0_min.insert(0, -10.0)
    
            self.b0_x0 = Entry(self.parent)
            self.b0_x0.grid(row = 3, column = 2)
            self.b0_x0.insert(0, 0.03)
    
            self.b0_max = Entry(self.parent)
            self.b0_max.grid(row = 3, column = 3)
            self.b0_max.insert(0, 10.03)
    
            #----------------------------------------------------------
    
            self.b1_label = Label(self.parent, text = "b1:")
            self.b1_label.grid(row = 4, column = 0, sticky = W)
    
            self.b1_min = Entry(self.parent)
            self.b1_min.grid(row = 4, column = 1)
            self.b1_min.insert(0, -10.0)
    
            self.b1_x0 = Entry(self.parent)
            self.b1_x0.grid(row = 4, column = 2)
            self.b1_x0.insert(0, 0.03)
    
            self.b1_max = Entry(self.parent)
            self.b1_max.grid(row = 4, column = 3)
            self.b1_max.insert(0, 10.0)
    
            #----------------------------------------------------------
    
            self.b2_label = Label(self.parent, text = "b2:")
            self.b2_label.grid(row = 5, column = 0, sticky = W)
    
            self.b2_min = Entry(self.parent)
            self.b2_min.grid(row = 5, column = 1)
            self.b2_min.insert(0, -10.0)
    
            self.b2_x0 = Entry(self.parent)
            self.b2_x0.grid(row = 5, column = 2)
            self.b2_x0.insert(0, 0.03)
    
            self.b2_max = Entry(self.parent)
            self.b2_max.grid(row = 5, column = 3)
            self.b2_max.insert(0, 10.5)
    
            #----------------------------------------------------------
    
            self.b3_label = Label(self.parent, text = "b3:")
            self.b3_label.grid(row = 6, column = 0, sticky = W)
    
            self.b3_min = Entry(self.parent)
            self.b3_min.grid(row = 6, column = 1)
            self.b3_min.insert(0, -10.0)
    
            self.b3_x0 = Entry(self.parent)
            self.b3_x0.grid(row = 6, column = 2)
            self.b3_x0.insert(0, 0.03)
    
            self.b3_max = Entry(self.parent)
            self.b3_max.grid(row = 6, column = 3)
            self.b3_max.insert(0, 10.0)

        if (self.hr_model == 'NS'):
        
            self.tau1_label = Label(self.parent, text = "Tau1:")
            self.tau1_label.grid(row = 1, column = 0, sticky = W)
    
            self.tau1_min = Entry(self.parent)
            self.tau1_min.grid(row = 1, column = 1)
            self.tau1_min.insert(0, 0.0001)
    
            self.tau1_x0 = Entry(self.parent)
            self.tau1_x0.grid(row = 1, column = 2)
            self.tau1_x0.insert(0, 1.0)
    
            self.tau1_max = Entry(self.parent)
            self.tau1_max.grid(row = 1, column = 3)
            self.tau1_max.insert(0, 10.0)
            
            #----------------------------------------------------------
    
            self.b0_label = Label(self.parent, text = "b0:")
            self.b0_label.grid(row = 3, column = 0, sticky = W)
    
            self.b0_min = Entry(self.parent)
            self.b0_min.grid(row = 3, column = 1)
            self.b0_min.insert(0, -10.0)
    
            self.b0_x0 = Entry(self.parent)
            self.b0_x0.grid(row = 3, column = 2)
            self.b0_x0.insert(0, 0.03)
    
            self.b0_max = Entry(self.parent)
            self.b0_max.grid(row = 3, column = 3)
            self.b0_max.insert(0, 10.03)
    
            #----------------------------------------------------------
    
            self.b1_label = Label(self.parent, text = "b1:")
            self.b1_label.grid(row = 4, column = 0, sticky = W)
    
            self.b1_min = Entry(self.parent)
            self.b1_min.grid(row = 4, column = 1)
            self.b1_min.insert(0, -10.0)
    
            self.b1_x0 = Entry(self.parent)
            self.b1_x0.grid(row = 4, column = 2)
            self.b1_x0.insert(0, 0.03)
    
            self.b1_max = Entry(self.parent)
            self.b1_max.grid(row = 4, column = 3)
            self.b1_max.insert(0, 10.0)
    
            #----------------------------------------------------------
    
            self.b2_label = Label(self.parent, text = "b2:")
            self.b2_label.grid(row = 5, column = 0, sticky = W)
    
            self.b2_min = Entry(self.parent)
            self.b2_min.grid(row = 5, column = 1)
            self.b2_min.insert(0, -10.0)
    
            self.b2_x0 = Entry(self.parent)
            self.b2_x0.grid(row = 5, column = 2)
            self.b2_x0.insert(0, 0.03)
    
            self.b2_max = Entry(self.parent)
            self.b2_max.grid(row = 5, column = 3)
            self.b2_max.insert(0, 10.5)
    
            #----------------------------------------------------------
    
        if (self.hr_model == 'CIR'):
        
            self.r0_label = Label(self.parent, text = "r0:")
            self.r0_label.grid(row = 1, column = 0, sticky = W)
    
            self.r0_min = Entry(self.parent)
            self.r0_min.grid(row = 1, column = 1)
            self.r0_min.insert(0, 0.0001)
    
            self.r0_x0 = Entry(self.parent)
            self.r0_x0.grid(row = 1, column = 2)
            self.r0_x0.insert(0, 1.0)
    
            self.r0_max = Entry(self.parent)
            self.r0_max.grid(row = 1, column = 3)
            self.r0_max.insert(0, 10.0)
            
            #----------------------------------------------------------
    
            self.kappa_label = Label(self.parent, text = "kappa:")
            self.kappa_label.grid(row = 3, column = 0, sticky = W)
    
            self.kappa_min = Entry(self.parent)
            self.kappa_min.grid(row = 3, column = 1)
            self.kappa_min.insert(0, -10.0)
    
            self.kappa_x0 = Entry(self.parent)
            self.kappa_x0.grid(row = 3, column = 2)
            self.kappa_x0.insert(0, 0.03)
    
            self.kappa_max = Entry(self.parent)
            self.kappa_max.grid(row = 3, column = 3)
            self.kappa_max.insert(0, 10.03)
    
            #----------------------------------------------------------
    
            self.theta_label = Label(self.parent, text = "theta:")
            self.theta_label.grid(row = 4, column = 0, sticky = W)
    
            self.theta_min = Entry(self.parent)
            self.theta_min.grid(row = 4, column = 1)
            self.theta_min.insert(0, -10.0)
    
            self.theta_x0 = Entry(self.parent)
            self.theta_x0.grid(row = 4, column = 2)
            self.theta_x0.insert(0, 0.03)
    
            self.theta_max = Entry(self.parent)
            self.theta_max.grid(row = 4, column = 3)
            self.theta_max.insert(0, 10.0)
    
            #----------------------------------------------------------
    
            self.sigma_label = Label(self.parent, text = "sigma:")
            self.sigma_label.grid(row = 5, column = 0, sticky = W)
    
            self.sigma_min = Entry(self.parent)
            self.sigma_min.grid(row = 5, column = 1)
            self.sigma_min.insert(0, -10.0)
    
            self.sigma_x0 = Entry(self.parent)
            self.sigma_x0.grid(row = 5, column = 2)
            self.sigma_x0.insert(0, 0.03)
    
            self.sigma_max = Entry(self.parent)
            self.sigma_max.grid(row = 5, column = 3)
            self.sigma_max.insert(0, 10.5)
    
            #----------------------------------------------------------




        """
        
        if (model_fit == 'SVE'): #->SVE
        
            #bound_min_sve = [0.0001,  0.0001, -10.00, -10.050, -10.00, -10.00]
    
            bound_min = [0.0001,  0.0001, -10.00, -10.00, -10.00, -10.00]
            bound_max = [10.0,     50.0,  10.03,   10.0,   10.5,   10.0]
            x0  = [1.0,       10.0,   0.03,   0.03,   0.03,   0.03]
            n_par = 6
    
    
        elif (model_fit == 'CIR'): #->CIR
    
            bound_min = [  -0.1,  0.1, 0.001, 0.001]
            bound_max = [ 10.00, 10.0, 10.00, 1.000]
            x0  = [0.0001,  1.0, 0.015,  0.01]
            n_par = 4
    
    
        elif (model_fit == 'NS'): #->NS
    
            bound_min = [0.0001,  -10.0, -10.0, -10.0]
            bound_max = [100,     +10.0, +10.0, +10.0]
            x0        = [5.0001,  0.03, 0.03, 0.03]
            n_par = 4
    


        self.dose_label = Label(self.parent, text = "a1:")
        self.dose_entry = Entry(self.parent)
        self.dose_label.grid(row = 0, column = 0, sticky = W)
        self.dose_entry.grid(row = 1, column = 0)
        """

        """
        self.modified_label = Label(self.parent, text = "Date Modified:")
        self.modified_entry = Entry(self.parent)
        self.modified_label.grid(row = 1, column = 0, sticky = W)
        self.modified_entry.grid(row = 1, column = 1)
        """

        self.submit_button = Button(self.parent, text = "Validate", command = self.insert_data)
        self.submit_button.grid(row = 7, column = 2, pady=10, padx=50, sticky = W)
        #self.exit_button = Button(self.parent, text = "Exit", command = self.parent.quit)
        #self.exit_button.grid(row = 0, column = 3)

        """
        # Set the treeview
        self.tree = ttk.Treeview( self.parent, columns=('Dose', 'Modification date'))
        self.tree.heading('#0', text='Item')
        self.tree.heading('#1', text='Dose')
        self.tree.heading('#2', text='Modification Date')
        self.tree.column('#1', stretch=YES)
        self.tree.column('#2', stretch=YES)
        self.tree.column('#0', stretch=YES)
        self.tree.grid(row=4, columnspan=4, sticky='nsew')
        self.treeview = self.tree
        # Initialize the counter
        self.i = 0
        """


    def insert_data(self):
        """
        Insertion method.
        """
        

        model_fit = self.hr_model
        
        if (model_fit == 'SVE'): #->SVE

            tau1_min = float(self.tau1_min.get())
            tau2_min = float(self.tau2_min.get())
            b0_min = float(self.b0_min.get())
            b1_min = float(self.b1_min.get())
            b2_min = float(self.b2_min.get())
            b3_min = float(self.b3_min.get())

            tau1_max = float(self.tau1_max.get())
            tau2_max = float(self.tau2_max.get())
            b0_max = float(self.b0_max.get())
            b1_max = float(self.b1_max.get())
            b2_max = float(self.b2_max.get())
            b3_max = float(self.b3_max.get())

            tau1_x0 = float(self.tau1_x0.get())
            tau2_x0 = float(self.tau2_x0.get())
            b0_x0 = float(self.b0_x0.get())
            b1_x0 = float(self.b1_x0.get())
            b2_x0 = float(self.b2_x0.get())
            b3_x0 = float(self.b3_x0.get())

            bound_min = [tau1_min,  tau2_min, b0_min, b1_min, b2_min, b3_min]
            bound_max = [tau1_max,  tau2_max, b0_max, b1_max, b2_max, b3_max]
            self.x0        = [tau1_x0,  tau2_x0, b0_x0, b1_x0, b2_x0, b3_x0]
            n_par = 6

            
            
        elif (model_fit == 'NS'): #->SVE


            tau1_min = self.tau1_min.get()
            b0_min = self.b0_min.get()
            b1_min = self.b1_min.get()
            b2_min = self.b2_min.get()

            tau1_max = self.tau1_max.get()
            b0_max = self.b0_max.get()
            b1_max = self.b1_max.get()
            b2_max = self.b2_max.get()

            tau1_x0 = self.tau1_x0.get()
            b0_x0 = self.b0_x0.get()
            b1_x0 = self.b1_x0.get()
            b2_x0 = self.b2_x0.get()

            bound_min = [tau1_min,  b0_min, b1_min, b2_min]
            bound_max = [tau1_max,  b0_max, b1_max, b2_max]
            self.x0        = [tau1_x0,  b0_x0, b1_x0, b2_x0]
            n_par = 4


        elif (model_fit == 'CIR'): #->SVE


            r0_min = self.r0_min.get()
            kappa_min = self.kappa_min.get()
            theta_min = self.theta_min.get()
            sigma_min = self.sigma_min.get()

            r0_max = self.r0_max.get()
            kappa_max = self.kappa_max.get()
            theta_max = self.theta_max.get()
            sigma_max = self.sigma_max.get()

            r0_x0 = self.r0_x0.get()
            kappa_x0 = self.kappa_x0.get()
            theta_x0 = self.theta_x0.get()
            sigma_x0 = self.sigma_x0.get()

            bound_min = [r0_min,  kappa_min, theta_min, sigma_min]
            bound_max = [r0_max,  kappa_max, theta_max, sigma_max]
            self.x0        = [r0_x0,  kappa_x0, theta_x0, sigma_x0]
            n_par = 4


        else: #->SVE

            print 'Modello non previsto!!!'
            #x0  = [1.0,       10.0,   0.03,   0.03,   0.03,   0.03]

   
        self.x_bnd = []
        for i in range(0, n_par):
    
            bndTmp = []
            
            b_min = bound_min[i]
            b_max = bound_max[i]
            
            bndTmp.append(b_min)
            bndTmp.append(b_max)
            
            self.x_bnd.append(bndTmp)

        
            """
            bound_min = [0.0001,  0.0001, -10.00, -10.00, -10.00, -10.00]
            bound_max = [10.0,     50.0,  10.03,   10.0,   10.5,   10.0]
            x0  = [1.0,       10.0,   0.03,   0.03,   0.03,   0.03]
            n_par = 6
            """


        
        
        
        
        
        self.destroy()
        self.master.destroy()
        
        
        """
        def set_prms_for_fit(model_fit):
        
        if (model_fit == 'SVE'): #->SVE
        
            bound_min = [0.0001,  0.0001, -10.00, -10.00, -10.00, -10.00]
            bound_max = [10.0,     50.0,  10.03,   10.0,   10.5,   10.0]
            x0  = [1.0,       10.0,   0.03,   0.03,   0.03,   0.03]
            n_par = 6
        
        
        elif (model_fit == 'CIR'): #->CIR
        
            bound_min = [  -0.1,  0.1, 0.001, 0.001]
            bound_max = [ 10.00, 10.0, 10.00, 1.000]
            x0  = [0.0001,  1.0, 0.015,  0.01]
            n_par = 4
        
        
        elif (model_fit == 'NS'): #->NS
        
            bound_min = [0.0001,  -10.0, -10.0, -10.0]
            bound_max = [100,     +10.0, +10.0, +10.0]
            x0        = [5.0001,  0.03, 0.03, 0.03]
            n_par = 4
        
        
        x_bnd = []
        for i in range(0, n_par):
        
            bndTmp = []
            
            b_min = bound_min[i]
            b_max = bound_max[i]
            
            bndTmp.append(b_min)
            bndTmp.append(b_max)
            
            x_bnd.append(bndTmp)
        
        
        return x0, x_bnd
        """
        
        
        #self.treeview.insert('', 'end', text="Item_"+str(self.i), values=(self.dose_entry.get()+" mg", self.modified_entry.get()))
        # Increment counter
        #self.i = self.i + 1         







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



#from pyxll import xl_func, xl_app, xl_menu, xl_macro, xlcAlert
from xls_utils import *
from Tkinter import *
import tkMessageBox
from db_qrys import getCurvesListFromDb, getDatesListFromDb
from sc_elab.core.SwpCurve import *
from win32com.client import constants as const



#def donothing():
#    tkMessageBox.showinfo("Nothing To do", "bye bye")


class W_bondFitting(LabelFrame):
    def __init__(self, master = None, curveL = [], type = "SWP"):
        c_date = None
        LabelFrame.__init__(self, master)
        self.master = master
        #self.geometry("800x600")
        #self.master.geometry("400x500")
        

        if (type != "BOND"):
            self.config(text="Available curves for bootstrapping:")
        else:
            self.config(text="Available portfolios for bond fitting:")
        
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
    
    
    
    def close_window(self):
        self.destroy()


    def sel(self):


        hr_model_map = {}
        hr_model_map[0] = 'SVE'
        hr_model_map[1] = 'NS'
        hr_model_map[2] = 'CIR'


        hr_model_tmp   = (str(self.variable1.get()).strip(""))[1]
        
        
        self.hr_model = hr_model_map[int(hr_model_tmp)]
              
        root = Tk()
        self.window_prms =  W_setParameters(root, self.hr_model)
        
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

        if type != "BOND":
                self.config(text="Set Bootstrap Options:", width=400, height=200)                
        else:
                self.config(text="Set Bond Fitting Options:", width=400, height=200)                
            
        #--- Boot swaps rates

        if type == "SWP":

            T1 = Label(self,height=1, width=30, text = "Bootstrap method (filling):").grid(row=0, sticky = "e")
            self.variable1 = StringVar(self)
            self.variable1.set("(0) Costant Fwd Rate")  # default value
            w1 = OptionMenu(self, self.variable1, "(0) Costant Fwd py Rate", "(1) Linear py Rate")
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

        elif type == "CDS":
            
            T1 = Label(self,height=1, width=30, text = "Hazard Rate interpolation model:").grid(row=0, sticky = "e")
            self.variable1 = StringVar(self)
            self.variable1.set("(0) Linear")  # default value
            w1 = OptionMenu(self, self.variable1,  "(0) LIN", "(1) AVD", "(2) SVE", "(3) CIR", "(4) NS")
            w1.grid(row=0, column=1)
            w1.config(width=30)
            
            
            T6 = Label(self, height=1, width=30, text="Bench. interpolation model:").grid(row=1, sticky="e")
            self.variable6 = StringVar(self)
            self.variable6.set("(0) Linear")  # default value
            w6 = OptionMenu(self, self.variable6, "(0) LIN", "(1) AVD", "(2) SVE", "(3) CIR", "(4) NS")
            w6.grid(row=1, column=1)
            w6.config(width=30)

            
            # ---
            T7 = Label(self, height=1, width=30, text="Bootstrap Method for HR:").grid(row=2, sticky="e")
            self.variable7 = StringVar(self)
            self.variable7.set("(0) Constant Hazard Rate")  # default value
            w7 = OptionMenu(self, self.variable7, "(0) Constant Hazard Rate", "(0) Linear Credit spread")
            w7.grid(row=2, column=1)
            w7.config(width=30)
            # ---

        elif type == "BOND":
            
            T1 = Label(self,height=1, width=30, text = "Hazard Rate interpolation model:").grid(row=0, sticky = "e")
            self.variable1 = StringVar(self)
            self.variable1.set("(0) Svensson")  # default value
            w1 = OptionMenu(self, self.variable1,  "(0) SVE", "(1) NS", "(2) CIR")
            w1.grid(row=0, column=1)
            w1.config(width=30)
            
            #fileMenu = Menu(self)
            
            #w1.add_cascade(label="File", menu=fileMenu)

            #T7 = Label(self,height=1, width=30, text = "Model evaluation:").grid(row=0, sticky = "e")
            #self.variable2 = StringVar(self)
            #self.variable2.set("(0) Linear")  # default value
            #w7 = OptionMenu(self, self.variable7,  "(0) RMV", "(1) RFV")
            #w7.grid(row=0, column=1)
            #w7.config(width=30)
            
            
            T6 = Label(self, height=1, width=30, text="Risk free interpolation model:").grid(row=1, sticky="e")
            self.variable6 = StringVar(self)
            self.variable6.set("(0) Linear")  # default value
            w6 = OptionMenu(self, self.variable6, "(0) LIN", "(1) AVD", "(2) SVE", "(3) CIR", "(4) NS")
            w6.grid(row=1, column=1)
            w6.config(width=30)

            
            # ---
            
            T7 = Label(self, height=1, width=30, text="Bond evaluation model:").grid(row=2, sticky="e")
            self.variable7 = StringVar(self)
            self.variable7.set("(0) Recovery market value")  # default value
            w7 = OptionMenu(self, self.variable7, "(0) RMV", "(1) RFV")
            w7.grid(row=2, column=1)
            w7.config(width=30)
            
            
            


        self.pack(fill="both", expand="yes")

        #questa istruzione va messa dopo il pack altrimenti non vienne intercettato il valore dal .get() successivo
        if type == "SWP": self.variable5.set("C://")
        
        #filemenu = Menu(menubar, tearoff=0)

        #B0 = Menu(self, text="Params",  tearoff=0)

        B1 = Button(self, text="Submit",  width=20, command=self.sel).grid(row=5, column=1, sticky ='e')

        self.pack(fill="both", expand="yes")
        
