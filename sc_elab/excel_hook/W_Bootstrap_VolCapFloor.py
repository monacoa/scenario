from Tkinter import *
import ttk
from sc_elab.excel_hook.connection import Connection


def getDatesListFromDbVolCapFloor():

    con = Connection()

    cursor = con.db_data()

    # Controlla la query
    qry_to_execute = '''SELECT distinct Data FROM DProTS_master 
                        WHERE TipoDato='VCapFloor'
                        ORDER by Data desc
                     '''

    cursor.execute(qry_to_execute)
    result = cursor.fetchall()

    date_list   = []
    for record in result:
        date_list.append(record[0])
    return date_list


# Finestra di selezione dei dati
class W_VolCapFloorDate(LabelFrame):
    def  __init__(self, master = None,parent=None):
        c_date = None
        if master:
            self.master = master
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
        dl = getDatesListFromDbVolCapFloor()
        for d in dl:
            #date = d.date()
            self.mylist.insert(END, d)

        self.mylist.pack(side=LEFT, fill='y')
        self.bar.pack(side=LEFT, fill='y')

        self.mylist.config(yscrollcommand=self.bar.set)
        self.bar.config(command=self.mylist.yview)
        # -----------
        # cretae button
        self.btn2 = Button(self, text="Cancel",  command=self.close_window)
        self.btn2.pack(side=BOTTOM, fill = 'x')
        # cretae button
        self.btn1 = Button(self, text="Select", command= lambda:self.selected_date())
        self.btn1.pack(side=BOTTOM, fill='x')

    def close_window(self):
        self.master.destroy()

    def selected_date(self):
        #recupero la data selezionata
        c_date = self.mylist.get(ACTIVE)
        #self.date = c_date
        self.date = c_date#.datetime()
        self.master.destroy()


# Finestra di selezione della Currency
class W_CurrencySelection(LabelFrame):
    def __init__(self, parent = None, currencies = None, contributors=None):
        # viene chiusa la finestra aperta
        if parent:
            self.master = parent.master
            parent.destroy()
        LabelFrame.__init__(self, self.master)
        self.config(text="Currency selection:")

        # Scrivo i menu di scelta
        Label(self, height=1, width=30, text="Currency").grid(row=0, sticky="e")
        self.curr = StringVar()
        self.curr.set(currencies[0])  # default value
        w1 = ttk.Combobox(self)
        w1.config(textvariable=self.curr, state="readonly", values=currencies)
        # w2.grid(row=6, column=2, rowspan=1, columnspan=3, pady=2, sticky=W + E + N + S)
        w1.grid(row=0, column=1, sticky=W + E + N + S)
        w1.config(width=30)

        Label(self, height=1, width=30, text="Contributor").grid(row=1, sticky="e")
        self.cont = StringVar()
        self.cont.set(contributors[0])  # default value
        w2 = ttk.Combobox(self)
        w2.config(textvariable=self.cont, state="readonly", values=contributors)
        #w2.grid(row=6, column=2, rowspan=1, columnspan=3, pady=2, sticky=W + E + N + S)
        w2.grid(row=1, column=1, sticky=W + E + N + S)
        w2.config(width=30)

        # Scrivo il pulsante
        B1 = Button(self, text="Submit", width=20, command=self.sel).grid(row=2, column=1, sticky='e')
        self.pack(fill="both", expand="yes")
        # create scrollbar


    def close_window(self):
        self.destroy()
        self.master.destroy()

    def sel(self):
        self.close_window()



# Finestra di selezione curva dei fattori di sconto e dati sulle volatilita'
def Bootstrap_BVol_menu(volsdata = None, discount_curves = None):

    list_choices=[]

    def close_window():
        list_choices.append(1)
        list_choices.append(curve_disc.get())
        list_choices.append(volatilities.get())
        c.destroy()
        c.master.destroy()

    def close_without_selection():
        list_choices.append(0)
        c.destroy()
        c.master.destroy()

    root=Tk()

    root.title('Curve and volatilities data selection')
    root.geometry('400x140')
    # Create the grid and the outer content frame
    c = ttk.Frame(root, padding=(5,5,12,0))
    c.grid(column=0, row=0, sticky=(N,W,E,S))
    root.grid_columnconfigure(0, weight=1)
    root.grid_rowconfigure(0, weight=1)

    LabelCurve = ttk.Label(c, text='Curve')
    LabelCurve.grid(row=1, column=0, rowspan=1, columnspan=1, pady=2, sticky=W + E + N + S)

    curve_disc = ttk.Combobox(c)
    CurveName = StringVar()
    CurveName.set(discount_curves[0])
    curve_disc.config(textvariable=CurveName, state="readonly", values=discount_curves)
    curve_disc.grid(row=1, column=1, rowspan=1, columnspan=1, pady=2, sticky=W + E + N + S)

    LabelOptions = ttk.Label(c, text='Volatilities')
    LabelOptions.grid(row=2, column=0, rowspan=1, columnspan=1, pady=2, sticky=W + E + N + S)

    volatilities = ttk.Combobox(c)
    VolName = StringVar()
    VolName.set(volsdata[0])
    volatilities.config(textvariable=VolName, state="readonly", values=volsdata)
    volatilities.grid(row=2, column=1, rowspan=1, columnspan=1, pady=2, sticky=W + E + N + S)

    # Bottoni
    SubmitButton = ttk.Button(c, text='Submit', command = close_window)
    SubmitButton.grid(row=3, column=0, rowspan=1, columnspan=1, pady=2, sticky=W + E + N + S)

    CancelButton = ttk.Button(c, text='Cancel', command = close_without_selection)
    CancelButton.grid(row=3, column=1, rowspan=1, columnspan=1, pady=2, sticky=W + E + N + S)

    root.mainloop()

    return list_choices