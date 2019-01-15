from Tkinter import *
from sc_elab.excel_hook.connection import Connection


def getDatesListFromDbSwaptions():

    con = Connection()

    cursor = con.db_data()

    qry_to_execute = '''SELECT distinct Data FROM DProTS_master 
                        WHERE TipoDato='VSwaption'
                        ORDER by Data desc
                     '''

    cursor.execute(qry_to_execute)
    result = cursor.fetchall()

    date_list   = []
    for record in result:
        date_list.append(record[0])
    return date_list


class W_SwaptionsDate(LabelFrame):
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
        dl = getDatesListFromDbSwaptions()
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


class W_SwaptionsOptionsPrint(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.master = master

        self.print_type = StringVar()
        self.print_type.set('matrix')

        Label(self, text="""Choose format print :""", justify=LEFT, padx=100).pack()

        self.rb_calib1 = Radiobutton(self, text='matrix', justify='left', variable=self.print_type,
                                     value='matrix').pack(anchor=W)
        self.rb_calib2 = Radiobutton(self, text='table', justify='left', variable=self.print_type,
                                     value='table').pack(anchor=W)

        # create button
        self.btn2 = Button(self, text="Cancel", command=lambda: self.close_window())
        self.btn2.pack(side=BOTTOM, fill='x')
        # create button
        self.btn1 = Button(self, text="Select", command=lambda: self.load_curve())
        self.btn1.pack(side=BOTTOM, fill='x')
        self.pack()

    def load_curve(self):
        print "option print", self.print_type.get()
        self.master.destroy()

    def close_window(self):
        self.master.destroy()

