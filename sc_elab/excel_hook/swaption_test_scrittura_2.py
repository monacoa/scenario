#from pyxll import xlcAlert
import numpy as np
import pandas as pd
from DEF_intef import nameSheetScaricoSwaption
import mysql.connector
from Tkinter import *
import datetime
import tkMessageBox
from db_qrys import getCurvesListFromDb, getDatesListFromDb, getProvidersFromDb


def getDatesListFromDb(db_dict):

    db = mysql.connector.connect(host="localhost", user="root", password="DatabaseRepl1ca",
                             database="db_mercato_matteo")
    cursor = db.cursor()

    qry_to_execute = '''SELECT distinct Data FROM db_mercato_matteo.dprots_master 
                       WHERE TipoDato='VSwaption' 
                        ORDER by Data desc
                       '''


    cursor.execute(qry_to_execute)
    result = cursor.fetchall()

    date_list   = []
    for record in result:
        date_list.append(record[0])
    return date_list


class W_curveDate(LabelFrame):
    def  __init__(self, master = None,parent=None, db_dict = {}):
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
        dl = getDatesListFromDb(db_dict)
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
        print 'qui io sono'
        self.master.destroy()


def write_matrix(Sheet,xla):
#def write_matrix():

    db_credential = {}
        #db_credential['host'] = "10.103.65.195"
        #db_credential['user'] = "xxx3"
        #db_credential['pwd'] = "yyy2"


    db_credential['host'] = "localhost"
    db_credential['user'] = "root"
    db_credential['pwd'] = "DatabaseRepl1ca"
    db_name = "db_mercato_matteo"

    db = mysql.connector.connect(host=db_credential['host'],
                                 user=db_credential['user'],
                                 passwd=db_credential['pwd'],
                                 database=db_name)

    cursor = db.cursor()

    root = Tk()
    app = W_curveDate(root, db_dict= db_credential)
    root.mainloop()

    ref_date = datetime.date(day=int(app.date[-2:]), month=int(app.date[5:7]), year=int(app.date[:4]))
    print 'qui io sono parte 2 la vendetta'


    #Sheet = nameSheetScaricoSwaptionname
    #file_new_data = 'input/files_caricamento_datastream/test_scarico_swaption.xlsx'

    qry_to_execute= '''
                 SELECT DProCFS.Tenor, DProCFS.MaturityInt, DProTS_Master.ValoreMid
                 FROM DProCFS, DProTS_Master
                 WHERE DProTS_Master.BloombergTicker = DProCFS.BloombergTicker
                 AND(DProCFS.TipoDato = 'VSwaption') 
                 and DPROTS_Master.Data= '%s' ''' %(ref_date)

    # AND(DProTS_Master.Data ='2008-07-10')



    print 'sono arrivato fino a qui'

    res = pd.read_sql(qry_to_execute, db)
    res3 = res.pivot_table(index='Tenor', columns='MaturityInt', values='ValoreMid',fill_value=-1)


    lunghezza = len(res3.index)
    larghezza = len(res3.columns)

    #a=np.random.standard_normal((20,14))
    #df = pd.DataFrame(a)
    #print df
    lunghezza = res3.shape[0]
    larghezza = res3.shape[1]
    #appoggio= df.index.values
    for i in xrange(0,res3.shape[0]):
        # i indica le righe
        appoggio7=res3.index.values[i].__float__()
        xla.Cells(i+2,1).Value=res3.index.values[i].__float__()
        for j in xrange(0, res3.shape[1]):
                #j indica le colonne
                #appoggio5=res3.columns.values[j].__float__()
                #appoggio6=res3.iloc[i,j]
            xla.Cells(1, j+2).Value = res3.columns.values[j].__float__()
            xla.Cells(i+2, j+2).Value = res3.iloc[i,j]
                #xla.Cells(j+2, i).Value = res3.columns.values[j].__float__()
                #xla.Cells(i+2, j+2).Value = res3.iloc[i,j]

#if __name__ == "__main__":
#       write_matrix()

