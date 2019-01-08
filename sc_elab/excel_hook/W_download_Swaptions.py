import numpy as np
from Tkinter import *
from sc_elab.excel_hook.xls_Calibration import findCalibrationPos,writeResultPandas
from sc_elab.excel_hook.xls_utils import drawBox, formatTestataCurva
from win32com.client import constants as const
import datetime
from DEF_intef import nameSheetScaricoSwaption

from sc_elab.excel_hook.connection import Connection

##### lettura dati da interfaccia

def getDatesListFromDbSwaptions():

    con = Connection()

    cursor = con.db_data()

    qry_to_execute = '''SELECT distinct Data FROM dprots_master 
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


##### scrittura dati su foglio excel

def intestazioneSwaptions( xla, rng,  attributi, nCols = 2, title= 'Matrix Swaption'):

    nRows           = len(attributi.keys())
    topLeftRow      = rng.Row
    topLeftCol      = rng.Column
    drawBox            (xla, 3,topLeftRow, topLeftCol,topLeftRow + nRows, topLeftCol + nCols - 1, 0)
    formatTestataCurva (xla, topLeftRow, topLeftCol, nCols, title)

    kk = attributi.keys()
    kk.sort()
    i = 0
    for k in kk:
        xla.Cells(topLeftRow + 1+ i, topLeftCol).Value   = k[3:]
        xla.Cells(topLeftRow + 1+ i, topLeftCol+1).Value = attributi[k]
        if isinstance(attributi[k],datetime.datetime):
            xla.Cells(topLeftRow + 1 + i, topLeftCol + 1).NumberFormat = "gg/MM/aaaa"
        xla.Cells(topLeftRow + 1 + i, topLeftCol + 1).HorizontalAlignment = const.xlCenter
        i+=1

    rangeStart = xla.Range(xla.Cells(topLeftRow + nRows + 1, topLeftCol),xla.Cells(topLeftRow + nRows + 1, topLeftCol))
    return rangeStart


def writeSwaptionsResOnXls(data, xla, ref_date, option_print):

    r = findCalibrationPos(xla, nameSheetScaricoSwaption)

    #---
    #mi posiziono nella prima cella utile per scrivere i risultati del fitting
    #---

    if option_print == 'matrix':
        Attributi = \
            {     "1. Date ref"    : ref_date
                , "2. Tipo Dato"   : 'VSwaption'
                , "3. Valore"      : 'MID'
                , "4. Contributor" : ''
                , "5. Currency"    : ''
                , "6. Rows"        : 'Tenor'
                , "7. Columns"     : 'Maturity'
                  }
    else:
        Attributi = \
            {     "1. Date ref"    : ref_date
                , "2. Tipo Dato"   : 'VSwaption'
                , "3. Valore"      : 'MID'
                , "4. Contributor" : ''
                , "5. Currency"    : ''
                  }

    r = intestazioneSwaptions(xla = xla, rng = r, attributi = Attributi)
    r = writeResultPandas(xla = xla , rng = r, df = data)


def write_Swaptions(xla, res, ref_date, option_print = 'matrix'):
    res = res.astype('float')
    if option_print == 'matrix':

        res3 = res.pivot_table(index='Tenor', columns='MaturityInt', values='ValoreMid',fill_value = -1)

        #cambio il format della matrice
        res3.columns = res3.columns.astype(np.float)
        res3.index = res3.index.astype(np.float)
        res3.reset_index(level='Tenor', inplace=True)

        writeSwaptionsResOnXls(res3, xla, ref_date, option_print)

    else:
        res['Usage'] = 'Y'
        writeSwaptionsResOnXls(res,  xla, ref_date, option_print)
