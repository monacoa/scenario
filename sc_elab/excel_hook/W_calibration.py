from Tkinter import *
import tkMessageBox
import numpy as np
import pandas as pd
import datetime

def FQ(label):
    print ('------------- FIN QUI TUTTO OK  %s ----------' %(label))
    sys.exit()

def model_parameters(value):
    if value == 'CIR':
        dict = {'r0'    :{'sv':'0.03','min':'0.000001', 'max':'0.5' ,'fix':0},
                'k'     :{'sv':'0.5' ,'min':'0.001'   , 'max':'10.0','fix':0},
                'sigma' :{'sv':'0.05','min':'0.005'   , 'max':'1.0' ,'fix':0},
                'theta' :{'sv':'0.03','min':'0.001'   , 'max':'0.3' ,'fix':0}
                }

        # mi servono ordinati
        names = ['r0', 'k', 'theta', 'sigma']
        attribute = ['sv', 'min', 'max', 'fix']

    elif value == 'VSCK':
        dict = {'r0'    :{'sv':'0.03', 'min': '-0.5'  , 'max':'0.5'  ,'fix':0},
                'k'     :{'sv':'0.03', 'min': '0.0001', 'max':'10.0' ,'fix':0},
                'theta' :{'sv':'0.03', 'min': '0.0001', 'max': '10.0','fix':0},
                'sigma' :{'sv':'0.1' , 'min': '0.0001', 'max': '10.0','fix':0}
                }

        # mi servono ordinati
        names = ['r0', 'k', 'theta', 'sigma']
        attribute = ['sv', 'min', 'max', 'fix']

    elif value == 'Jarrow Yildirim':
        dict ={'aN'   :{'sv':'0.07' , 'min': '-0.8'  , 'max': '5.0' ,'fix':0},
              'aR'    :{'sv':'0.02' , 'min': '-0.8'  , 'max': '5.0' ,'fix':0},
              'sigmaR':{'sv':'0.009', 'min': '0.0001', 'max': '0.01','fix':0},
              'sigmaN':{'sv':'0.009', 'min': '0.0001', 'max': '0.01','fix':0},
              'sigmaI':{'sv':'0.009', 'min': '0.0001', 'max': '0.01','fix':0},
              'rhoNR' :{'sv':'0.2'  , 'min': '-0.5'  , 'max': '0.5' ,'fix':0},
              'rhoNI' :{'sv':'0.2'  , 'min': '-0.5'  , 'max': '0.5' ,'fix':0},
              'rhoRI' :{'sv':'0.2'  , 'min': '-0.5'  , 'max': '0.5' ,'fix':0}
              }
        names = ['aN', 'aR', 'sigmaR', 'sigmaN', 'sigmaI', 'rhoNR', 'rhoNI', 'rhoRI']
        attribute = ['sv', 'min', 'max', 'fix']

    elif value == 'G2++':
        dict ={'a'      :{'sv':'0.1'  , 'min': '0.0001', 'max': '5.0','fix':0},
              'sigma'   :{'sv':'0.01'  , 'min': '0.0001', 'max': '0.7','fix':0},
              'b'       :{'sv': '0.1', 'min': '0.0001', 'max': '5.0','fix':0},
              'eta'     :{'sv': '0.01', 'min': '0.0001', 'max': '0.7','fix':0},
              'rho'     :{'sv': '0.1' , 'min': '-0.7'  , 'max': '0.7','fix':0}
              }
        names = ['a', 'sigma', 'b', 'eta', 'rho']
        attribute = ['sv', 'min', 'max', 'fix']

    else:
        dict ={}
        names = []
        attribute = []

    return dict,names,attribute


class W_calib_models(Frame):
    def close_window(self):
        self.res = 0
        self.master.destroy()


    def __init__(self, master = None,nameWorkbook = None ,nameWorksheet = None):
        Frame.__init__(self, master)
        self.master = master
        self.master.title("Calibration models")

        self.nameWorkbook = nameWorkbook
        self.nameWorksheet = nameWorksheet

        self.model = StringVar()
        self.model.set(0)

        Label(self,text="""Choose your calibrator :""",justify=LEFT,padx=20).pack()

        self.calib_avaible = ['CIR',
                         'VSCK',
                         'Jarrow Yildirim',
                         'G2++']

        for name_calib in self.calib_avaible:
            Radiobutton(self,
                        text=name_calib,
                        padx=20,
                        variable=self.model,
                        value=name_calib).pack(anchor=W)

        self.res = 1
        # create button
        self.btn2 = Button(self, text="Cancel",  command=lambda:self.close_window())
        self.btn2.pack(side=BOTTOM, fill = 'x')
        # create button
        self.btn1 = Button(self, text="Select", command= lambda:self.selected_calib(self.model.get()))
        self.btn1.pack(side=BOTTOM, fill='x')
        self.pack()

    def selected_calib(self, value):
        #recupero la calibrazione selezionata
        print 'Il modello selezionato per la calibrazione:', value

        if value in self.calib_avaible:
            self.newWindow = W_calib_menu(parent = self , model = value,nameWorkbook=self.nameWorkbook,nameWorksheet=self.nameWorksheet)
        else:
            tkMessageBox.showinfo("bye bye", "Nothing Selected!")
            self.close_window()


import ttk
class W_calib_menu(LabelFrame):

    def close_window(self):
        self.res = 0
        self.destroy()
        self.master.destroy()

    def actionTAB(self):
        if self.set_mkt_ts.get() == 'MKT':
            self.nb.tab(0,state = "normal")
            self.nb.tab(1,state = "disabled")
            self.nb.select(0)
        else:
            self.nb.tab(0,state = "disabled")
            self.nb.tab(1,state = "normal")
            self.nb.select(1)

    def __init__(self, parent = None , model = "",nameWorkbook = None ,nameWorksheet = None):
        #elimino la finestra precedente
        if parent:
            self.master = parent.master
            parent.destroy()

        font9 = "-family {Segoe UI} -size 10 -weight bold -slant roman"  \
            " -underline 0 -overstrike 0"

        #########################################################################
        # Inizializzo l'oggetto
        #########################################################################

        LabelFrame.__init__(self, self.master)
        self.grid()
        self.master.title("Configure calibration")


        self.nameWorkbook = nameWorkbook
        self.nameWorksheet = nameWorksheet

        # da alimentare con il nome del foglio da leggere
        # e dallo sheet
        objectOnSheetDictionary = readSheetObject(self.nameWorkbook,self.nameWorksheet)
        if objectOnSheetDictionary[1].empty:
            objectOnSheet = pd.DataFrame()
            tmpCurve = tmpOptions = tmpTS = []

        else:
            objectOnSheet = readFeaturesObject(objectOnSheetDictionary)
            tmpCurve = objectOnSheet.loc[objectOnSheet.TypeObject == 'Curve', 'Name'].tolist()
            tmpOptions = objectOnSheet.loc[objectOnSheet.TypeObject == 'Option', 'Name'].tolist()
            tmpTS = objectOnSheet.loc[objectOnSheet.TypeObject == 'TS', 'Name'].tolist()

        #########################################################################
        # Titolo della form
        #########################################################################

        label_tit = Label(self,font=font9,text='Model: ' + model,
                               anchor='nw',justify='left').grid(row = 0, column = 0, rowspan = 1, columnspan = 2,pady=5, sticky = W+E+N+S)

        #########################################################################
        # Area di scelta tra MKT e TS
        #########################################################################

        self.set_mkt_ts = StringVar()
        self.rb_type1 = Radiobutton(self,text='Market',justify='left',variable=self.set_mkt_ts,value='MKT',command=self.actionTAB)
        self.rb_type1.grid(row = 1, column = 1, rowspan = 1, columnspan = 1, pady =5, sticky = W+E+N+S)

        self.rb_type2 = Radiobutton(self,text='Time Series',justify='left',variable=self.set_mkt_ts,value='TS',command=self.actionTAB)
        self.rb_type2.grid(row = 1, column = 2, rowspan = 1, columnspan = 1, pady = 5, sticky = W+E+N+S)
        self.set_mkt_ts.set('MKT')

        #########################################################################
        # Inizilizzo l'area dei TAB: MKT, TS, PARAMS
        #########################################################################

        self.nb = ttk.Notebook(self)
        self.nb.grid(row = 3, column = 1, rowspan = 7, columnspan = 5, sticky = W+E+N+S)

        self.nb_t0 = Frame(self.nb)
        self.nb.add(self.nb_t0,text="MKT",compound="left")

        self.nb_t1 = Frame(self.nb)
        self.nb.add(self.nb_t1, text="TS",compound="left",state = "disabled")

        self.nb_t2 = Frame(self.nb)
        self.nb.add(self.nb_t2, text="PARAMS",compound="left")

        #########################################################################
        # Area MKT
        #########################################################################

        #### Area Calibration Type

        Label2 = Label(self.nb_t0,text='Calibration Type')
        Label2.grid(row = 0, column = 1, rowspan = 1, columnspan = 1, pady =5, sticky = W+N+S)

        self.mkt_calibration_type = StringVar()

        self.rb_calib1 = Radiobutton(self.nb_t0,text='Options',justify='left',variable=self.mkt_calibration_type,value='OPT')
        self.rb_calib1.grid(row = 1, column = 1, rowspan = 1, columnspan = 1, pady =2, sticky = W+N+S)

        self.rb_calib2 = Radiobutton(self.nb_t0,text='Curve',justify='left',variable=self.mkt_calibration_type,value='CURVE')
        self.rb_calib2.grid(row = 2, column = 1, rowspan = 1, columnspan = 1, pady =2, sticky = W+N+S)

        self.rb_calib3 = Radiobutton(self.nb_t0,text='Curve + Options',justify='left',variable=self.mkt_calibration_type,value='CURVE_OPT')
        self.rb_calib3.grid(row = 3, column = 1, rowspan = 1, columnspan = 1, pady =2, sticky = W+N+S)

        if model in ['CIR','VSCK']:
            self.rb_calib1.config(state = "disabled")
            # self.rb_calib3.config(state = "disabled")
            self.mkt_calibration_type.set('CURVE')
        elif model in ['G2++']:
            self.rb_calib1.config(state = 'disabled')
            self.rb_calib2.config(state = 'disabled')
            self.mkt_calibration_type.set('CURVE_OPT')
        else:
            self.mkt_calibration_type.set('OPT')

        #### Separatore Verticale

        self.TSeparator1 = ttk.Separator(self.nb_t0,orient="vertical")
        self.TSeparator1.grid(row = 1, column = 2, rowspan = 4, columnspan = 1, padx = 15, sticky = W+E+N+S)

        #### Area Loss Function

        Label1 = Label(self.nb_t0,text='Loss Function')
        Label1.grid(row = 0, column = 3, rowspan = 1, columnspan = 1, pady =5, sticky = W+N+S)

        self.loss_function_type = IntVar()

        self.rb_loss1 = Radiobutton(self.nb_t0,text= "Euclidean absolute",justify='left',variable=self.loss_function_type,value=1)
        self.rb_loss1.grid(row = 1, column = 3, rowspan = 1, columnspan = 1, pady =1, sticky = W+N+S)

        self.rb_loss2 = Radiobutton(self.nb_t0,text= "Euclidean relative",justify='left',variable=self.loss_function_type,value=2)
        self.rb_loss2.grid(row = 2, column = 3, rowspan = 1, columnspan = 1, pady =1, sticky = W+N+S)

        self.rb_loss3 = Radiobutton(self.nb_t0,text= "Manhattan absolute",justify='left',variable=self.loss_function_type,value=3)
        self.rb_loss3.grid(row = 3, column = 3, rowspan = 1, columnspan = 1, pady =1, sticky = W+N+S)

        self.rb_loss4 = Radiobutton(self.nb_t0,text= "Manhattan relative",justify='left',variable=self.loss_function_type,value=4)
        self.rb_loss4.grid(row = 4, column = 3, rowspan = 1, columnspan = 1, pady =1, sticky = W+N+S)

        self.loss_function_type.set(1)

        #### Alimentazione della curva di input da foglio Excel
        Label(self.nb_t0, text='').grid(row = 5, column = 0, rowspan = 1, columnspan = 1, sticky = W+E+N+S)

        Label3 = Label(self.nb_t0, text='Curve')
        Label3.grid(row = 6, column = 1, rowspan = 1, columnspan = 1, pady =2, sticky = W+E+N+S)

        self.cb_curve = ttk.Combobox(self.nb_t0)
        self.NameCurve = StringVar()
        self.cb_curve.config(textvariable = self.NameCurve, state = "readonly", values = tmpCurve)
        self.cb_curve.grid(row = 6, column = 2, rowspan = 1, columnspan = 3, pady =2, sticky = W+E+N+S)

        #### Alimentazione delle opzioni di input da foglio Excel

        Label4 = Label(self.nb_t0,text='Options')
        Label4.grid(row = 7, column = 1, rowspan = 1, columnspan = 1, pady =2, sticky = W+E+N+S)

        self.cb_options = ttk.Combobox(self.nb_t0)
        self.NameOption = StringVar()
        self.cb_options.config(textvariable = self.NameOption, state = "readonly", values = tmpOptions)
        self.cb_options.grid(row = 7, column = 2, rowspan = 1, columnspan = 3, pady =2, sticky = W+E+N+S)

        #########################################################################
        # Area TS
        #########################################################################

        ####
        Label(self.nb_t1, text='').grid(row=0, column=0, rowspan=1, columnspan=1, sticky=W + E + N + S)

        Label13 = Label(self.nb_t1, text='Time Series')
        Label13.grid(row=1, column=1, rowspan=1, columnspan=1, pady=2, sticky=W + N + S)

        self.cb_ts = ttk.Combobox(self.nb_t1)
        self.NameTS = StringVar()
        self.cb_ts.config(textvariable=self.NameTS, state="readonly", values=tmpTS)
        self.cb_ts.grid(row=1, column=2, rowspan=1, columnspan=1, pady=2, sticky=W + N + S)
        self.cb_ts.bind("<<ComboboxSelected>>",
                        lambda _: self.list_combobox_ts(dictObject=objectOnSheetDictionary, tableObject=objectOnSheet))

        Label(self.nb_t1, text='').grid(row=2, column=0, rowspan=1, columnspan=1, sticky=W + E + N + S)

        ####
        Label14 = Label(self.nb_t1, text='Date min')
        Label14.grid(row=3, column=1, rowspan=1, columnspan=1, pady=2, sticky=W + N + S)

        self.te_dateMin = ttk.Combobox(self.nb_t1)
        self.te_dateMin.grid(row=3, column=2, rowspan=1, columnspan=1, pady=2, sticky=W + N + S)
        self.te_dateMin.configure(width=15)

        ####
        Label15 = Label(self.nb_t1, text='Date max')
        Label15.grid(row=4, column=1, rowspan=1, columnspan=1, pady=2, sticky=W + N + S)

        self.te_dateMax = ttk.Combobox(self.nb_t1)
        self.te_dateMax.grid(row=4, column=2, rowspan=1, columnspan=1, pady=2, sticky=W + N + S)
        self.te_dateMax.configure(width=15)

        #########################################################################
        # Area PARAMS
        #########################################################################

        # carico i parametri per il modello selezionato
        self.param_dict, self.params_names, self.params_attribute = model_parameters(model)

        # creo la tabella in base al modello
        self.create_range_parameter(self.nb_t2, self.params_names, self.params_attribute)

        #########################################################################
        # Bottoni finali
        #########################################################################
        self.res = 1

        self.bSubmit = Button(self, text='Submit', command=lambda: self.go_to_calib(dictObject=objectOnSheetDictionary,
                                                                                    tableObject=objectOnSheet))
        self.bSubmit.grid(row=10, column=4, rowspan=2, columnspan=1, sticky=W + E + N + S)

        self.bCancel = Button(self, text='Cancel', command=lambda: self.close_window())
        self.bCancel.grid(row=10, column=5, rowspan=2, columnspan=1, sticky=W + E + N + S)

    def update_param(self):
        for y in self.params_names:
            # controllo se il valore del fix sia pari a 1
            # in caso pongo min=max=sv
            i = 0
            for x in self.params_attribute:
                if x == 'fix':
                    self.param_dict[y][x] = eval("self.{}[-1].get()".format(y))
                else:
                    if eval("self.{}[-1].get()".format(y)) == 0:
                        self.param_dict[y][x] = eval("self.{}[i].get()".format(y))
                        i = i + 1
                    else:
                        self.param_dict[y][x] = eval("self.{}[0].get()".format(y))

    def createInsertRange(self, var, row_nr, frame, par_name):
        Label(frame, height=1, width=5, text=par_name).grid(row=row_nr, column=0)

        p0 = Entry(frame, textvariable=var[0])
        p0.grid(row=row_nr, column=1)
        p0.config(width=7)

        p2 = Entry(frame, textvariable=var[1])
        p2.grid(row=row_nr, column=2)
        p2.config(width=7)

        p3 = Entry(frame, textvariable=var[2])
        p3.grid(row=row_nr, column=3)
        p3.config(width=7)

        p4 = Checkbutton(frame, variable=var[3])
        p4.grid(row=row_nr, column=4)

    def create_range_parameter(self, frame, params_names, params_attribute):
        i00 = Label(frame, height=1, width=10, text="").grid(row=0, column=0)
        i01 = Label(frame, height=1, width=10, text="N sample").grid(row=1, column=0)

        self.nTime = IntVar()
        nTimeEntry = Entry(frame, textvariable=self.nTime)
        nTimeEntry.grid(row=1, column=1)
        nTimeEntry.config(width=7)
        self.nTime.set('1')

        row_idx = 1
        ### create table
        i0 = Label(frame, height=1, width=10, text="").grid(row=row_idx + 1, column=0)
        i1 = Label(frame, height=1, width=10, text="Parameter").grid(row=row_idx + 2, column=0)
        i1 = Label(frame, height=1, width=10, text="Starting value").grid(row=row_idx + 2, column=1)
        i2 = Label(frame, height=1, width=10, text="Min").grid(row=row_idx + 2, column=2)
        i3 = Label(frame, height=1, width=10, text="Max").grid(row=row_idx + 2, column=3)
        i4 = Label(frame, height=1, width=10, text="Fixed").grid(row=row_idx + 2, column=4)

        # l'ordine di questi due array definisce l'ordine con cui si interroga il dizionario
        row_idx = row_idx + 3

        for y in params_names:
            # definisco un array con il nome
            setattr(self, y, [])
            for x in params_attribute:
                if x == 'fix':
                    eval("self.{}.append(IntVar())".format(y))
                    eval("self.{}[-1].set(self.param_dict[y][x])".format(y))
                else:
                    eval("self.{}.append(StringVar())".format(y))
                    eval("self.{}[-1].set(self.param_dict[y][x])".format(y))

            self.createInsertRange(eval("self.{}".format(y)), row_idx, frame, y)
            row_idx = row_idx + 1


    def list_combobox_ts(self, dictObject, tableObject):

        if self.NameTS.get() != '':
            # lettura delle date della TS selezionata
            tmp = tableObject.loc[tableObject.Name == self.NameTS.get(), 'keys'].values[0]
            v = np.array([isinstance(dictObject[tmp][0].values[i], datetime.datetime) for i in
                          xrange(dictObject[tmp][0].values.__len__())])

            # creazione della lista con le date disponibili
            date_list_tmp = dictObject[tmp][0].values[np.array(v)]
            date_list = [date_list_tmp[i].strftime("%d/%m/%Y") for i in xrange(date_list_tmp.__len__())]

            # selezione della data minima
            self.te_dateMin['values'] = date_list
            self.te_dateMin.current(0)

            # selezione della data massima
            date_list_tmp[::-1].sort()
            date_list = [date_list_tmp[i].strftime("%d/%m/%Y") for i in xrange(date_list_tmp.__len__())]
            self.te_dateMax['values'] = date_list
            self.te_dateMax.current(0)


    def go_to_calib(self, dictObject, tableObject):
        # recupero la calibrazione selezionata
        self.update_param()
        print 'Letto tutte le configurazioni'

        if self.set_mkt_ts.get() == 'MKT':
            print 'Calibratore per il mercato'

            if self.mkt_calibration_type.get() == 'CURVE':
                if self.NameCurve.get() == "":
                    tkMessageBox.showinfo("Error", "Nessuna curva selezionata.")
                else:
                    tmp = tableObject.loc[tableObject.Name == self.NameCurve.get(), 'keys'].values[0]
                    self.CurveChosen = dictObject[tmp]
                    self.master.destroy()
            elif self.mkt_calibration_type.get() == 'OPT':
                if self.NameOption.get() == "":
                    tkMessageBox.showinfo("Error", "Nessuna opzione selezionata.")
                else:
                    tmp = tableObject.loc[tableObject.Name == self.NameOption.get(), 'keys'].values[0]
                    self.OptionChosen = dictObject[tmp]
                    self.master.destroy()
            else:
                if self.NameCurve.get() == "":
                    tkMessageBox.showinfo("Error", "Nessuna curva selezionata.")
                elif self.NameOption.get() == "":
                    tkMessageBox.showinfo("Error", "Nessuna opzione selezionata.")
                else:
                    tmp_curve = tableObject.loc[tableObject.Name == self.NameCurve.get(), 'keys'].values[0]
                    self.CurveChosen = dictObject[tmp_curve]
                    tmp_opt = tableObject.loc[tableObject.Name == self.NameOption.get(), 'keys'].values[0]
                    self.OptionChosen = dictObject[tmp_opt]
                    self.master.destroy()


        if self.set_mkt_ts.get() == 'TS':
            print 'Calibratore per il TS'

            if self.NameTS.get() == "":
                tkMessageBox.showinfo("Error", "Nessuna Time Series selezionata.")
            else:
                tmp = tableObject.loc[tableObject.Name == self.NameTS.get(), 'keys'].values[0]
                self.TSChosen = dictObject[tmp]

                self.TS_dateMIN = datetime.datetime.strptime(self.te_dateMin.get(), "%d/%m/%Y")
                self.TS_dateMAX = datetime.datetime.strptime(self.te_dateMax.get(), "%d/%m/%Y")

                self.master.destroy()


def readSheetObject_clearRows(input_dict = None):

    k = 0
    obj_sheet = {}
    indicator_modified = False
    for key_considering in input_dict.keys():
        df = input_dict[key_considering]
        tmp = df.isna().sum(1)
        tmp = tmp[tmp == df.shape[1]].index.__array__()
        i = 0

        if tmp.size == 0 :
            k = k +1
            obj_sheet[k] = df
        elif tmp[0] > 0:
            tmp = np.insert(tmp,0,-1)

        for idx_row in tmp:
            indicator_modified = True
            if i == (tmp.shape[0] - 1):
                df_new = df.loc[(idx_row + 1):df.shape[0], :]
                if not df_new.empty:
                    k = k + 1
                    df_new.reset_index(drop=True,inplace= True)
                    obj_sheet[k] = df_new

            elif (idx_row + 1) != tmp[i + 1]:
                df_new = df.loc[(idx_row + 1):(tmp[i + 1] - 1), :]
                if not df_new.empty:
                    k = k + 1
                    df_new.reset_index(drop=True,inplace= True)
                    obj_sheet[k] = df_new

            i = i + 1

    return obj_sheet , indicator_modified


def readSheetObject_clearColumns(input_dict = None):

    k = 0
    obj_sheet = {}
    indicator_modified = False
    for key_considering in input_dict.keys():
        df = input_dict[key_considering]
        tmp = df.isna().sum(0)
        tmp = tmp[tmp == df.shape[0]].index.__array__()
        i = 0

        if tmp.size == 0 :
            k = k +1
            obj_sheet[k] = df
        elif tmp[0] > 0:
            tmp = np.insert(tmp,0,-1)

        for idx_col in tmp:
            indicator_modified = True
            if i == (tmp.shape[0] - 1):
                df_new = df.loc[:, (idx_col + 1):df.shape[0]]
                if not df_new.empty:
                    k = k + 1
                    df_new.columns = np.arange(0, df_new.shape[1])
                    obj_sheet[k] = df_new

            elif (idx_col + 1) != tmp[i + 1]:
                df_new = df.loc[:, (idx_col + 1):(tmp[i + 1] - 1)]
                if not df_new.empty:
                    k = k + 1
                    df_new.columns = np.arange(0, df_new.shape[1])
                    obj_sheet[k] = df_new

            i = i + 1

    return obj_sheet , indicator_modified


def readSheetObject(workbook_path,sheet_name):

    # carico il foglio
    df = pd.read_excel(workbook_path,sheet_name,header = None)

    # estraggo tutti gli oggetti presenti nel foglio
    obj_sheet = {1: df}
    obj_sheet , indicator_row_work = readSheetObject_clearRows(obj_sheet)
    obj_sheet , indicator_col_work = readSheetObject_clearColumns(obj_sheet)

    while indicator_row_work or indicator_col_work:
        obj_sheet, indicator_row_work = readSheetObject_clearRows(obj_sheet)
        obj_sheet, indicator_col_work = readSheetObject_clearColumns(obj_sheet)

    return obj_sheet


def readFeaturesObject(input_dict):
    col = ['keys', 'TypeObject', 'Name']
    element_on_sheet = pd.DataFrame(columns=col)

    for k in input_dict.keys():
        item = input_dict[k]
        if u'CurveType' in item.loc[:, 0].values:
            element_on_sheet = element_on_sheet.append({'keys': k,
                                                        'TypeObject': 'Curve',
                                                        'Name': item.loc[0, 0]}, ignore_index=True)

        elif u'OptionType' in item.loc[:, 0].values:
            element_on_sheet = element_on_sheet.append({'keys': k,
                                                        'TypeObject': 'Option',
                                                        'Name': item.loc[0, 0]}, ignore_index=True)

        elif u'TSType' in item.loc[:, 0].values:
            element_on_sheet = element_on_sheet.append({'keys': k,
                                                        'TypeObject': 'TS',
                                                        'Name': item.loc[0, 0]}, ignore_index=True)

        elif u'MatrixType' in item.loc[:, 0].values:
            element_on_sheet = element_on_sheet.append({'keys': k,
                                                        'TypeObject': 'Matrix',
                                                        'Name': item.loc[0, 0]}, ignore_index=True)



    return element_on_sheet