from pyxll import xlcAlert
from win32com.client import constants as const
from W_calibration import readSheetObject, readSheetObject_clearColumns, readSheetObject_clearRows

import pandas as pd
import numpy as np
from Tkinter import *
import tkMessageBox
import ttk
from DEF_intef import nameSheetBootstrap


def popup_messagebox(msg):
    xlcAlert(msg)
#----
def drawBox(xla, spessore , rTopLeft = 0, cTopLeft = 0,rBottomRight=0,cBottomRight=0, Colore=0):
    if (rTopLeft <= 0) or(cTopLeft <= 0) or (rBottomRight <= 0)or (cBottomRight <= 0):
        msg = "Le coordinate del box devono essere maggiori di zero"
        popup_messagebox(msg)
        sys.exit()

    RR = xla.Range(xla.Cells(rTopLeft, cTopLeft), xla.Cells(rBottomRight, cBottomRight))
    RR.Borders(const.xlDiagonalDown).LineStyle = const.xlNone
    RR.Borders(const.xlDiagonalUp).LineStyle = const.xlNone
    RR.Borders(const.xlEdgeLeft).LineStyle = const.xlContinuous
    RR.Borders(const.xlEdgeLeft).Weight = spessore
    RR.Borders(const.xlEdgeLeft).ColorIndex = Colore

    RR.Borders(const.xlEdgeTop).LineStyle = const.xlContinuous
    RR.Borders(const.xlEdgeTop).Weight = spessore
    RR.Borders(const.xlEdgeTop).ColorIndex = Colore

    RR.Borders(const.xlEdgeBottom).LineStyle = const.xlContinuous
    RR.Borders(const.xlEdgeBottom).Weight = spessore
    RR.Borders(const.xlEdgeBottom).ColorIndex = Colore

    RR.Borders(const.xlEdgeRight).LineStyle = const.xlContinuous
    RR.Borders(const.xlEdgeRight).Weight = spessore
    RR.Borders(const.xlEdgeRight).ColorIndex = Colore


def drawLine (xla, rTopLeft, cTopLeft, rBottomRight, cBottomRight, hor, spessore):
    (xla.Range(xla.Cells(rTopLeft, cTopLeft), xla.Cells(rBottomRight, cBottomRight))).Select()
    if (hor =="o"):
        xla.Selection.Borders(const.xlEdgeBottom).LineStyle = const.xlContinuous
        xla.Selection.Borders(const.xlEdgeBottom).Weight = spessore
        xla.Selection.Borders(const.xlEdgeBottom).ColorIndex = const.xlAutomatic

    else:
        xla.Selection.Borders(const.xlEdgeRight).LineStyle = const.xlContinuous
        xla.Selection.Borders(const.xlEdgeRight).Weight = spessore
        xla.Selection.Borders(const.xlEdgeRight).ColorIndex = const.xlAutomatic
#----------
def formatTestataCurva(xla,nRiga,nColonna,nLarghezza,testo):

    (xla.Range(xla.Cells(nRiga, nColonna), xla.Cells(nRiga, nColonna + nLarghezza - 1))).Select()
    xla.Selection.HorizontalAlignment   = const.xlCenter
    xla.Selection.VerticalAlignment     = const.xlBottom
    xla.Selection.WrapText              = False
    xla.Selection.Orientation           = 0
    xla.Selection.AddIndent             = False
    xla.Selection.IndentLevel           = 0
    xla.Selection.ShrinkToFit           = False
    xla.Selection.ReadingOrder          = const.xlContext
    xla.Selection.MergeCells            = False
    xla.Selection.Merge()
    xla.Selection.Font.ColorIndex       = 2
    xla.Selection.Font.Bold             = True
    xla.Selection.Interior.ColorIndex   = 55
    xla.Selection.Interior.Pattern      = const.xlSolid
    xla.Selection.Value                 = testo

"""
def formatTestataCurvaCDS(xla,nRiga,nColonna,nLarghezza,testo, idx_elab):

    (xla.Range(xla.Cells(nRiga, nColonna - (idx_elab - 1)), xla.Cells(nRiga, nColonna - (idx_elab - 1)))).Select()
    
    #xla.Cells(nRiga, nColonna + nLarghezza - 1).Value = 'XXX'
    #xla.Cells(nRiga, nColonna + nLarghezza - 1).Value = 'XXX'
    
    xla.Selection.HorizontalAlignment   = const.xlCenter
    xla.Selection.VerticalAlignment     = const.xlBottom
    xla.Selection.WrapText              = False
    xla.Selection.Orientation           = 0
    xla.Selection.AddIndent             = False
    xla.Selection.IndentLevel           = 0
    xla.Selection.ShrinkToFit           = False
    xla.Selection.ReadingOrder          = const.xlContext
    xla.Selection.MergeCells            = False
    xla.Selection.Merge()
    xla.Selection.Font.ColorIndex       = 2
    xla.Selection.Font.Bold             = True
    xla.Selection.Interior.ColorIndex   = 55
    xla.Selection.Interior.Pattern      = const.xlSolid
    xla.Selection.Value                 = testo
    
    (xla.Range(xla.Cells(nRiga, nColonna), xla.Cells(nRiga, nColonna + nLarghezza - 1))).Select()
"""

#----------

def readCurvesNames(xla, s, rangeStart, direzione, distanza, offset = 0):

    r = xla.Range(rangeStart)
    curveL = []

    if direzione.lower() == "o":
        i = 0
        while (r.Value != None):
            i += 1
            if offset > 0:
                r = xla.Range(xla.Cells(row + offset, col + 1), xla.Cells(row + offset, col + 1))
            nomeCurva = r.Value
            curveL.append((nomeCurva,i))
            nCols     = r.Columns.Count
            row       = r.Row
            col       = r.Column
            r         = xla.Range(xla.Cells(row, col + distanza), xla.Cells(row, col + (nCols - 1) + distanza))
    else:
        j = 1
        nomeCurva = r.Value

        if offset > 0:
            row = r.Row
            col = r.Column
            r_tmp = xla.Range(xla.Cells(row + offset, col + 1), xla.Cells(row + offset, col + 1))
            nomeCurva = r_tmp.Value

        curveL.append((nomeCurva, j))

        while (r.Value != None):
            nCols = r.Columns.Count
            row   = r.Row
            col   = r.Column
            r = xla.Range(xla.Cells(row + 1, col), xla.Cells(row + 1, col))
            # porto avanti ancora per controllare

            if (r.Value == None):
                
                
                r = xla.Range(xla.Cells(row  + distanza, col), xla.Cells(row + distanza, col))
                
                #print 'r.Value C: ', r.Value

                if r.Value != None:
                    j += 1
                    nomeCurva = r.Value
                    # ----
                    if offset > 0:
                        r_tmp       = xla.Range(xla.Cells(row + offset + distanza, col + 1), xla.Cells(row + offset +  distanza, col + 1))
                        nomeCurva   = r_tmp.Value
                    # ----
                    curveL.append((nomeCurva, j))

    return curveL



def readCurvesParmsNames (xla, s, rangeStart):
    r = xla.Range(rangeStart)
    curveDict = {}
    # --- inizializzo la prima
    nomeCurva = r.Value
    j = 1
    nomeCurva = nomeCurva+"-"+str(j)
    curveDict[nomeCurva] = []
    rp = xla.Range(xla.Cells(r.Row, "I"), xla.Cells(r.Row, "I"))
    p_i = 0
    while (rp.Value != None):
        print "RP address", rp.Row, rp.Column
        p_i += 1
        nome = rp.Value
        if nome[0:3] == "PIL": distanzaO = 4
        else:                  distanzaO = 7
        curveDict[nomeCurva].append((nome, p_i))
        row = rp.Row
        col = rp.Column
        rp = xla.Range(xla.Cells(row, col + distanzaO), xla.Cells(row, col + distanzaO))
    # --- proseguo con le successive
    while (r.Value != None):
        print "loop su R", r.Value
        r = xla.Range(xla.Cells(r.Row + 1, r.Column), xla.Cells(r.Row + 1, r.Column))
        if (r.Value == None):
            r = xla.Range(xla.Cells(r.Row + 1, r.Column), xla.Cells(r.Row + 1, r.Column))

            if r.Value != None: # ho trovato una nuova curva
                j += 1
                nomeCurva = r.Value
                nomeCurva = nomeCurva + "-" + str(j)
                curveDict[nomeCurva] = []
                # --- prendo i relativi parametri disponibili
                rp  = xla.Range(xla.Cells(r.Row, "I"), xla.Cells(r.Row, "I"))
                p_i = 0
                while (rp.Value != None):
                    p_i += 1
                    nome = rp.Value
                    curveDict[nomeCurva].append((nome, p_i))
                    row  = rp.Row
                    col  = rp.Column
                    if nome[0:3] == "PIL": distanzaO = 4
                    else                 : distanzaO = 7
                    rp = xla.Range(xla.Cells(row, col + distanzaO), xla.Cells(row, col + distanzaO))

    print curveDict
    return curveDict


def findRigthPlaceBootCurveSeg(xla, r, distCurve, dir="O"):
    rOut = None
    if dir == "v" :
        if (r.Value == None): return r
        nCols = r.Columns.Count
        row   = r.Row
        col   = r.Column
        j = 0
        while (r.Value != None):
            j += 1
            r = xla.Range(xla.Cells(row + j, col), xla.Cells(row + j, col))
            #porto avanti ancora per controllare
            if (r.Value == None):  r = xla.Range(xla.Cells(row + j + distCurve , col), xla.Cells(row + j+distCurve, col))
        r = xla.Range(xla.Cells(row + j + distCurve, col), xla.Cells(row+j+distCurve, col))

    else:
        while (r.Value != None):
            nCols = r.Columns.Count
            row =r.Row
            col =r.Column
            #r = xla.Range(xla.Cells(row, col + distCurve), xla.Cells(row, col + (nCols-1) + distCurve))
            r = xla.Range(xla.Cells(row, col + distCurve), xla.Cells(row, col + distCurve))

    rOut = r
    if (rOut == None):
        msg = "Unable to compute the output range for your curve"
        print msg
        sys.exit()

    return rOut


def findRigthPlaceBootCurveSeg_m(xla, r, distCurve, dir="O"):
    rOut = None
    if dir == "v" :
        if (r.Value == None): return r
        nCols = r.Columns.Count
        row   = r.Row
        col   = r.Column
        j = 0
        while (r.Value != None):
            j += 1
            r = xla.Range(xla.Cells(row + j, col), xla.Cells(row + j, col))
            #porto avanti ancora per controllare
            if (r.Value == None):  r = xla.Range(xla.Cells(row + j + distCurve , col), xla.Cells(row + j+distCurve, col))
        r = xla.Range(xla.Cells(row + j + distCurve, col), xla.Cells(row+j+distCurve, col))

    else:
        while (r.Value != None):
            nCols = r.Columns.Count
            row =r.Row
            col =r.Column
            r = xla.Range(xla.Cells(row, col + distCurve), xla.Cells(row, col + (nCols-1) + distCurve))
            #r = xla.Range(xla.Cells(row, col + distCurve), xla.Cells(row, col + distCurve))

    rOut = r
    if (rOut == None):
        msg = "Unable to compute the output range for your curve"
        print msg
        sys.exit()

    return rOut


def allSheet(wb):
    names = []
    for s in wb.Worksheets:
        names.append(s.Name)

    return names


def findCalibrationPos(xla, nameSheet):
    rangeInitial = "B2"
    distanza = 2
    sheet = xla.ActiveWorkbook.Sheets(nameSheet)
    r = sheet.Range(rangeInitial)
    col = r.Column
    i = 0

    rStart = xla.Range(xla.Cells(r.Row + distanza, col), xla.Cells(r.Row + distanza, col))

    while r.Value != None or rStart.Value != None:
        i = 1
        r = xla.Range(xla.Cells(r.Row + 1, col), xla.Cells(r.Row + 1, col))
        rStart = xla.Range(xla.Cells(r.Row + distanza, col), xla.Cells(r.Row + distanza, col))

    if i == 0:
        rStart = sheet.Range(rangeInitial)

    return rStart


def writeResultPandas( xla, rng , df, flagPrintColumns = True):

    nRows           = df.shape[0]
    nCols           = df.shape[1]
    topLeftRow      = rng.Row
    topLeftCol      = rng.Column
    drawBox(xla, 3 , topLeftRow, topLeftCol,topLeftRow + nRows - (1 - int(flagPrintColumns)), topLeftCol + nCols - 1, 0)

    if flagPrintColumns == True:
        for j in xrange(0,nCols):
            xla.Cells(topLeftRow, topLeftCol + j).Font.Bold = True
            xla.Cells(topLeftRow, topLeftCol + j).HorizontalAlignment = const.xlCenter
            xla.Cells(topLeftRow , topLeftCol + j).Value = df.columns.values[j]

        # scrittura della linea sotto i nomi delle colonne
        RR = xla.Range(xla.Cells(topLeftRow, topLeftCol), xla.Cells(topLeftRow, topLeftCol + nCols - 1))
        RR.Borders(const.xlEdgeBottom).LineStyle = const.xlContinuous
        RR.Borders(const.xlEdgeBottom).Weight = const.xlThin
        RR.Borders(const.xlEdgeBottom).ColorIndex = 0

        topLeftRow = topLeftRow + 1


    for i in xrange(0,nRows):
        for j in xrange(0,nCols):
            xla.Cells(topLeftRow + i, topLeftCol+j).HorizontalAlignment = const.xlCenter
            xla.Cells(topLeftRow + i, topLeftCol+j).Value   = df.iloc[i,j]

    rangeStart = xla.Range(xla.Cells(topLeftRow + nRows + 1, topLeftCol),xla.Cells(topLeftRow + nRows + 1, topLeftCol))
    return rangeStart


def readCurveDiscoutFitting(input_dict):
   col = ['keys', 'TypeObject', 'Name']
   element_on_sheet = pd.DataFrame(columns=col)

   for k in input_dict.keys():
      item = input_dict[k]
      if u'Discount Factors' in item.loc[:,1].values:
          element_on_sheet=element_on_sheet.append({'keys': k,
                                                    'TypeObject': 'Discount Curve',
                                                    'Name': item.loc[0,0]}, ignore_index=True)
      if u'Interp. Model' in item.loc[:, 0].values:
          element_on_sheet = element_on_sheet.append({'keys': k,
                                                      'TypeObject': 'Fitting',
                                                      'Name': item.loc[0, 0]}, ignore_index=True)

   return element_on_sheet


class ChooseAvaiableBootstrapCurve(Frame):

    def close_window(self):
        self.res_disc = self.curve_disc.get()
        self.res_parm = self.curve_param.get()

        self.master.destroy()

    def close_without_selection(self):
        self.res_disc = None
        self.res_parm = None

        self.master.destroy()

    def __init__(self, master= None , discount_curves = None, prm = None):

        Frame.__init__(self, master)

        self.master = master
        self.master.title('Curve data selection')

        Label(self.master, text='Curve').grid(row=1, column=0, rowspan=1, columnspan=1, pady=2, sticky=W + E + N + S)

        self.curve_disc = ttk.Combobox(self.master)
        CurveDisc = StringVar()
        self.curve_disc.config(textvariable=CurveDisc, state="readonly", values=discount_curves)
        self.curve_disc.grid(row=1, column=1, rowspan=1, columnspan=1, pady=2, sticky=W + E + N + S)

        Label(self.master, text='FittingParams').grid(row=2, column=0, rowspan=1, columnspan=1, pady=2, sticky=W + E + N + S)

        self.curve_param = ttk.Combobox(self.master)
        CurveParam = StringVar()
        self.curve_param.config(textvariable=CurveParam, state="readonly", values=prm)
        self.curve_param.grid(row=2, column=1, rowspan=1, columnspan=1, pady=2, sticky=W + E + N + S)

        # Bottoni
        SubmitButton = ttk.Button(self.master, text='Submit', command = self.close_window)
        SubmitButton.grid(row=3, column=0, rowspan=1, columnspan=1, pady=2, sticky=W + E + N + S)

        CancelButton = ttk.Button(self.master, text='Cancel', command = self.close_without_selection)
        CancelButton.grid(row=3, column=1, rowspan=1, columnspan=1, pady=2, sticky=W + E + N + S)


def CurveBootstrapedFromXls(book):

    # leggo gli oggetti sul foglio in input
    try:
        avaiable_curve = readSheetObject(book, nameSheetBootstrap)
    except:
        root = Tk()
        root.withdraw()
        msg0 = "Non e' presente il foglio: {} !".format(nameSheetBootstrap)
        tkMessageBox.showinfo("Attenzione!", msg0)
        root.destroy()
        return

    # leggo i nomi degli oggeti Discout e Fitting
    obj_curve = readCurveDiscoutFitting(avaiable_curve)

    # lista delle possibili scelte
    disc_curves_choices = obj_curve.loc[obj_curve.loc[:,'TypeObject'] == 'Discount Curve','Name'].tolist()
    params_choices      = obj_curve.loc[obj_curve.loc[:,'TypeObject'] == 'Fitting','Name'].tolist()


    # scelta dell'utente tra la lista dei fattori di sconto e dei paramentri
    root=Tk()
    choices = ChooseAvaiableBootstrapCurve(master = root, discount_curves = disc_curves_choices, prm = params_choices)
    root.mainloop()

    if choices.res_disc == None:
        root=Tk()
        tkMessageBox.showinfo('Salutation', 'Au revoir')
        root.destroy()
        return

    # lettura della curva selezionata
    sel_df  = int(obj_curve.loc[obj_curve.loc[:,'Name'] == choices.res_disc,'keys'])

    # lettura della data della curva
    df_date_ref = avaiable_curve[sel_df].loc[(avaiable_curve[sel_df].loc[:,0]=='Date Ref'),1].values[0]

    # selezione dei nodi solo con la Y
    curve = avaiable_curve[sel_df].loc[(avaiable_curve[sel_df].loc[:, 2] == 'Y'), [0, 1]]

    # dizionario di output dei valori di interesse
    output = {}
    output['curve_dates'] = curve.loc[:, 0].dt.date.tolist()
    output['curve_df_val'] = curve.loc[:, 1].astype(float).tolist()

    if choices.res_parm != '':

        # lettura dei parametri
        sel_prm = int(obj_curve.loc[obj_curve.loc[:, 'Name'] == choices.res_parm, 'keys'])
        obj_param = avaiable_curve[sel_prm]

        # lettura della data della curva
        prm_date_ref = obj_param.loc[(obj_param.loc[:, 0] == 'Date Ref'), 1].values[0]

        # verifica che la data di riferimento sia uguale tra i parametri e i fattori di sconto
        if prm_date_ref != df_date_ref:
            root = Tk()
            tkMessageBox.showinfo('Warning', "La data di riferimento della curva dei fattori di sconto e' diversa rispetto a quella dei parametri.")
            root.destroy()
            return

        # elimino le intestazioni: in base alle date o ai tempi
        row_sep = obj_param.loc[obj_param.loc[:, 0] == 'Date',].index[0]
        colnames = obj_param.iloc[row_sep, :].tolist()

        # elimino anche il primo valore che corrisponde alla data di riferimento # caso solo del modello LINEARE ?!
        obj_param = obj_param.loc[(row_sep + 2):obj_param.shape[0], :]
        obj_param.columns = colnames

        out_param ={}
        out_param['Dates'] = obj_param.loc[:,'Date'].dt.date.tolist()
        out_param['a']     = []
        out_param['b']     = []

        a = obj_param.loc[:,'a'].values
        b = obj_param.loc[:,'b'].values

        # metto nello stesso formato del Database
        # il formato di questi due vettori e' veramente ORRIBILE

        for i in xrange(obj_param.shape[0]):
            out_param['a'].append(np.array([a[i]]))
            out_param['b'].append(np.array([b[i]]))

        output['date_ref']     = prm_date_ref.date()
        output['params_model'] = 'LIN'
        output['params']       = out_param

    return output

