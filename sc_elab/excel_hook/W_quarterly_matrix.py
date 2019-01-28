import numpy  as np
import pandas as pd
import ttk
from sc_elab.excel_hook.xls_utils import allSheet, findCalibrationPos, drawBox, formatTestataCurva, writeResultPandas
from win32com.client import constants as const

from Tkinter import *

class W_dim_matrix(Frame):
    def __init__(self, master = None):

        Frame.__init__(self, master)
        self.master = master
        self.master.title("Matrix Dimension")


        Label(self.master, text='Write the dimension of the matrix template').grid(row = 0, column = 0, rowspan = 1, columnspan = 3, sticky = W+E+N+S)

        self.dimMatrix = IntVar()
        self.dimMatrix.set(4)
        Entry(self.master,textvariable = self.dimMatrix,width=7,justify='center').grid(row = 1, column = 0, rowspan = 1, columnspan = 3, sticky = W+E+N+S,padx=10,pady=(20,2))

        self.butt_proceed = Button(self.master, text="Go",command = lambda: self.elaborate_request()).grid(row=10, column=0,columnspan = 3,sticky = W+E,padx=20,pady=(20,2))
        self.butt_cancel  = Button(self.master, text="Cancel",command = lambda: self.close_window()).grid(row=11, column=0,columnspan = 3,sticky = W+E,padx=20,pady=2)

    def elaborate_request(self):
        self.master.destroy()

    def close_window(self):
        self.master.destroy()



class W_select_matrix(Frame):
    def __init__(self, master = None , ListMatrix = []):

        Frame.__init__(self, master)
        self.master = master
        self.master.title("Matrix Selection")


        Label(self.master, text='Avaiable matrix').grid(row = 0, column = 0, rowspan = 1, columnspan = 3, sticky = W+E+N+S)

        self.NameMatrix = StringVar()

        self.cb_matrix_avaible = ttk.Combobox(self.master,textvariable = self.NameMatrix, state = "readonly", values = ListMatrix)
        self.cb_matrix_avaible.grid(row = 1, column = 0, rowspan = 1, columnspan = 3, sticky = W+E+N+S,padx=10,pady=(20,2))

        self.butt_proceed = Button(self.master, text="Go",command = lambda: self.elaborate_request()).grid(row=10, column=0,columnspan = 3,sticky = W+E,padx=20,pady=(20,2))
        self.butt_cancel  = Button(self.master, text="Cancel",command = lambda: self.close_window()).grid(row=11, column=0,columnspan = 3,sticky = W+E,padx=20,pady=2)

    def elaborate_request(self):
        self.MatrixNameChosen = self.NameMatrix.get()
        self.master.destroy()

    def close_window(self):
        self.master.destroy()


class W_matrix_trim(Frame):

    def __init__(self, master = None , matrix = None):

        Frame.__init__(self, master)
        self.grid()
        self.master.title("Matrix Elaborate")

        for r in xrange(1):
            self.master.rowconfigure(r, weight=1)
        for c in xrange(5):
            self.master.columnconfigure(c, weight=1)

        # self.matrix_input
        self.matrix_input = matrix.values
        self.nRows, self.nCols = matrix.shape

        ############## Definizione dello spazio master
        FrameInput = Frame(master)
        FrameInput.grid(row = 0, column = 0, rowspan = 2, columnspan = 1, sticky = W+E+N+S)
        Label(FrameInput, height=1, width=20, text='INPUT MATRIX',font=('Arial', 10, 'bold')).grid(row=0, column=0, columnspan = self.nCols,sticky = W)

        ttk.Separator(master,orient="vertical").grid(row=0, column=1,sticky="ns")

        FrameParam = Frame(master)
        FrameParam.grid(row = 0, column = 2, rowspan = 2, columnspan = 1, sticky = W+E+N+S)
        Label(FrameParam, height=1, width=20, text='OPTIONS',font=('Arial', 10, 'bold')).grid(row=0, column=0, columnspan = 2,sticky = W)

        ttk.Separator(master,orient="vertical").grid(row=0, column=3,sticky="ns")

        FrameOutput = Frame(master)
        FrameOutput.grid(row = 0, column = 4, rowspan = 2, columnspan = 1, sticky = W+E+N+S)
        Label(FrameOutput, height=1, width=20, text='OUTPUT MATRIX',font=('Arial', 10, 'bold')).grid(row=0, column=0, columnspan = self.nCols,sticky = W)


        ############# INPUT
        # Matrice di transizione di input
        r,c = self.create_matrix_range(frame = FrameInput, row_idx = 1,col_idx = 0, matrix_print = self.matrix_input, name_variable = 'IN_ATM', title = 'Annual Transition Matrix', state = DISABLED)

        # Autovalori della matrice di transizione di input
        self.M_diag_iput = np.diag(np.linalg.eig(self.matrix_input)[0])
        r,c = self.create_matrix_range(frame = FrameInput, row_idx = r+1,col_idx = 0, matrix_print = self.M_diag_iput, name_variable = 'IN_ATM_E', title = 'Eigenvalues - Annual Transition Matrix', state = NORMAL)

        ############# PARAMETRI
        self.ck1 = IntVar()
        Checkbutton(FrameParam,justify='left',text='Fix Eigenvalues positive in Input matrix',variable = self.ck1,anchor="w",state=DISABLED).grid(row=2, column=0,columnspan = 1,pady=(10,2),sticky = W+E)
        self.ck1.set(1)

        self.ck2 = IntVar()
        Checkbutton(FrameParam,justify='left',text='Floor at 0 for Quarterly Transition Matrix',variable = self.ck2, anchor="w").grid(row=3, column=0,columnspan = 1,pady=2,sticky = W+E)
        self.ck2.set(1)

        self.ck3 = IntVar()
        Checkbutton(FrameParam,justify='left',text='Total amount each row sum equal to 100% \n for Quarterly Transition Matrix',variable = self.ck3, anchor="w").grid(row=4, column=0,columnspan = 1,pady=2,sticky = W+E)
        self.ck3.set(1)


        Label(FrameParam, text="Confidence (%) between Annual Output matrix \n and Input matrix ", justify='left').grid(row=5, column=0,pady=(20,10),sticky = W+E)
        self.ck4 = DoubleVar()
        self.ck4.set(1.0)
        Entry(FrameParam,textvariable = self.ck4,width=7,justify='center').grid(row=5, column=1,columnspan = 1,pady=2,sticky = W+E)


        self.butt_proceed = Button(FrameParam, text="Elaborate",command = lambda: self.elaborate_request()).grid(row=10, column=0,columnspan = 2,sticky = W+E,padx=10,pady=(20,2))
        self.butt_cancel  = Button(FrameParam, text="Cancel",command = lambda: self.close_window(False)).grid(row=11, column=0,columnspan = 2,sticky = W+E,padx=10,pady=2)
        self.butt_return  = Button(FrameParam, text="Save and Quit",command = lambda: self.close_window(True)).grid(row=12, column=0,columnspan = 2,sticky = W+E,padx=10,pady=2)

        ############# OUTPUT
        #self.M_out_trim = np.zeros(self.matrix_input.shape)
        self.update_matrix_output(np.zeros(self.matrix_input.shape))

        # Matrice di transizione di output - Trimestrale
        r,c = self.create_matrix_range(frame = FrameOutput, row_idx = 1,col_idx = 0, matrix_print = self.M_out_trim, name_variable = 'QTM', title = 'Quarterly Transition Matrix',state = DISABLED)
        # Autovalori della Matrice di transizione di output - Trimestrale
        r,c_fix = self.create_matrix_range(frame = FrameOutput, row_idx = r+1,col_idx = 0, matrix_print = self.M_out_trim_diag, name_variable = 'QTM_E', title = 'Eigenvalues - Quarterly Transition Matrix',state = DISABLED)

        # Matrice di transizione di output - Annuale
        r,c = self.create_matrix_range(frame = FrameOutput, row_idx = 1,col_idx = c_fix, matrix_print = self.M_out_ann, name_variable = 'ATM', title = 'Recalculated Annual Transition Matrix',state = NORMAL)
        # Autovalori della Matrice di transizione di output - Annuale
        r,c = self.create_matrix_range(frame = FrameOutput, row_idx = r+1,col_idx = c_fix, matrix_print = self.M_out_ann_diag, name_variable = 'ATM_E', title = 'Eigenvalues - Recalculated Annual Transition Matrix',state = DISABLED)



    def update_matrix_output(self, m):
        self.M_out_trim = m.copy()
        self.M_out_trim_diag = np.diag(np.linalg.eig(self.M_out_trim)[0])
        self.M_out_ann = np.dot(self.M_out_trim, np.dot(self.M_out_trim, np.dot(self.M_out_trim,self.M_out_trim)))
        self.M_out_ann_diag = np.diag(np.linalg.eig(self.M_out_ann)[0])



    def create_matrix_range(self, frame, row_idx, col_idx, matrix_print, name_variable, title, state):
        Label(frame, height=1, width=7, text="").grid(row=row_idx, column=col_idx)
        Label(frame, height=1, width=30, text=title).grid(row=row_idx+ 1, column=col_idx+1, columnspan = self.nCols)

        row_idx = row_idx + 2
        col_idx = col_idx + 1
        for i in xrange(0,self.nRows):
            for j in xrange(0, self.nCols):
                y = name_variable + str(i) + str(j)
                exec("self.{} = StringVar()".format(y))
                eval("self.{}.set(str(round(100 * matrix_print[i,j],2)))".format(y))
                x = 'E_' + name_variable + str(i) + str(j)
                exec("self.{} = Entry(frame, textvariable= self.{} ,width=7,justify='center',state=state)".format(x,y))
                if matrix_print[i,j] < 0:
                    eval("self.{}.config(background = 'red')".format(x))

                # divido altrimenti la variabile ha valore None
                eval("self.{}.grid(row=row_idx + i , column=col_idx + j)".format(x) )

        Label(frame, height=1, width=7, text="").grid(row=row_idx + i + 1, column=col_idx + j + 1)

        return row_idx + i + 1 , col_idx + j + 1



    def elaborate_request(self):
        print 'Eigenvalues positive:' + str(self.ck1.get())
        print 'Floor at 0 in Quarterly Transition Matrix:' + str(self.ck2.get())
        print 'Total amount each row sum equal to 100%:' + str(self.ck3.get())

        w, v = np.linalg.eig(self.matrix_input)
        # pongo gli autovalori maggiori di 0
        w_floor = np.clip(np.array(w), 0, None)
        # matrice diagonale trimestrale
        w_4 = np.power(w_floor, 0.25)

        diag = np.diag(w_4)
        a_trim = np.dot(v, np.dot(diag, np.linalg.inv(v)))

        #np.dot(a_trim, np.dot(a_trim, np.dot(a_trim, a_trim))) - mat

        if self.ck2.get() == 1:
            a_trim_floor = np.clip(a_trim, 0, None)
        else:
            a_trim_floor = a_trim

        base = np.array(np.sum(a_trim_floor, axis=1))

        if self.ck3.get() == 1:
            a_trim_new = a_trim_floor / base[:, None]
        else:
            a_trim_new = a_trim_floor

        # aggiorno la finestra
        self.update_matrix_output(a_trim_new)

        self.diff = np.abs(self.M_out_ann - self.matrix_input)

        # aggiorno le matrici da salvare in output
        self.update_matrix_show('QTM',self.M_out_trim)
        self.update_matrix_show('QTM_E',self.M_out_trim_diag)
        self.update_matrix_show('ATM',self.M_out_ann)
        self.update_matrix_show('ATM_E',self.M_out_ann_diag)



    def update_matrix_show(self,nv,matrix):
        for i in xrange(0,self.nRows):
            for j in xrange(0, self.nCols):
                y = nv + str(i) + str(j)
                eval("self.{}.set(str(round(100*matrix[i,j],2)))".format(y))

                # modifico lo sfondo di tutte le matrici perche' le altre sono disabilitate
                #  e non cambiano il loro comportamento

                x = 'E_' + nv + str(i) + str(j)
                if self.diff[i,j] > (self.ck4.get()/100.):
                    eval("self.{}.config(background = 'red')".format(x))
                else:
                    eval("self.{}.config(background = 'green')".format(x))


    def close_window(self, flag_save):
        if flag_save == True:
            self.flag_save = True
        else:
            self.flag_save = False

        self.master.destroy()


def writeQuarterlyMatrixResOnXls(xla, W_class, book, nameSheet, nameMatrix):

    allSheetInBook = allSheet(book)

    # -------------- controllo l'esistenza del foglio  ----------------
    if not (nameSheet in allSheetInBook):
        s = book.Sheets.Add()
        s.Name = nameSheet
    else:
        s = book.Sheets(nameSheet)
        s.Activate()
    # -----------------------------------------------------------------

    r = findCalibrationPos(xla, nameSheet)

    #---
    #mi posiziono nella prima cella utile per scrivere i risultati del fitting
    #---

    matrix_quarterly = pd.DataFrame(W_class.M_out_trim)
    matrix_yearly    = pd.DataFrame(W_class.M_out_ann)

    Attributi_q = \
        {     "0. MatrixType"              : 'Transition Matrix'
            , "1. Time Reference"          : 'Quarterly'
            , "2. Eigenvalues positive"    : str(W_class.ck1.get() == 1)
            , "3. Floor at 0 in Quarterly" : str(W_class.ck2.get() == 1)
            , "4. Amount row sum = 100%"   : str(W_class.ck3.get() == 1)
        }

    Attributi_a = \
        {     "0. MatrixType"             : 'Transition Matrix'
            , "1. Time Reference"         : 'Yearly Recalculated'
            , "2. Eigenvalues positive"   : str(W_class.ck1.get() == 1)
            , "3. Floor at 0 in Quarterly": str(W_class.ck2.get() == 1)
            , "4. Amount row sum = 100%"  : str(W_class.ck3.get() == 1)
        }

    row = r.Row
    col = r.Column
    r = intestazioneElabMatrix(xla = xla, rng = r, attributi = Attributi_q , title = nameMatrix + '_Quarterly')
    r = writeResultPandas(xla = xla , rng = r, df = matrix_quarterly, flagPrintColumns = False)

    r = xla.Range(xla.Cells(row,col + matrix_quarterly.shape[1] + 1),xla.Cells(row,col + matrix_quarterly.shape[1] + 1))
    r = intestazioneElabMatrix(xla = xla, rng = r, attributi = Attributi_a , title = nameMatrix + '_Yearly_Recalculated')
    r = writeResultPandas(xla = xla , rng = r, df = matrix_yearly, flagPrintColumns = False)
    xla.Cells.ColumnWidth = 20


def intestazioneElabMatrix( xla, rng,  attributi, nCols = 2, title= 'Calibration'):

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
        xla.Cells(topLeftRow + 1 + i, topLeftCol + 1).HorizontalAlignment = const.xlCenter
        i+=1

    rangeStart = xla.Range(xla.Cells(topLeftRow + nRows + 1, topLeftCol),xla.Cells(topLeftRow + nRows + 1, topLeftCol))

    return rangeStart


def writeTemplateQuarterlyMatrixInput(xla, nameSheet, dimMatrix):

    r = findCalibrationPos(xla, nameSheet)

    #---
    #mi posiziono nella prima cella utile per scrivere i risultati del fitting
    #---

    mat = pd.DataFrame(np.zeros((dimMatrix,dimMatrix)))

    Attributi = \
        {     "0. MatrixType"      : 'Transition Matrix'
        }

    row = r.Row
    col = r.Column
    r = intestazioneElabMatrix(xla = xla, rng = r, attributi = Attributi , title = 'Template Elaborate Input Matrix')
    r = writeResultPandas(xla = xla , rng = r, df = mat, flagPrintColumns = False)
    xla.Cells.ColumnWidth = 18
