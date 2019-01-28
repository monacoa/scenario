from win32com.client import constants as const
from Tkinter import *
import tkMessageBox


class W_template (Frame):
    def __init__(self, master = None):
        Frame.__init__(self, master)
        self.master = master

        tmp_type = StringVar()
        tmp_type.set('CCDS')

        Label(self,text="""Choose template :""",justify=LEFT,padx=100).pack()

        self.rb1 = Radiobutton(self,text='Bond',justify='left',variable=tmp_type,value='BOND').pack(anchor=W)
        self.rb2 = Radiobutton(self,text='Curve, CDS, CapFloorSwaptions, Opzioni',justify='left',variable=tmp_type,value='CCDS').pack(anchor=W)
        self.rb3 = Radiobutton(self,text='TS'   ,justify='left',variable=tmp_type,value='TS',state="disabled").pack(anchor=W)
        self.rb4 = Radiobutton(self,text='Emittente',justify='left',variable=tmp_type,value='EMT',state="disabled").pack(anchor=W)

        self.template = False

        # create button
        self.btn2 = Button(self, text="Cancel",  command=lambda:self.close_window())
        self.btn2.pack(side=BOTTOM, fill = 'x')
        # create button
        self.btn1 = Button(self, text="Select", command= lambda:self.download_template(tmp_type.get()))
        self.btn1.pack(side=BOTTOM, fill='x')
        self.pack()

    def download_template(self, tmp_type):

        if tmp_type == 'CCDS':
            self.template = ['Dati','Curve','CDS','CapFloorSwaptions','Opzioni']
            self.close_window()

        elif tmp_type == 'BOND':
            self.template = ['Bond_master']
            self.close_window()

        elif tmp_type == 'TS':
            self.template = False
            self.close_window()

        elif tmp_type == 'EMT':
            self.close_window()

        else:
            self.donothing()

    def close_window(self):
        self.master.destroy()

    def donothing(self):
        tkMessageBox.showinfo("Nothing To do", "bye bye")
        self.master.destroy()


def writeTemplate(xla,wb,nameSheet,ttt):
    s = wb.Sheets.Add()
    s.Name = nameSheet

    i = 1
    for k in ttt.keys():
        xla.Cells(1,i).Value = k
        i = i + 1

    formatTemplate(xla = xla, nRow = 15, nCol = i-1)


def formatTemplate(xla,nRow,nCol):

    #xla.Cells.Columns.ColumnWidth = 20
    s = xla.Cells.Columns.AutoFit()

    xla.Range(xla.Cells(1,1),xla.Cells(1,nCol)).Select()
    xla.Selection.HorizontalAlignment   = const.xlCenter
    xla.Selection.VerticalAlignment     = const.xlBottom
    xla.Selection.WrapText              = False
    xla.Selection.Orientation           = 0
    xla.Selection.AddIndent             = False
    xla.Selection.IndentLevel           = 0
    xla.Selection.ShrinkToFit           = False
    xla.Selection.ReadingOrder          = const.xlContext
    xla.Selection.MergeCells            = False
    xla.Selection.Font.ColorIndex       = 2
    xla.Selection.Font.Bold             = True
    xla.Selection.Interior.ColorIndex   = 50
    xla.Selection.Interior.Pattern      = const.xlSolid

    for i in xrange(2,nRow + 1):

        xla.Range(xla.Cells(i, 1), xla.Cells(i, nCol)).Select()

        if (i % 2) == 0:
            xla.Selection.Interior.ColorIndex = 2
            xla.Selection.Interior.Pattern = const.xlSolid
        else:
            xla.Selection.Interior.ColorIndex = 35
            #xla.Selection.Interior.Pattern = const.xlSolid



def ask_question(header,msg):
    root = Tk()
    root.withdraw()
    answer = tkMessageBox.askquestion(header, msg, icon='warning')
    root.destroy()

    return answer
