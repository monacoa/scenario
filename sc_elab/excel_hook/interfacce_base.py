

from pyxll import xl_func, xl_app, xl_menu, xl_macro, xlcAlert
from msilib import Control




import Tkinter
import tkMessageBox



@xl_menu("AAAA")
def test_tk():
    top = Tkinter.Tk()
    
    def helloCallBack():
       tkMessageBox.showinfo( "Hello Python", "Hello World")
    
    B = Tkinter.Button(top, text ="Hello", command = helloCallBack)
    
    B.pack()
    top.mainloop()


#@xl_func
#def hello(name):
#    return "Hello, %s" % name



@xl_menu("win32com test", sub_menu="More Examples")
def win32com_menu_test(control):
    # get the current selected range and set some text
    selection = xl_app().Selection
    selection.Value = "Hello!"
    xlcAlert("Some text has been written to the current cell")



@xl_macro
def checkbox_example():
    xl = xl_app()
    check_box = xl.ActiveSheet.CheckBoxes(xl.Caller)
    if check_box.Value == 1:
        xl.Range("checkbox_output").Value = "CHECKED"
    else:
        xl.Range("checkbox_output").Value = "Click the check box"


"""
#@xl_func
def checkbox_example(control):
    xl = xl_app()
    check_box = xl.ActiveSheet.CheckBoxes(xl.Caller)
    if check_box.Value == 1:
        xl.Range("checkbox_output").Value = "CHECKED"
    else:
        xl.Range("checkbox_output").Value = "Click the check box"
"""


@xl_func
def button_example(control):

    xl          = xl_app()
    range       = xl.Range("PLUTO")
    range_p     = xl.Range("PIPPO")


    range.Value = range.Value + 1
    
    #range.Interior.ColorIndex = 6
    #range_p.Interior.ColorIndex = 6
    range_p.Style = "Normal"
    
    range_p.Style.Font.Name = "Verdana"
    range_p.Style.Font.Size = 16
    range_p.Style.Font.Color = System.Drawing.ColorTranslator.ToOle(System.Drawing.Color.Red);

    #range_p.Style.Font.Color = System.Drawing.ColorTranslator.ToOle(System.Drawing.Color.Red)
    #range_p.Style.Interior.Color = System.Drawing.ColorTranslator.ToOle(System.Drawing.Color.Gray)
    #range_p.Style.Interior.Pattern = Excel.XlPattern.xlPatternSolid

    
    #ws = wb.Worksheets("Sheet1")
    """
    for i in range (1,21):
        xl.Cells(i,1).Value = i
        #xl.Cells(i,1).Interior.ColorIndex = i
    """


@xl_func
def on_text_button(control):
    xl = xl_app()
    xl.Selection.Value = "This text was added by Comollis."


@xl_menu("on_text_button2")
def on_text_button2(control):
    xl = xl_app()    
    #out = sommaro()   
    #xl.Selection.Value = "%s" %(out)
    xl.Selection.Value = "AAA" 
    
    #return a + 2*b

@xl_func
def sommaro(a= 1, b = 2):
    
    return a + 2*b

    """
    import win32com.client as win32
    excel = win32.gencache.EnsureDispatch('Excel.Application')
    wb = excel.Workbooks.Add()
    ws = wb.Worksheets("Sheet1")
    for i in range (1,21):
        ws.Cells(i,1).Value = i
        ws.Cells(i,1).Interior.ColorIndex = i
    wb.SaveAs('cell_color.xlsx')
    excel.Application.Quit()
    """


@xl_func
def sommario(control):
    xl = xl_app()
    xl.Selection.Value = "Son Mario"

# Connection example: Windows, without a DSN, using the Windows SQL Server driver
# Connection example: Windows, without a DSN, using the Windows SQL Server driver
#def on_text_button(a = 1):
    
#    return 'pippo'
    #xl = xl_app()
    
    #out = sommaro()
    
    #xl.Selection.Value = "%s" %(out)
