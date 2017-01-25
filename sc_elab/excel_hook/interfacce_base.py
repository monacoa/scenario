

from pyxll import xl_func, xl_app, xl_menu, xl_macro, xlcAlert
from msilib import Control

#@xl_func
#def hello(name):
#    return "Hello, %s" % name



@xl_menu("win32com test", sub_menu="More Examples")
def win32com_menu_test(control):
    # get the current selected range and set some text
    selection = xl_app().Selection
    selection.Value = "Hello!"
    xlcAlert("Some text has been written to the current cell")


@xl_func
def button_example(control):
    xl = xl_app()
    range = xl.Range("PLUTO")
    range.Value = range.Value + 1


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



#def on_text_button(a = 1):
    
#    return 'pippo'
    #xl = xl_app()
    
    #out = sommaro()
    
    #xl.Selection.Value = "%s" %(out)
