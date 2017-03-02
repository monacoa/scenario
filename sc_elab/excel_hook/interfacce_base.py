from pyxll import xl_func, xl_app, xl_menu, xl_macro, xlcAlert
from xls_utils import *
from Tkinter import *
import tkMessageBox
from db_qrys import getCurvesListFromDb, getDatesListFromDb
from sc_elab.core.SwpCurve import *
from win32com.client import constants as const



#=======================================================================================================================
# punto di ingresso per load curve
# =======================================================================================================================

from W_swapCurve import W_curveType
from xls_swapCurve import writeCurveOnXls
from DEF_intef import nameSheetCurve

@xl_func
def load_swap_curve_from_db(control):
    nameSheet = nameSheetCurve
    xla = xl_app()
    book = xla.ActiveWorkbook
    #-----
    #creo foglio nameSheetCurve se non esiste
    try:
        s = book.Sheets(nameSheet)
        s.Activate()
    except:
        s = book.Sheets.Add()
        s.Name = nameSheet
   #------------------
    root = Tk()
    app  = W_curveType(root)
    root.mainloop()

    curve_des = app.new_window.new_window.curve
    curve_date= app.new_window.date

    cc = Curve()
    cc.ref_date = datetime.date(day=int(curve_date[-2:]), month=int(curve_date[5:7]), year=int(curve_date[:4]))
    cc.description= curve_des
    cc.loadDataFromDB()
    cc.init_finalize()

    writeCurveOnXls(cc, nameSheet, xla)


#=======================================================================================================================
# punto di ingresso per bootstrap
#=======================================================================================================================
from xls_bootCurve import writeBootstrapResOnXls
from xls_swapCurve import readCurveFromXls

from W_bootstrapCurve import W_bootstrapSelection
import ctypes
@xl_func
def bootstrap_from_xls(control):
    nameSheet = nameSheetCurve
    xla = xl_app()
    book = xla.ActiveWorkbook
    try:
        s = book.Sheets(nameSheet)
        s.Activate()
    except:
        root = Tk()
        msg = "Missing input sheet for Swap Curves in your workbook... \nNothing to do for me!"
        tkMessageBox.showinfo("Warning!", msg)
        root.destroy()
        return

    rangeStart = "B2"
    distance = 5

    curveL = readCurvesNames(xla,s,rangeStart,"o", distance)
    root = Tk()
    #root.wm_withdraw()
    W = W_bootstrapSelection(root, curveL)
    root.mainloop()

    curveDes = W.curve
    curvePos = W.pos

    #opt
    opt_swaps     = (str(W.new_window.variable1.get()).strip(""))[1]
    opt_fut_gap   = (str(W.new_window.variable2.get()).strip(""))[1]
    opt_conv_adj  = (str(W.new_window.variable3.get()).strip(""))[1]
    opt_out_cap   = (str(W.new_window.variable4.get()).strip(""))[1]
    opt_path_graph=  W.new_window.variable5.get()

    str_boot_opt = opt_swaps+","+opt_fut_gap + "," + opt_conv_adj + "," + opt_out_cap

    data_opt                    = {}

    data_opt['GapFutures']      = opt_fut_gap
    data_opt['SwapGapMethod']   = opt_swaps
    data_opt['RegimeOutput']    = opt_out_cap
    data_opt['Convexity']       = opt_conv_adj
    data_opt['MakeGraph']       = True
    data_opt['SaveGraph']       = True
    data_opt['FutureTenor']     = 90
    data_opt['Path']            = opt_path_graph

    curve    = readCurveFromXls(xla, curveDes, curvePos, nameSheet)
    boot_out = curve.bootstrap(data_opt)
    writeBootstrapResOnXls(curve, xla, str_boot_opt,boot_out)


#=======================================================================================================================
# punto di ingresso per fitting
#=======================================================================================================================

from W_fittingCurve import W_fittingType
from xls_bootCurve import readBootstrappedCurveFromXls
from xls_fittingCurve import writeFittingBootResOnXls, writeFittingPyResOnXls
from DEF_intef import nameSheetBootstrap
@xl_func
def fitting_from_xls(control):
    nameSheet = nameSheetBootstrap
    xla = xl_app()
    book = xla.ActiveWorkbook
    try:
        s = book.Sheets(nameSheet)
        s.Activate()
    except:
        root = Tk()
        msg = "Missing sheet  for Bootstrap Elab in your workbook... \nNothing to do for me!"
        tkMessageBox.showinfo("Warning!", msg)
        root.destroy()
        return

    root = Tk()
    # root.wm_withdraw()
    #W = W_fittingSelection(root, curveL)
    W = W_fittingType(root)
    root.mainloop()

    curveDes = W.new_window.curve
    curvePos = W.new_window.pos
    fit_type = W.fit_type

    print "des:", curveDes
    print "pos:", curvePos,  type(curvePos)
    print "TYPE", fit_type

    opt_dict = {}
    opt_dict['interp']          = (str(W.new_window.new_window.variable1.get()).strip(""))[1]
    opt_dict['opt_fwd_tenor']   = (str(W.new_window.new_window.variable2.get()).strip(""))[1]
    opt_dict['opt_path_graph']  =  W.new_window.new_window.variable5.get()
    opt_dict['fit_type']        = fit_type
    print "**********************options:", opt_dict

    if fit_type == "boot":
        Bcurve = readBootstrappedCurveFromXls(xla, curveDes, curvePos, nameSheet)
        res = Bcurve.fittingFromBoot(opt_dict)
        writeFittingBootResOnXls(Bcurve, xla, opt_dict, res, curvePos)
    else:
        nameSheet = nameSheetCurve
        Bcurve = readCurveFromXls(xla, curveDes, curvePos, nameSheet)
        Bcurve.show()
        res = Bcurve.fittingFromPY(opt_dict)
        writeFittingPyResOnXls(Bcurve, xla, opt_dict, res, curvePos)