from pyxll import xl_func, xl_app, xl_menu, xl_macro, xlcAlert
from xls_utils import *
from Tkinter import *
import tkMessageBox
from db_qrys import getCurvesListFromDb, getDatesListFromDb
from sc_elab.core.SwpCurve import *
from sc_elab.core.BondPortfolio import *

from win32com.client import constants as const



#=======================================================================================================================
# punto di ingresso per load curve
# =======================================================================================================================

from W_swapCurve import W_curveType,W_curveDate
from W_bondPortfolio import W_bondType,W_bondDate

from xls_swapCurve import writeCurveOnXls
from xls_bondPortfolio import writePortfoliOnXls

from DEF_intef import nameSheetCurve, nameSheetCDS, nameSheetBond

@xl_func
def load_swap_curve_from_db(control):
    root = Tk()
    app  = W_curveType(root)
    root.mainloop()

    curve_des = app.new_window.new_window.curve
    curve_date= app.new_window.date
    curve_type = app.new_window.type

    print "descrizione:", curve_des
    print "data:", curve_date
    print "type:", curve_type


    xla = xl_app()
    book = xla.ActiveWorkbook
    # -----
    # creo foglio nameSheetCurve se non esiste
    if curve_type == "SWP":
        nameSheet = nameSheetCurve
        try:
            s = book.Sheets(nameSheet)
            s.Activate()
        except:
            s = book.Sheets.Add()
            s.Name = nameSheet
            # ------------------
        cc = Curve()
    elif curve_type == "CDS":
        nameSheet = nameSheetCDS
        try:
            s = book.Sheets(nameSheet)
            s.Activate()
        except:
            s = book.Sheets.Add()
            s.Name = nameSheet
            # ------------------
        cc = CdsCurve()

    cc.ref_date = datetime.date(day=int(curve_date[-2:]), month=int(curve_date[5:7]), year=int(curve_date[:4]))
    cc.description= curve_des

    if curve_type == "CDS":
        cc.ratingProvider = app.new_window.new_window.new_window.rating.get()
        cc.sectorProvider = app.new_window.new_window.new_window.sector.get()

    cc.loadDataFromDB()
    if curve_type == "SWP": cc.init_finalize()

    writeCurveOnXls(cc, nameSheet, xla, curve_type)

#=======================================================================================================================
# punto di ingresso per load CDS
# =======================================================================================================================

@xl_func
def load_cds_curve_from_db(control):
    root = Tk()
    curve_type = "CDS"

    app =  W_curveDate(master = root, parent = None, type = curve_type)
    root.mainloop()

    curve_des = app.new_window.curve
    curve_date= app.date

    xla  = xl_app()
    book = xla.ActiveWorkbook
    # -----
    # creo foglio nameSheetCurve se non esiste
    # -----
    nameSheet = nameSheetCDS
    try:
        s = book.Sheets(nameSheet)
        s.Activate()
    except:
        s = book.Sheets.Add()
        s.Name = nameSheet
        # ------------------
    cc = CdsCurve()

    cc.ref_date       = datetime.date(day=int(curve_date[-2:]), month=int(curve_date[5:7]), year=int(curve_date[:4]))
    cc.description    = curve_des

    cc.ratingProvider = app.new_window.new_window.rating.get()
    cc.sectorProvider = app.new_window.new_window.sector.get()
    cc.loadDataFromDB()
    writeCurveOnXls(cc, nameSheet, xla, curve_type)

#=======================================================================================================================
# punto di ingresso per load BondData
# =======================================================================================================================

@xl_func
def load_bond_data_from_db(control):
    
    root = Tk()
    data_type = "BOND"
    
    

    app =  W_bondDate(master = root, parent = None, type = data_type)
    root.mainloop()

    bond_des = app.new_window.bond
    bond_date= app.date

    xla  = xl_app()
    book = xla.ActiveWorkbook
    # -----
    # creo foglio nameSheetCurve se non esiste
    # -----
    nameSheet = nameSheetBond
    try:
        s = book.Sheets(nameSheet)
        s.Activate()
    except:
        s = book.Sheets.Add()
        s.Name = nameSheet
        # ------------------
    cc = BondPortfolio()

    cc.ref_date       = datetime.date(day=int(bond_date[-2:]), month=int(bond_date[5:7]), year=int(bond_date[:4]))
    cc.description    = bond_des
    

    #cc.ratingProvider = app.new_window.new_window.rating.get()
    #cc.sectorProvider = app.new_window.new_window.sector.get()
    cc.loadDataFromDB(bond_date)
    writePortfoliOnXls(cc, nameSheet, xla, data_type)


#=======================================================================================================================
# punto di ingresso per bootstrap curve SWAP
#=======================================================================================================================

from xls_bootCurve import writeBootstrapResOnXls
from xls_bootCurve import writeCDSBootstrapRes1OnXls, writeCDSBootstrapRes2OnXls
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
    W = W_bootstrapSelection(root, curveL, "SWP")
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

    curve        = readCurveFromXls(xla, curveDes, curvePos, nameSheet)
    codeL, codeR = curve.getCurveCode()
    
    boot_out     = curve.bootstrap(data_opt)
    #print "risultati bootstrap:", boot_out
    print 'YYYYYYYYYYYYYYYYYY'
    

    if boot_out == None:
        # significa che ho intercettato un errore!
        root = Tk()
        root.withdraw()
        msg = "Unable to perform curve bootstrap!"
        tkMessageBox.showinfo("ERROR!", msg)
        root.destroy()
        return

    writeBootstrapResOnXls(curve, xla, str_boot_opt,boot_out, codeL, codeR)


#=======================================================================================================================
# punto di ingresso per fitting
#=======================================================================================================================

from W_fittingCurve import W_fittingType
from xls_bootCurve import readBootstrappedCurveFromXls
from xls_fittingCurve import writeFittingBootResOnXls, writeFittingPyResOnXls
from DEF_intef import nameSheetBootstrap

@xl_func
def fitting_from_xls(control):

    root = Tk()
    # root.wm_withdraw()
    W = W_fittingType(root)
    root.mainloop()

    curveDes = W.new_window.curve
    curvePos = W.new_window.pos
    fit_type = W.fit_type


    opt_dict = {}
    opt_dict['interp']          = (str(W.new_window.new_window.variable1.get()).strip(""))[1]
    opt_dict['opt_fwd_tenor']   = (str(W.new_window.new_window.variable2.get()).strip(""))[1]
    opt_dict['opt_path_graph']  =  W.new_window.new_window.variable5.get()
    opt_dict['fit_type']        = fit_type

    opt_dict['bound_min_sve'] = []
    opt_dict['bound_max_sve'] = []
    opt_dict['sve_params']    = ['const1', 'const2', 'beta0', 'beta1', 'beta2', 'beta3']
    opt_dict['bound_min_sve'].append( float(W.new_window.new_window.SVE_c1min.get().strip("")))
    opt_dict['bound_min_sve'].append( float(W.new_window.new_window.SVE_c2min.get().strip("")))
    opt_dict['bound_min_sve'].append( float(W.new_window.new_window.SVE_b0min.get().strip("")))
    opt_dict['bound_min_sve'].append( float(W.new_window.new_window.SVE_b1min.get().strip("")))
    opt_dict['bound_min_sve'].append( float(W.new_window.new_window.SVE_b2min.get().strip("")))
    opt_dict['bound_min_sve'].append( float(W.new_window.new_window.SVE_b3min.get().strip("")))

    opt_dict['bound_max_sve'].append( float(W.new_window.new_window.SVE_c1max.get().strip("")))
    opt_dict['bound_max_sve'].append( float(W.new_window.new_window.SVE_c2max.get().strip("")))
    opt_dict['bound_max_sve'].append( float(W.new_window.new_window.SVE_b0max.get().strip("")))
    opt_dict['bound_max_sve'].append( float(W.new_window.new_window.SVE_b1max.get().strip("")))
    opt_dict['bound_max_sve'].append( float(W.new_window.new_window.SVE_b2max.get().strip("")))
    opt_dict['bound_max_sve'].append( float(W.new_window.new_window.SVE_b3max.get().strip("")))

    opt_dict['bound_min_ns'] = []
    opt_dict['bound_max_ns'] = []
    opt_dict['ns_params']    = ['const1', 'beta0', 'beta1', 'beta2']
    
    opt_dict['bound_min_ns'].append( float(W.new_window.new_window.NS_c1min.get().strip("")))
    opt_dict['bound_min_ns'].append( float(W.new_window.new_window.NS_b0min.get().strip("")))
    opt_dict['bound_min_ns'].append( float(W.new_window.new_window.NS_b1min.get().strip("")))
    opt_dict['bound_min_ns'].append( float(W.new_window.new_window.NS_b2min.get().strip("")))

    opt_dict['bound_max_ns'].append( float(W.new_window.new_window.NS_c1max.get().strip("")))
    opt_dict['bound_max_ns'].append( float(W.new_window.new_window.NS_b0max.get().strip("")))
    opt_dict['bound_max_ns'].append( float(W.new_window.new_window.NS_b1max.get().strip("")))
    opt_dict['bound_max_ns'].append( float(W.new_window.new_window.NS_b2max.get().strip("")))

    
    opt_dict['bound_min_cir'] = []
    opt_dict['bound_max_cir'] = []
    opt_dict['cir_params'] =['r0', 'kappa', 'theta', 'sigma']
    opt_dict['bound_min_cir'].append( float(W.new_window.new_window.CIR_r0min.get().strip("")))
    opt_dict['bound_min_cir'].append( float(W.new_window.new_window.CIR_kmin.get().strip("")))
    opt_dict['bound_min_cir'].append( float(W.new_window.new_window.CIR_thetamin.get().strip("")))
    opt_dict['bound_min_cir'].append( float(W.new_window.new_window.CIR_sigmamin.get().strip("")))

    opt_dict['bound_max_cir'].append( float(W.new_window.new_window.CIR_r0max.get().strip("")))
    opt_dict['bound_max_cir'].append( float(W.new_window.new_window.CIR_kmax.get().strip("")))
    opt_dict['bound_max_cir'].append( float(W.new_window.new_window.CIR_thetamax.get().strip("")))
    opt_dict['bound_max_cir'].append( float(W.new_window.new_window.CIR_sigmamax.get().strip("")))

    #---
    xla = xl_app()
    book = xla.ActiveWorkbook

    if fit_type == "boot":
        nameSheet   = nameSheetBootstrap
        Bcurve      = readBootstrappedCurveFromXls(xla, curveDes, curvePos, nameSheet)
        res         = Bcurve.fittingFromBoot(opt_dict)
        writeFittingBootResOnXls(Bcurve, xla, opt_dict, res, curvePos)
    else:
        nameSheet   = nameSheetCurve
        Bcurve      = readCurveFromXls(xla, curveDes, curvePos, nameSheet)
        Bcurve.show()
        res         = Bcurve.fittingFromPY(opt_dict)
        writeFittingPyResOnXls(Bcurve, xla, opt_dict, res, curvePos)

#=======================================================================================================================
# punto di ingresso per bond fitting
#=======================================================================================================================

from W_fittingCurve import W_fittingType
from xls_bootCurve import readBootstrappedCurveFromXls
from xls_fittingCurve import writeFittingBootResOnXls, writeFittingPyResOnXls
from DEF_intef import nameSheetBootstrap

@xl_func
def bond_fitting_from_xls(control):

    root = Tk()
    W = W_fittingType(root)
    root.mainloop()

    curveDes = W.new_window.curve
    curvePos = W.new_window.pos
    fit_type = W.fit_type


    opt_dict = {}
    opt_dict['interp']          = (str(W.new_window.new_window.variable1.get()).strip(""))[1]
    opt_dict['opt_fwd_tenor']   = (str(W.new_window.new_window.variable2.get()).strip(""))[1]
    opt_dict['opt_path_graph']  =  W.new_window.new_window.variable5.get()
    opt_dict['fit_type']        = fit_type

    opt_dict['bound_min_sve'] = []
    opt_dict['bound_max_sve'] = []
    opt_dict['sve_params']    = ['const1', 'const2', 'beta0', 'beta1', 'beta2', 'beta3']
    opt_dict['bound_min_sve'].append( float(W.new_window.new_window.SVE_c1min.get().strip("")))
    opt_dict['bound_min_sve'].append( float(W.new_window.new_window.SVE_c2min.get().strip("")))
    opt_dict['bound_min_sve'].append( float(W.new_window.new_window.SVE_b0min.get().strip("")))
    opt_dict['bound_min_sve'].append( float(W.new_window.new_window.SVE_b1min.get().strip("")))
    opt_dict['bound_min_sve'].append( float(W.new_window.new_window.SVE_b2min.get().strip("")))
    opt_dict['bound_min_sve'].append( float(W.new_window.new_window.SVE_b3min.get().strip("")))

    opt_dict['bound_max_sve'].append( float(W.new_window.new_window.SVE_c1max.get().strip("")))
    opt_dict['bound_max_sve'].append( float(W.new_window.new_window.SVE_c2max.get().strip("")))
    opt_dict['bound_max_sve'].append( float(W.new_window.new_window.SVE_b0max.get().strip("")))
    opt_dict['bound_max_sve'].append( float(W.new_window.new_window.SVE_b1max.get().strip("")))
    opt_dict['bound_max_sve'].append( float(W.new_window.new_window.SVE_b2max.get().strip("")))
    opt_dict['bound_max_sve'].append( float(W.new_window.new_window.SVE_b3max.get().strip("")))

    opt_dict['bound_min_ns'] = []
    opt_dict['bound_max_ns'] = []
    opt_dict['ns_params']    = ['const1', 'beta0', 'beta1', 'beta2']
    
    opt_dict['bound_min_ns'].append( float(W.new_window.new_window.NS_c1min.get().strip("")))
    opt_dict['bound_min_ns'].append( float(W.new_window.new_window.NS_b0min.get().strip("")))
    opt_dict['bound_min_ns'].append( float(W.new_window.new_window.NS_b1min.get().strip("")))
    opt_dict['bound_min_ns'].append( float(W.new_window.new_window.NS_b2min.get().strip("")))

    opt_dict['bound_max_ns'].append( float(W.new_window.new_window.NS_c1max.get().strip("")))
    opt_dict['bound_max_ns'].append( float(W.new_window.new_window.NS_b0max.get().strip("")))
    opt_dict['bound_max_ns'].append( float(W.new_window.new_window.NS_b1max.get().strip("")))
    opt_dict['bound_max_ns'].append( float(W.new_window.new_window.NS_b2max.get().strip("")))

    
    opt_dict['bound_min_cir'] = []
    opt_dict['bound_max_cir'] = []
    opt_dict['cir_params'] =['r0', 'kappa', 'theta', 'sigma']
    opt_dict['bound_min_cir'].append( float(W.new_window.new_window.CIR_r0min.get().strip("")))
    opt_dict['bound_min_cir'].append( float(W.new_window.new_window.CIR_kmin.get().strip("")))
    opt_dict['bound_min_cir'].append( float(W.new_window.new_window.CIR_thetamin.get().strip("")))
    opt_dict['bound_min_cir'].append( float(W.new_window.new_window.CIR_sigmamin.get().strip("")))

    opt_dict['bound_max_cir'].append( float(W.new_window.new_window.CIR_r0max.get().strip("")))
    opt_dict['bound_max_cir'].append( float(W.new_window.new_window.CIR_kmax.get().strip("")))
    opt_dict['bound_max_cir'].append( float(W.new_window.new_window.CIR_thetamax.get().strip("")))
    opt_dict['bound_max_cir'].append( float(W.new_window.new_window.CIR_sigmamax.get().strip("")))

    #---
    xla = xl_app()
    book = xla.ActiveWorkbook

    if fit_type == "boot":
        nameSheet   = nameSheetBootstrap
        Bcurve      = readBootstrappedCurveFromXls(xla, curveDes, curvePos, nameSheet)
        res         = Bcurve.fittingFromBoot(opt_dict)
        writeFittingBootResOnXls(Bcurve, xla, opt_dict, res, curvePos)
    else:
        nameSheet   = nameSheetCurve
        Bcurve      = readCurveFromXls(xla, curveDes, curvePos, nameSheet)
        Bcurve.show()
        res         = Bcurve.fittingFromPY(opt_dict)
        writeFittingPyResOnXls(Bcurve, xla, opt_dict, res, curvePos)

#=======================================================================================================================
# punto di ingresso per SAVE
#=======================================================================================================================

from W_saveCurve        import  W_saveType
from db_saveCurve       import  saveInterpolationParmsOnDb
from db_saveCurve       import  saveZcDfOnDB, deletingDbCurves
from xls_fittingCurve   import  readInterpolatedCurveFromXls

@xl_func
def save_from_xls(control):
    root = Tk()
    #root.wm_withdraw()
    W = W_saveType(root)
    root.mainloop()

    type = W.saveType

    if type == "Boot":
        des         = W.new_window.curve
        pos         = W.new_window.pos
        DF          = W.new_window.new_window.var1.get()
        ZC          = W.new_window.new_window.var2.get()
        xla         = xl_app()
        nameSheet   = nameSheetBootstrap
        res,codes   = saveZcDfOnDB(xla, nameSheet, des, pos, DF, ZC)
        #---
        root2 = Tk()
        root2.withdraw()
        if (res):
            msg = "Bootstrap results are on DB, well done Comollis!"
            tkMessageBox.showinfo("YES WE CAN!", msg)
        else:
            ans = tkMessageBox.askquestion("Unable to save Bootstrap results because they're already on DB.", "DELETING... Are You Sure?", icon='warning')
            if ans =='yes':
                print "ho risposto yes, entro in delete"
                deletingDbCurves(codes)
                print "ora entro in save"
                r, cd = saveZcDfOnDB(xla, nameSheet, des, pos, DF, ZC)
                if not r:
                    msg = "Something's wrong..... SEPPUKU!!!!!"
                    tkMessageBox.showinfo("x@!#!", msg)

                else:
                    msg = "Bootstrap results are on DB, well done Comollis!"
                    tkMessageBox.showinfo("YES WE CAN!", msg)
            else:
                msg = "Unable to save Bootstrap results because they're already on DB... Please delete IT before!!"
                tkMessageBox.showinfo("x@!#!", msg)
        root2.destroy()
    elif type == 'FitFromBoot':
        pos_curve   =  W.new_window.pos_curve
        des_curve   =  W.new_window.curve
        pos_parms   =  W.new_window.new_window.pos_parms
        des_parms   =  W.new_window.new_window.parms
        xla         =  xl_app()
        nameSheet   =  nameSheetBootstrap
        int_curve   =  readInterpolatedCurveFromXls(xla,nameSheet, pos_curve, pos_parms)
        # ---
        saveInterpolationParmsOnDb(int_curve)
        # ---
    else:
        msg = "Unable to save your selection...\n but well done Comollis, anywhere!! \n;)"
        tkMessageBox.showinfo("OOOOPS!", msg)

# ==========================================
# punto d'ingresso per bootstrap CDS
# ==========================================

@xl_func
def bootstrap_cds_from_xls(control):
    nameSheet = nameSheetCDS
    xla       = xl_app()
    book      = xla.ActiveWorkbook
    try:
        s = book.Sheets(nameSheet)
        s.Activate()
    except:
        root = Tk()
        msg = "Missing input sheet for CDS Curves in your workbook... \nNothing to do for me!"
        tkMessageBox.showinfo("Warning!", msg)
        root.destroy()
        return

    rangeStart = "B2"
    distance = 5

    curveL = readCurvesNames(xla,s,rangeStart,"o", distance)

    root = Tk()
    W = W_bootstrapSelection(root, curveL = curveL, type = "CDS")
    root.mainloop()

    curveDes = W.curve
    curvePos = W.pos
    

    opt_boot_meth   = (str(W.new_window.variable1.get()).strip(""))[1]
    opt_rf_interp   = (str(W.new_window.variable6.get()).strip(""))[1]
    opt_hr_interp   = (str(W.new_window.variable7.get()).strip(""))[1]

    str_boot_opt = opt_boot_meth+","+opt_rf_interp + "," + opt_hr_interp
    data_opt                    = {}
    
    data_opt['hr_bootMethod']  = opt_boot_meth #0 = LCS, 1 = CHR
    data_opt['bench_interp']   = opt_rf_interp
    data_opt['hr_interp']      = opt_hr_interp
    
    curve_xl        = readCurveFromXls(xla, curveDes, curvePos, nameSheet, "CDS")

    opt_download = {}
    
    interp_rf_model = curve_xl.mapCodeModelInv(opt_rf_interp)
    


    opt_download['refDate'] = curve_xl.ref_date
    opt_download['valuta'] = curve_xl.curr
    opt_download['rating'] = curve_xl.rating
    opt_download['seniority'] = curve_xl.seniority

    opt_download['tipo_modello'] = interp_rf_model


    codeBenchList = ['%LS', '%DS', '%LFS', '%DFS']

    for codeTmp in codeBenchList:

        opt_download['codeSeg'] = codeTmp
        flag_loaded = curve_xl.loadBenchDataFromDB(opt_download)
        if (flag_loaded == 1):
            break
    
    
    
    if flag_loaded == 0:
        # significa che ho intercettato un errore!
        root = Tk()
        root.withdraw()
        msg0 = "Curva benchmark associata al modello %s non presente alla data del %s, cambia modello o data!!" %(interp_rf_model, curve_xl.ref_date)
        tkMessageBox.showinfo("Attenzione!!", msg0)

        root.destroy()
        return

    

    data_opt['DataRef']        = curve_xl.ref_date
    data_opt['RecoveryRate'] = float(curve_xl.recovery)/100.0
    data_opt['opt_path_graph']  =  'C:\\'
    
    yy = str(curve_xl.ref_date.year) 
    gg = str(curve_xl.ref_date.day)
    mm = str(curve_xl.ref_date.month)
    yy = yy[1:]


    curve_xl.ref_date

    opt_boot_meth_s = curve_xl.mapBootCDS(int(opt_boot_meth))
    
    codice_curva = 'CSN3' + curve_xl.curr + 'HR' +'BLM' + '0' + yy + mm + gg + 'CFRIL' + 'SW' + opt_boot_meth_s

    curve_xl.cds_boot_method    = data_opt['hr_bootMethod']
    curve_xl.rf_interp_type     = data_opt['bench_interp']  
    curve_xl.recovery           = data_opt['RecoveryRate']


    # -------------- elaborazione ---------------------------
    boot_out     = curve_xl.bootstrap(data_opt)
    print 'boot_out: ', boot_out
    
    #boot_out['CodiceCurva'] = 'XXX'

    
    if boot_out == None:
        # significa che ho intercettato un errore!
        root = Tk()
        root.withdraw()
        msg = "Unable to perform CDS bootstrap!"
        tkMessageBox.showinfo("ERROR!", msg)
        root.destroy()
        return
    
    writeCDSBootstrapRes1OnXls(curve_xl, xla, str_boot_opt, boot_out, codice_curva)
    writeCDSBootstrapRes2OnXls(curve_xl, xla, str_boot_opt, boot_out, codice_curva)
    
   
