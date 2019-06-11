from pyxll import xl_func, xl_app
from xls_utils import *
from sc_elab.core.SwpCurve import *
from sc_elab.core.BondPortfolio import *

from sc_elab.core import funzioni_bond_fitting as bf




#=======================================================================================================================
# punto di ingresso per load curve
# =======================================================================================================================

from W_swapCurve import W_curveType,W_curveDate
from W_bondPortfolio import W_bondDate, W_bondFitting

from xls_swapCurve import writeCurveOnXls
from xls_bondPortfolio import writePortfoliOnXls,readPortfolioFromXls

from DEF_intef import nameSheetCurve, nameSheetCDS, nameSheetBond


@xl_func
def load_swap_curve_from_db(control):
    root = Tk()
    app  = W_curveType(root)
    root.mainloop()

    curve_des  = app.new_window.new_window.curve
    curve_date = app.new_window.date
    curve_type = app.new_window.type



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

    print 'curve_des:',curve_des

    xla  = xl_app()
    book = xla.ActiveWorkbook

    # -----
    # creo foglio nameSheetCurve se non esiste
    # -----
    nameSheet = nameSheetCDS

    allSheetInBook = allSheet(book)
    # -------------- controllo l'esistenza del foglio  ----------------
    if not (nameSheetElabMatrix in allSheetInBook):
        s = book.Sheets.Add()
        s.Name = nameSheet
    else:
        s = book.Sheets(nameSheet)
        s.Activate()
    # -----------------------------------------------------------------

    for x in xrange(0,len(curve_des)):
        cc = CdsCurve()

        cc.ref_date       = datetime.date(day=int(curve_date[-2:]), month=int(curve_date[5:7]), year=int(curve_date[:4]))
        cc.description    = curve_des[x]

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

from xls_bootCurve import writeBootstrapResOnXls,writeCDSBootstrapRes1OnXls, writeCDSBootstrapRes2OnXls 
from xls_bootCurve import writeCDSBootstrapRes1OnXls_surv_m, writeCDSBootstrapRes1OnXls_md_m, writeCDSBootstrapRes1OnXls_hr_m 
from xls_bootCurve import writeCDSBootstrapRes2OnXls_zc_m, writeCDSBootstrapRes2OnXls_py_m

from xls_swapCurve import readCurveFromXls

from xls_bondPortfolio import writeBondFittingRes1OnXls, writeBondFittingRes2OnXls, writeBondFittingRes3OnXls, writeBondFittingRes4OnXls

from W_bootstrapCurve import W_bootstrapSelection


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
    curvePos = int(W.pos[0])


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

    #curve        = readPortfolioFromXls(xla, curveDes, curvePos, nameSheet)
    curve        = readCurveFromXls(xla, curveDes, curvePos, nameSheet, type = "SWP")
    codeL, codeR = curve.getCurveCode()
    
    boot_out     = curve.bootstrap(data_opt)
    

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
    
    opt_dict['MakeGraph'] = True

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


    nameSheet = nameSheetBond
    xla       = xl_app()
    book      = xla.ActiveWorkbook
    try:
        s = book.Sheets(nameSheet)
        s.Activate()
    except:
        root = Tk()
        msg = "Non e presente il foglio dei dati relativi ai Bond!!"
        tkMessageBox.showinfo("Warning!", msg)
        root.destroy()
        return

    rangeStart = "B2"
    distance = 2
    

    curveL = readCurvesNames(xla,s,rangeStart,"v", distance, 0)
    root = Tk()
    #W = W_bootstrapSelection(root, curveL = curveL, type = "BOND") # da modificare 
    W = W_bondFitting(root, curveL = curveL, type = "BOND") # da modificare
    
    #W2 = W_setParameters(root)

    root.mainloop()
    

    #par_bnd_min = W.new_window.new_window2.bound_min
    par_bnd = W.new_window.window_prms.x_bnd
    par_x0      = W.new_window.window_prms.x0

    curveDes = W.curve
    curvePos = W.pos

    hr_model_tmp   = (str(W.new_window.variable1.get()).strip(""))[1]
    rf_interp_tmp   = (str(W.new_window.variable6.get()).strip(""))[1]
    bond_model_tmp   = (str(W.new_window.variable7.get()).strip(""))[1]

    hr_model_map = {}
    hr_model_map[0] = 'SVE'
    hr_model_map[1] = 'NS'
    hr_model_map[2] = 'CIR'

    bond_model_map = {}
    bond_model_map[0] = 'RMV'
    bond_model_map[1] = 'RFV'
    
    hr_model   =   hr_model_map[int(hr_model_tmp)]
    bond_model = bond_model_map[int(bond_model_tmp)]

    str_elab_opt = hr_model + "," + rf_interp_tmp + "," + bond_model

    bf_options_elab = {}
    bf_options_elab['BondModel']   = bond_model
    bf_options_elab['HRateModel']  = hr_model
    bf_options_elab['MKTRef']      = 'de'
    bf_options_elab['MakeGraph']      = True
    bf_options_elab['MakeDump']      = False
    
    
    
    # ------------ lettura dati da foglio excel ------------------------------
    portfolio_xl        = readPortfolioFromXls(xla, curveDes, curvePos, nameSheet)
    
    bf_options_elab['RR'] = portfolio_xl.recoveryRate
    bf_options_elab['DataRef'] =  datetime.datetime.fromordinal(portfolio_xl.ref_date.toordinal())

    if bf_options_elab['RR'] == None:
        bf_options_elab['RR'] = 0.35
        portfolio_xl.recoveryRate = 0.35

        root = Tk()
        root.withdraw()
        msg0 = "Valorizzare correttamente il Recovery Rate.\n Assegnato il valore pari a %s   "%(bf_options_elab['RR'])
        tkMessageBox.showinfo("Attenzione!!", msg0)
        root.destroy()

    if (bf_options_elab['RR'] >0.99 or bf_options_elab['RR'] < 0):
        # significa che ho intercettato un errore!
        root = Tk()
        root.withdraw()
        msg0 = "Elaborazione non avviata: livello recovery rate (RR = %s) non compreso in [0,1]!" %(bf_options_elab['RR'])
        tkMessageBox.showinfo("Attenzione!!", msg0)
        
        root.destroy()
        return

    portfolio_xl.BondModel = bond_model
    portfolio_xl.HRateModel = hr_model

    curve_rf  = Curve()
    curve_cds = CdsCurve()

    interp_rf_model = curve_cds.mapCodeModelInv(rf_interp_tmp)

    # ------------ dati to download_curve risk free ------------------------------

    opt_curve_rf_download = {}

    opt_curve_rf_download['refDate'] = portfolio_xl.ref_date
    opt_curve_rf_download['valuta'] = portfolio_xl.curr
    opt_curve_rf_download['rating'] = portfolio_xl.rating
    opt_curve_rf_download['seniority'] = portfolio_xl.seniority

    opt_curve_rf_download['tipo_modello'] = interp_rf_model
    
    codeBenchList = ['%LS', '%DS', '%LFS', '%DFS']
    

    for codeTmp in codeBenchList:

        opt_curve_rf_download['codeSeg'] = codeTmp
        flag_loaded_rf = curve_cds.loadBenchDataFromDB_Bond(opt_curve_rf_download)
        if (flag_loaded_rf == 1):
            break
    
    if flag_loaded_rf == 0:
        # significa che ho intercettato un errore!
        root = Tk()
        root.withdraw()
        
        msg0 = "Curva benchmark associata al modello {} non presente alla data del {}, cambia modello o data!".format(interp_rf_model, portfolio_xl.ref_date)
        tkMessageBox.showinfo("Attenzione!", msg0)

        root.destroy()
        return

    if flag_loaded_rf == 2:
        # significa che ho intercettato un errore!
        root = Tk()
        root.withdraw()

        msg0 = "Il modello %s non gestito." % (interp_rf_model)
        tkMessageBox.showinfo("Attenzione!", msg0)

        root.destroy()
        return

    data_zc_rf = {}
    data_zc_rf['Model'] = curve_cds.bench_model
    data_zc_rf['ValoreNodo'] = curve_cds.bench_values
    data_zc_rf['MatDate'] = curve_cds.bench_dates
    data_zc_rf['prms'] = curve_cds.bench_prms
    data_zc_rf['DiscountFactors'] = curve_cds.bench_df_val
    
    """
    data_zc_rf = {}    
    data_zc_rf['Model'] =   2
    data_zc_rf['ValoreNodo'] =  [0.001820906, 0.001820906, 0.001786091, 0.001803993, 0.001835487, 0.001995183, 0.002285516, 0.002663186, 0.00312566, 0.003658126, 0.004242336, 0.004854089, 0.005471713, 0.006076757, 0.006667271, 0.007231011, 0.007768326, 0.008283799, 0.008764564, 0.009211182, 0.009629121, 0.01002267, 0.010395239, 0.010741152, 0.011061247, 0.011358885, 0.011636853, 0.011897481, 0.012139466, 0.012363825, 0.01257249, 0.012767119, 0.012949141]
    data_zc_rf['MatDate'] = [datetime.date(2014, 9, 30), datetime.date(2015, 8, 30), datetime.date(2016, 9, 30), datetime.date(2017, 9, 30), datetime.date(2018, 9, 30), datetime.date(2019, 8, 30), datetime.date(2020, 9, 30), datetime.date(2021, 9, 30), datetime.date(2022, 9, 30), datetime.date(2023, 8, 30), datetime.date(2024, 9, 30), datetime.date(2025, 9, 30), datetime.date(2026, 9, 30), datetime.date(2027, 8, 30), datetime.date(2028, 9, 30), datetime.date(2029, 9, 30), datetime.date(2030, 9, 30), datetime.date(2031, 8, 30), datetime.date(2032, 9, 30), datetime.date(2033, 9, 30), datetime.date(2034, 9, 30), datetime.date(2035, 8, 30), datetime.date(2036, 9, 30), datetime.date(2037, 9, 30), datetime.date(2038, 9, 30), datetime.date(2039, 8, 30), datetime.date(2040, 9, 30), datetime.date(2041, 9, 30), datetime.date(2042, 9, 30), datetime.date(2043, 8, 30), datetime.date(2044, 9, 30), datetime.date(2045, 9, 30), datetime.date(2046, 9, 30)]
    data_zc_rf['prms'] =  {'Dates': [datetime.date(2014, 9, 30)], 'const2': [25.003520717137043], 'const1': [4.3432558798248619], 'beta2': [0.026062917739813479], 'beta1': [0.044350833251596708], 'beta0': [-0.034445432562544648], 'beta3': [0.12488641360204013]}
    data_zc_rf['DiscountFactors'] = []
    """
    
    
    data_zc_infl = {}
    data_ts_infl = {}
    
    data_zc_infl['Model'] = []
    data_zc_infl['prms'] = []
    
    indx_ref = portfolio_xl.portfolio_anag[0]['Indicizzazione']


    if indx_ref == 'CPTFEMU':
        flag_load_infl_curve = curve_cds.loadInflCurveFromDB(opt_curve_rf_download)
    
        if flag_load_infl_curve == 0:
            # significa che ho intercettato un errore!
            root = Tk()
            root.withdraw()
            msg0 = "Curva inflazione non presente alla data del %s, cambia data!!" %(curve_rf.ref_date)
            tkMessageBox.showinfo("Attenzione!!", msg0)
    
            root.destroy()
            return
    
        flag_load_infl_ts = curve_cds.loadInflTSFromDB(opt_curve_rf_download)
    
        if flag_load_infl_ts == 0:
            # significa che ho intercettato un errore!
            root = Tk()
            root.withdraw()
            msg0 = "Serie storica inflazione non presente nel DB!!"
            tkMessageBox.showinfo("Attenzione!!", msg0)
    
            root.destroy()
            return


        data_zc_infl['MatDate'] = curve_cds.infl_curve_dates
        data_zc_infl['DiscountFactors'] = curve_cds.infl_curve_df_val
        
        zc_infl_times, zc_infl_rates = bf.fromDf2Rates(data_zc_infl['MatDate'], data_zc_infl['DiscountFactors']) 
        data_zc_infl['ValoreNodo'] = zc_infl_rates

        data_ts_infl['MatDate'] = curve_cds.infl_ts_dates
        data_ts_infl['Values'] =  curve_cds.infl_ts_val

    else:


        data_zc_infl['MatDate'] = []
        data_zc_infl['DiscountFactors'] = []
        data_zc_infl['ValoreNodo'] = []
    
        data_ts_infl['MatDate'] = []
        data_ts_infl['Values'] =  []

    #x0, x_bnd = bf.set_prms_for_fit(hr_model)
    
    x0    = par_x0
    x_bnd = par_bnd
    
    
    dictPortfolio,resDict = bf.fromXLSToBondFittingPortfolio(portfolio_xl.portfolio_anag)

    if resDict == True:
        return

    # ------------ elaborazione ---------------------------
 
    
    flag_elab, res_elab = bf.compute_bond_fitting(bf_options_elab, dictPortfolio, data_zc_rf, data_zc_infl, data_ts_infl, x0, x_bnd)

    

    
    if flag_elab == None:
        # significa che ho intercettato un errore!
        root = Tk()
        root.withdraw()
        msg = "Unable to perform Bond fitting!"
        tkMessageBox.showinfo("ERROR!", msg)
        root.destroy()
        return
    
    
    codice_curva = portfolio_xl.description
 
    # -------------- scrittura risultati bond fitting ---------------------------
    
    writeBondFittingRes1OnXls(portfolio_xl, xla, str_elab_opt, res_elab, codice_curva)
    writeBondFittingRes2OnXls(portfolio_xl, xla, str_elab_opt, res_elab, codice_curva)
    writeBondFittingRes3OnXls(portfolio_xl, xla, str_elab_opt, res_elab, codice_curva)
    writeBondFittingRes4OnXls(portfolio_xl, xla, str_elab_opt, res_elab, codice_curva)

    #------------------------------------------------------------

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
            msg = "Bootstrap results are on DB!"
            tkMessageBox.showinfo("YES WE CAN!", msg)
        else:
            ans = tkMessageBox.askquestion("Unable to save Bootstrap results because they're already on DB.", "DELETING... Are You Sure?", icon='warning')
            if ans =='yes':
                deletingDbCurves(codes)
                r, cd = saveZcDfOnDB(xla, nameSheet, des, pos, DF, ZC)
                if not r:
                    msg = "Something's wrong!!!!!"
                    tkMessageBox.showinfo("x@!#!", msg)

                else:
                    msg = "Bootstrap results are on DB!"
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
        root = Tk()
        msg = "Unable to save your selection...!! \n;)"
        tkMessageBox.showinfo("OOOOPS!", msg)
        root.destroy()
        
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

    n_c = len( W.curve)
    flag_plot = True

    if (n_c > 1):
        flag_plot = False

    for i in range(0, n_c):

        rangeStart = "B2"
        distance = 5

        xla       = xl_app()
        book      = xla.ActiveWorkbook
        s = book.Sheets(nameSheet)
        s.Activate()

        curveL     = readCurvesNames(xla,s,rangeStart,"o", distance)

        curveDes = W.curve[i]
        curvePos = int(W.pos[i]) 
        
    
        opt_boot_meth   = (str(W.new_window.variable1.get()).strip(""))[1]
        opt_rf_interp   = (str(W.new_window.variable6.get()).strip(""))[1]
        opt_hr_interp   = (str(W.new_window.variable7.get()).strip(""))[1]
    
        str_boot_opt = opt_boot_meth+","+opt_rf_interp + "," + opt_hr_interp
        data_opt                    = {}
        
        data_opt['hr_bootMethod']  = opt_boot_meth
        data_opt['bench_interp']   = opt_rf_interp
        data_opt['hr_interp']      = opt_hr_interp
        data_opt['Emittente']      = curveDes
        
        
        
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
        curve_xl.hr_model           = data_opt['hr_interp']
        
        
    
        # -------------- elaborazione ---------------------------
        boot_out     = curve_xl.bootstrap(data_opt, flag_plot)
        
        
        if boot_out == None:
            # significa che ho intercettato un errore!
            root = Tk()
            root.withdraw()
            msg = "Unable to perform CDS bootstrap!"
            tkMessageBox.showinfo("ERROR!", msg)
            root.destroy()
            return
        
        if (n_c > 1):
        
        
            if (i==0):
                flag_first = True
            else:
                flag_first = False

            writeCDSBootstrapRes1OnXls_surv_m(curve_xl, xla, str_boot_opt, boot_out, flag_first, codice_curva, i)
            s = xla.Cells.Columns.AutoFit()

            #writeCDSBootstrapRes1OnXls_md_m(curve_xl, xla, str_boot_opt, boot_out, codice_curva)
            #s = xla.Cells.Columns.AutoFit()

        
            #writeCDSBootstrapRes1OnXls_hr_m(curve_xl, xla, str_boot_opt, boot_out, codice_curva)
            #s = xla.Cells.Columns.AutoFit()

            writeCDSBootstrapRes2OnXls_zc_m(curve_xl, xla, str_boot_opt, boot_out, flag_first, codice_curva, i)
            s = xla.Cells.Columns.AutoFit()

            #writeCDSBootstrapRes2OnXls_py_m(curve_xl, xla, str_boot_opt, boot_out, codice_curva)
            #s = xla.Cells.Columns.AutoFit()
        

        else:

            writeCDSBootstrapRes1OnXls(curve_xl, xla, str_boot_opt, boot_out, codice_curva)
            s = xla.Cells.Columns.AutoFit()
        
            writeCDSBootstrapRes2OnXls(curve_xl, xla, str_boot_opt, boot_out, codice_curva)
            s = xla.Cells.Columns.AutoFit()

# ==========================================
# punto d'ingresso per CALIBRAZIONE
# ==========================================

from sc_elab.excel_hook.W_calibration import W_calib_models
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from sc_elab.core.funzioni_calibrazioni import *
from sc_elab.excel_hook.xls_Calibration import *

@xl_func
def calibration_from_xls(control):

    nameSheet = nameSheetCalib

    xla = xl_app()
    book = xla.ActiveWorkbook
    wbName = str(book.FullName)
    book.Save()

    root = Tk()
    root.iconbitmap(default= imp.find_module('sc_elab')[1] + r'\\excel_hook\\fig\\icona.ico')

    allSheetInBook = allSheet(book)

    # -------------- controllo l'esistenza del foglio  ----------------
    if nameSheet in allSheetInBook:
        s = book.Sheets(nameSheet)
        s.Activate()
    else:
        msg = "Missing input sheet(%s) for Calibration in your workbook... \nNothing to do for me!" %nameSheet
        tkMessageBox.showinfo("Warning!", msg)
        root.destroy()
        return
    # -----------------------------------------------------------------

    W = W_calib_models(master = root, nameWorkbook= wbName, nameWorksheet=nameSheet)
    root.mainloop()

    if W.res > 0:
        W1 = W.newWindow

        if W1.res > 0:
            loss_function_type_tmp = W1.loss_function_type.get()

            if loss_function_type_tmp == 1:
                loss_function_type_power = 2
                loss_function_type_absrel = 'abs'

            elif  loss_function_type_tmp == 2:
                loss_function_type_power = 2
                loss_function_type_absrel = 'rel'

            elif loss_function_type_tmp == 3:
                loss_function_type_power = 1
                loss_function_type_absrel = 'abs'

            elif loss_function_type_tmp == 4:
                loss_function_type_power = 1
                loss_function_type_absrel = 'rel'


            model_dict = W1.param_dict

            model = W.model.get()
            type_data = W1.set_mkt_ts.get()

            if type_data == 'MKT':

                if (model in ['CIR','VSCK']):
                    mkt_value, mkt_to_fit, type_cap = preProcessignCurve(W1.CurveChosen)

                    x0_m  = []
                    x_bnd = []

                    l_bound = []
                    h_bound = []

                    for p_name in W1.params_names:
                        x0_m.append(float(W1.param_dict[p_name]['sv']))
                        x_bnd.append([float(W1.param_dict[p_name]['min']), float(W1.param_dict[p_name]['max'])])
                        l_bound.append(float(W1.param_dict[p_name]['min']))
                        h_bound.append(float(W1.param_dict[p_name]['max']))

                    if W.model.get() == 'CIR':
                        ff = minimize(loss_zc_model_cir, args = (mkt_to_fit ,loss_function_type_power,loss_function_type_absrel), x0 = x0_m, method='TNC', bounds=x_bnd)
                    else:
                        ff = minimize(loss_zc_model_vsck, args = (mkt_to_fit ,loss_function_type_power,loss_function_type_absrel), x0 = x0_m, method='TNC', bounds=x_bnd)

                    # creo la lista dei risultati ottimali
                    list_model_params_opt = []
                    list_model_params_opt.append(ff.x[0])
                    list_model_params_opt.append(ff.x[1])
                    list_model_params_opt.append(ff.x[2])
                    list_model_params_opt.append(ff.x[3])

                    if W.model.get() == 'CIR':
                        mkt_value['VALUE_OPT'] = compute_zc_cir_rate(list_model_params_opt, mkt_value["TIME"])
                    else:
                        mkt_value['VALUE_OPT'] = compute_zc_vsck_rate(list_model_params_opt, mkt_value["TIME"])


                    # converto i risultati in composto nel caso in cui in input lo siano
                    if type_cap == 'CMP':
                        mkt_value['VALUE_OPT'] = fromContinuousToCompost(mkt_value['VALUE_OPT'])

                    # calcolo il chi quadro
                    chi2 = computeCHI2(mkt=mkt_value["VALUE"], mdl=mkt_value['VALUE_OPT'])

                    # scrivo su foglio Excel
                    writeCalibrationResOnXls(type_data = type_data,
                                             model = model,
                                             W_class = W1,
                                             xla = xla,
                                             chi2 = chi2,
                                             opt_dict = list_model_params_opt,
                                             res = mkt_value,
                                             capitalization_type = type_cap)

                    # produco il grafico
                    mkt_value.set_index('TIME',inplace=True)
                    mkt_value.plot(style=['o', '-'])
                    plt.title('Calibration results')
                    plt.xlabel('Time')
                    plt.ylabel('Rate')
                    plt.show()

                if W1.NameOption.get() != "":
                    opt_total = W1.OptionChosen

            else:

                if (model in ['CIR','VSCK']):
                    # pre processing time series
                    mkt_to_fit = preProcessignTimeSeries(df = W1.TSChosen, dt_min= W1.TS_dateMIN , dt_max= W1.TS_dateMAX)

                    # creazione del set di parametri
                    x0_m  = []
                    x_bnd = []

                    l_bound = []
                    h_bound = []

                    for p_name in W1.params_names:
                        x0_m.append(float(W1.param_dict[p_name]['sv']))
                        x_bnd.append([float(W1.param_dict[p_name]['min']), float(W1.param_dict[p_name]['max'])])
                        l_bound.append(float(W1.param_dict[p_name]['min']))
                        h_bound.append(float(W1.param_dict[p_name]['max']))

                    # ottimizzazione
                    if W.model.get() == 'CIR':
                        ff = minimize(mle_cir, args=(mkt_to_fit), x0 = x0_m, method='TNC', bounds=x_bnd)
                    else:
                        ff = minimize(mle_vsck, args=(mkt_to_fit), x0 = x0_m, method='TNC', bounds=x_bnd)

                    # creo la lista dei risultati ottimali
                    list_model_params_opt = []
                    list_model_params_opt.append(mkt_to_fit.loc[0,'VALUE'])
                    list_model_params_opt.append(ff.x[1])
                    list_model_params_opt.append(ff.x[2])
                    list_model_params_opt.append(ff.x[3])

                    # creazione dei risultati
                    if W.model.get() == 'CIR':
                        model_value = generate_cir_perc(params = list_model_params_opt, data_to_fit = mkt_to_fit)
                    else:
                        model_value = generate_vsck_perc(params = list_model_params_opt, data_to_fit = mkt_to_fit)

                    # calcolo il chi quadro
                    chi2 = computeCHI2(mkt=model_value["VALUE"], mdl=model_value['MODEL VALUE MEAN'])

                    # scrivo su foglio Excel
                    writeCalibrationResOnXls(type_data = type_data,model = model, W_class = W1, xla = xla, chi2 = chi2, opt_dict = list_model_params_opt, res = model_value)

                    # produzione del grafico
                    model_value.set_index('DATE',inplace=True)
                    model_value.plot(style=['bo','r-', 'k--', 'k--'])
                    plt.title('Calibration results')
                    plt.xlabel('Date')
                    plt.ylabel('Value')
                    plt.show()

# ==========================================
# punto d'ingresso per Template Calibration
# ==========================================

@xl_func
def template_elaborate_calibration(control):

    xla = xl_app()
    wb = xla.ActiveWorkbook

    allSheetInBook = allSheet(wb)

    # -------------- controllo l'esistenza del foglio  ----------------
    if not (nameSheetCalib in allSheetInBook):
        s = wb.Sheets.Add()
        s.Name = nameSheetCalib
    else:
        s = wb.Sheets(nameSheetCalib)
        s.Activate()
    # -----------------------------------------------------------------

    '''
    possibilita di scegliere la tipologia di template
    
    root = Tk()
    root.iconbitmap(default= imp.find_module('sc_elab')[1] + r'\\excel_hook\\fig\\icona.ico')
    root.mainloop()
    '''

    writeTemplateCalibration(xla = xla, nameSheet = nameSheetCalib)

# ==========================================
# punto d'ingresso per Caricamento Dati
# ==========================================


from loading_data import test_load_nuovi_dati


@xl_func
def caricamento_dati(control):

    xla = xl_app()
    book = xla.ActiveWorkbook

    wbName = str(book.FullName)
    book.Save()

    #pathwb = 'C:/Users/scalambrinm/workspace/scenario/sc_elab/core/input/files_caricamento_datastream/Curva_Depositi_EUR_ICAP.xlsx'
    test_load_nuovi_dati(wbName)


# ==========================================
# punto d'ingresso per scarico Swaptions
# ==========================================

import pandas as pd

from W_download_Swaptions import W_SwaptionsDate, W_SwaptionsOptionsPrint
from xls_download_Swaptions import write_Swaptions
from DEF_intef import nameSheetScaricoSwaption
from sc_elab.excel_hook.connection import Connection

@xl_func
def writeSwaptions(control):

    nameSheet = nameSheetScaricoSwaption
    xla = xl_app()
    book = xla.ActiveWorkbook

    # -------------- controllo l'esistenza del foglio di input  ----------------
    try:
        s = book.Sheets(nameSheet)
        s.Activate()
    except:
        s = book.Sheets.Add()
        s.Name = nameSheet
    # -------------- apro la finestra di input della scelta  ----------------------------

    con = Connection()
    #necessario per inizializzare il db
    cursor = con.db_data()

    root = Tk()
    app = W_SwaptionsDate(root)
    root.mainloop()

    ref_date = datetime.date(day=int(app.date[-2:]), month=int(app.date[5:7]), year=int(app.date[:4]))

    qry_to_execute= '''
                 SELECT DProCFS.Tenor, DProCFS.MaturityInt, DProTS_master.ValoreMid
                 FROM DProCFS, DProTS_master
                 WHERE DProTS_master.BloombergTicker = DProCFS.BloombergTicker
                 AND(DProCFS.TipoDato = 'VSwaption') 
                 and DProTS_master.Data= '%s' ''' %(ref_date)

    res = pd.read_sql(qry_to_execute, con.db)


    #query che permette di trovare per una certa data il contributor e la currency. Attualmente
    #lo scarico permette di ottenere una sola superficie di swaption per contributor e per currency
    #pertanto attualmente non e prevista la possibilita di avere due superfici diverse alla stessa data
    qry_to_execute= '''
                    SELECT distinct DProCFS.Contributor,DProCFS.Currency from DProCFS,DProTS_master 
                    WHERE DProTS_master.BloombergTicker = DProCFS.BloombergTicker
                    AND (DProCFS.TipoDato = 'VSwaption') 
                    and DProTS_master.Data= '%s' ''' %(ref_date)

    res2 =pd.read_sql(qry_to_execute, con.db)
    if res2.shape[0] > 1:
        root = Tk()
        tkMessageBox.showinfo("Warning!", 'I dati selezionati hanno molteplici Contributor o Currency.')
        root.mainloop()
        return

    currency = res2['Currency'][0]
    contributor = res2['Contributor'][0]

    # Interrogo il campo deltaBp per capire se le Swaption sono shiftate. Qualora il campo deltaBp
    # per tutte le componenti della superficie di Swaption ad una certa data sia identico e pari
    # a None, allora restituisce No Shifted

    qry_to_execute = '''
                        SELECT distinct DProCFS.deltaBp from DProCFS,DProTS_master 
                        WHERE DProTS_master.BloombergTicker = DProCFS.BloombergTicker
                        AND (DProCFS.TipoDato = 'VSwaption') 
                        and DProTS_master.Data= '%s' ''' % (ref_date)

    res3 = pd.read_sql(qry_to_execute, con.db)


    if not res3['deltaBp'][0] and res3.shape[0] == 1 :
        tipo_modello = "No Shifted"
    else:
        tipo_modello = '''Shifted'''

    root = Tk()
    app = W_SwaptionsOptionsPrint(root)
    root.mainloop()

    write_Swaptions(xla, res, ref_date, currency, contributor, tipo_modello, option_print = app.print_type.get() )
    s = xla.Cells.Columns.AutoFit()


# ==========================================
# punto d'ingresso per TEMPLATE
# ==========================================

from sc_elab.excel_hook.createTemplate import writeTemplate, W_template, ask_question
from sc_elab.excel_hook.xls_utils import allSheet
from sc_elab.core.Tipologia_curva_dizionario import *
from sc_elab.core.db_data_structure_v0 import table_dict, table_dict_Dati
from Tkinter import *
import tkMessageBox
import imp

@xl_func
def create_Template(control):

    xla = xl_app()
    book = xla.ActiveWorkbook

    root = Tk()
    root.title('Menu Template')
    root.iconbitmap(default= imp.find_module('sc_elab')[1] + r'\\excel_hook\\fig\\icona.ico')
    app = W_template(root)
    root.mainloop()

    if app.template == False:
        return

    tmp_array = app.template
    allSheetInBook = allSheet(book)

    for tmp in tmp_array:
        if tmp == 'Dati':
            nameSheet = 'Dati'
            table = table_dict_Dati

            if not(nameSheet in allSheetInBook):
                writeTemplate(xla, book, nameSheet, table)
            else:
                answer = ask_question('Delete sheet', 'Vuoi eliminare il foglio %s ?' % nameSheet)
                if answer == 'yes':

                    t = book.Worksheets(nameSheet).Delete()
                    writeTemplate(xla, book, nameSheet, table)


        elif tmp == 'Bond_master':

            nameSheet = 'Dati'
            table = table_dict[tmp]

            if not(nameSheet in allSheetInBook):
                writeTemplate(xla, book, nameSheet, table)
            else:
                answer = ask_question('Delete sheet', 'Vuoi eliminare il foglio %s ?' % nameSheet)
                if answer == 'yes':

                    t = book.Worksheets(nameSheet).Delete()
                    # il foglio non e' stato eliminato
                    writeTemplate(xla, book, nameSheet, table)

        else:
            nameSheet = 'Anagrafica_' + tmp
            table = table_dict[corrisp_tabella[tmp]]

            if not(nameSheet in allSheetInBook):
                writeTemplate(xla, book, nameSheet, table)
            else:
                answer = ask_question('Delete sheet', 'Vuoi eliminare il foglio %s ?' % nameSheet)
                if answer == 'yes':

                    t = book.Worksheets(nameSheet).Delete()
                    writeTemplate(xla, book, nameSheet, table)





# ========================================================
# punto d'ingresso per ELABORATE TRANSITION MATRIX
# ========================================================

from sc_elab.excel_hook.W_quarterly_matrix import W_matrix_trim, W_select_matrix, writeQuarterlyMatrixResOnXls, W_dim_matrix, writeTemplateQuarterlyMatrixInput
from sc_elab.excel_hook.W_calibration import readSheetObject, readFeaturesObject

from sc_elab.excel_hook.xls_utils import allSheet
from DEF_intef import nameSheetElabMatrix, nameSheetElabMatrixResult

from Tkinter import *
import tkMessageBox
import imp

@xl_func
def elaborate_quarterly_matrix(control):

    xla = xl_app()
    wb = xla.ActiveWorkbook

    allSheetInBook = allSheet(wb)

    # -------------- controllo l'esistenza del foglio  ----------------
    if not (nameSheetElabMatrix in allSheetInBook):
        s = wb.Sheets.Add()
        s.Name = nameSheetElabMatrix
    else:
        s = wb.Sheets(nameSheetElabMatrix)
        s.Activate()
    # -----------------------------------------------------------------

    wb.Save()

    objectOnSheetDictionary = readSheetObject(workbook_path = str(wb.FullName), sheet_name = nameSheetElabMatrix )


    if objectOnSheetDictionary[1].empty:
        root = Tk()
        root.withdraw()
        tkMessageBox.showinfo("Warning!", 'Nessuna matrice presente nel foglio!')
        root.destroy()
        return

    else:
        objectOnSheet = readFeaturesObject(objectOnSheetDictionary)
        tmpMatrix = objectOnSheet.loc[objectOnSheet.TypeObject == 'Matrix', 'Name'].tolist()


    root = Tk()
    root.iconbitmap(default= imp.find_module('sc_elab')[1] + r'\\excel_hook\\fig\\icona.ico')
    app1 = W_select_matrix(master = root, ListMatrix = tmpMatrix)
    root.mainloop()

    tmp = objectOnSheet.loc[objectOnSheet.Name == app1.MatrixNameChosen, 'keys'].values[0]
    MatrixChosen = objectOnSheetDictionary[tmp].dropna(axis = 0).astype('float').copy()

    if MatrixChosen.shape[0] != MatrixChosen.shape[1]:
        root = Tk()
        root.withdraw()
        tkMessageBox.showinfo("Warning!", "La matrice deve essere quadrata. Errore!")
        root.destroy()
        return

    root = Tk()
    root.iconbitmap(default= imp.find_module('sc_elab')[1] + r'\\excel_hook\\fig\\icona.ico')
    app2 = W_matrix_trim(master = root, matrix = MatrixChosen)
    root.mainloop()

    if app2.flag_save == True:
        writeQuarterlyMatrixResOnXls(xla = xla, W_class = app2, book = wb, nameSheet = nameSheetElabMatrixResult, nameMatrix = app1.MatrixNameChosen)


@xl_func
def template_elaborate_matrix(control):

    xla = xl_app()
    wb = xla.ActiveWorkbook

    allSheetInBook = allSheet(wb)

    # -------------- controllo l'esistenza del foglio  ----------------
    if not (nameSheetElabMatrix in allSheetInBook):
        s = wb.Sheets.Add()
        s.Name = nameSheetElabMatrix
    else:
        s = wb.Sheets(nameSheetElabMatrix)
        s.Activate()
    # -----------------------------------------------------------------

    root = Tk()
    root.iconbitmap(default= imp.find_module('sc_elab')[1] + r'\\excel_hook\\fig\\icona.ico')
    app = W_dim_matrix(master = root)
    root.mainloop()

    writeTemplateQuarterlyMatrixInput(xla = xla, nameSheet = nameSheetElabMatrix, dimMatrix = app.dimMatrix.get())

# ===================================================================
# punto d'ingresso per scarico volatilita' implicite nei Cap e Floor
# ===================================================================

from W_Bootstrap_VolCapFloor import W_VolCapFloorDate, W_CurrencySelection
from xls_Bootstrap_VolCapFloor import write_VolCapFloor
from DEF_intef import nameSheetCapFloorVolatilities
from sc_elab.excel_hook.connection import Connection

@xl_func
def CapFloor_BVol_fromDBtoXls(control):

    nameSheet = nameSheetCapFloorVolatilities
    xla = xl_app()
    book = xla.ActiveWorkbook

    # -------------- controllo l'esistenza del foglio di input  ----------------
    try:
        s = book.Sheets(nameSheet)
        s.Activate()
    except:
        s = book.Sheets.Add()
        s.Name = nameSheet
    # -------------- apro la finestra di input della scelta  ----------------------------

    con = Connection()
    # necessario per inizializzare il db
    cursor = con.db_data()

    root = Tk()
    app = W_VolCapFloorDate(root)
    root.mainloop()

    ref_date = datetime.date(day=int(app.date[-2:]), month=int(app.date[5:7]), year=int(app.date[:4]))

    # query che permette di trovare per una certa data la currency
    qry_to_execute = '''
                    SELECT distinct DProCFS.Currency from DProCFS,DProTS_master
                    WHERE DProTS_master.BloombergTicker = DProCFS.BloombergTicker
                    AND (DProCFS.TipoDato = 'VCapFloor')
                    and DProTS_master.Data= '%s' ''' % (ref_date)

    currencies = pd.read_sql(qry_to_execute, con.db)['Currency'].tolist()

    # query che permette di trovare per una certa data il contributor
    qry_to_execute = '''
                        SELECT distinct DProCFS.Contributor from DProCFS,DProTS_master
                        WHERE DProTS_master.BloombergTicker = DProCFS.BloombergTicker
                        AND (DProCFS.TipoDato = 'VCapFloor')
                        and DProTS_master.Data= '%s' ''' % (ref_date)

    contributors = pd.read_sql(qry_to_execute, con.db)['Contributor'].tolist()

    # La seguente selezione genera un errore, ma tutto funziona
    root = Tk()
    app = W_CurrencySelection(root,currencies=currencies,contributors=contributors)
    root.mainloop()


    currency = app.curr.get()
    contributor = app.cont.get()

    # query di scarico dei dati da stampare su foglio Excel
    qry_to_execute = '''
                 SELECT DProCFS.MaturityInt, DProCFS.Strike, DProTS_master.ValoreMid
                 FROM DProCFS, DProTS_master
                 WHERE DProTS_master.BloombergTicker = DProCFS.BloombergTicker
                 AND(DProCFS.TipoDato = 'VCapFloor')
                 and DProCFS.Contributor = '%s'
                 and DProCFS.Currency = '%s' 
                 and DProTS_master.Data= '%s' 
                 order by DProCFS.MaturityInt ASC''' % (contributor, currency, ref_date)

    res = pd.read_sql(qry_to_execute, con.db)

    if res['ValoreMid'].shape[0]==0:
        tkMessageBox.showwarning('Attenzione', 'Non si dispone di dati per la coppia Contributor-Currency selezionata')

    # Interrogo il campo deltaBp per capire se i Cap e Floor sono shiftati. Qualora il campo deltaBp
    # per tutte le componenti della superficie di volatilita Cap e Floor ad una certa data sia identico e pari
    # a None, allora restituisce No Shifted

    qry_to_execute = '''
                        SELECT distinct DProCFS.deltaBp from DProCFS,DProTS_master 
                        WHERE DProTS_master.BloombergTicker = DProCFS.BloombergTicker
                        AND (DProCFS.TipoDato = 'VCapFloor') 
                        and DProTS_master.Data= '%s' ''' % (ref_date)

    res3 = pd.read_sql(qry_to_execute, con.db)

    if not res3['deltaBp'][0] and res3.shape[0] == 1:
        tipo_modello = "No Shifted"
        shift='0'
    elif res3.shape[0] > 1:
        root = Tk()
        root.withdraw()
        tkMessageBox.showwarning('Attenzione', 'Sono presenti shift diversi')
        shift = '-999'
        root.destroy()
    else:
        tipo_modello = 'Shifted'
        shift = res3['deltaBp'][0]

    # root = Tk()
    # app = W_SwaptionsOptionsPrint(root)
    # root.mainloop()

    write_VolCapFloor(xla, res, ref_date, currency, contributor, tipo_modello, shift)
    s = xla.Cells.Columns.AutoFit()



# ====================================================================
# punto d'ingresso per Bootstrap volatilita' implicite nei Cap e Floor
# ====================================================================

import pandas as pd

from W_Bootstrap_VolCapFloor import Bootstrap_BVol_menu
from xls_Bootstrap_VolCapFloor import writeBootstrapVolOnXls, readFeaturesDiscCurve
from DEF_intef import nameSheetCapFloorVolatilities, nameSheetBootstrap
from W_calibration import readSheetObject, readFeaturesObject
from sc_elab.core.anagrafica_dati import MaturityFromStringToYear
from sc_elab.core.funzioni_boot_cap_floor import Bootstrap_CapFloor_ATM
import tkMessageBox

@xl_func
def BootstrapCapFloorVol_on_xls(control):

    nameSheetVols = nameSheetCapFloorVolatilities
    nameSheetCurve= nameSheetBootstrap
    xla = xl_app()
    book = xla.ActiveWorkbook

    allSheetInBook = allSheet(book)

    # -------------- controllo l'esistenza del foglio di input per le volatilita' e scarico i dati -------------
    if nameSheetVols not in allSheetInBook:
        msg = "Missing input sheet(%s) for Bootstrap in your workbook... \nNothing to do for me!" % nameSheetVols
        root = Tk()
        tkMessageBox.showwarning("Warning!", msg)
        root.destroy()
        return

    # -------------- controllo l'esistenza del foglio di input per la curva e scarico i dati -------------
    if nameSheetCurve not in allSheetInBook:
        msg = "Missing input sheet(%s) for Bootstrap in your workbook... \nNothing to do for me!" % nameSheetCurve
        root = Tk()
        tkMessageBox.showwarning("Warning!", msg)
        root.destroy()
        return

    # scarico i dati
    volsdata = readSheetObject(str(book.FullName), nameSheetVols)
    discount_curves = readSheetObject(str(book.FullName), nameSheetCurve)

    # leggo le caratteristiche identificative di questi oggetti
    volsdata_list = readFeaturesObject(volsdata)
    disc_curves_list = readFeaturesDiscCurve(discount_curves)

    # creo le liste da passare alle combobox per la scelta di curva e volatilita'
    volsdata_choices = []
    for i in volsdata_list.index:
        if volsdata_list.loc[i, 'TypeObject'] == 'Option':
            volsdata_tmp = str(volsdata_list.loc[i, 'keys']) + ' ' + str(volsdata_list.loc[i, 'Name'])
            volsdata_choices.append(volsdata_tmp)

    disc_curves_choices = []
    for i in disc_curves_list.index:
        if disc_curves_list.loc[i, 'TypeObject'] == 'Discount Curve':
            discount_curve_tmp = str(disc_curves_list.loc[i, 'keys']) + ' ' + str(disc_curves_list.loc[i, 'Name'])
            disc_curves_choices.append(discount_curve_tmp)

    # controllo di avere i dati per il Bootstrap
    if len(volsdata_choices) == 0 or len(disc_curves_choices) == 0:
        root = Tk()
        tkMessageBox.showwarning('Warning', 'Missing data to perform the bootstrap...I cannot do anything!')
        root.destroy()
        return

    # apro la finestra di selezione dei dati

    choices = Bootstrap_BVol_menu(volsdata_choices, disc_curves_choices)
    if choices[0] == 0:
        root=Tk()
        tkMessageBox.showwarning('Salutation', 'Au-revoir')
        root.destroy()
        return

    selected_disc_curve = int(choices[1][:1])
    selected_vols = int(choices[2][:1])

    # Eseguo un controllo sulle date di riferimento
    refdate_vol   = volsdata[selected_vols].loc[(volsdata[selected_vols].loc[:,0]=='Date ref'),1].values[0]
    refdate_curve = discount_curves[selected_disc_curve].loc[(discount_curves[selected_disc_curve].loc[:,0]=='Date Ref'),1].values[0]
    if refdate_vol != refdate_curve:
        root=Tk()
        tkMessageBox.showinfo('Attenzione', 'Curve and volatilities reference date are not the same')
        root.destroy()

    # seleziono i dati che andranno in input alla funzione di bootstrap e li formatto come array di floats

    shift = float(volsdata[selected_vols].loc[volsdata[selected_vols][0] == 'Shift', 1])

    volsdata_noint = volsdata[selected_vols].loc[(volsdata[selected_vols].loc[:, 3] == 'Y'), [0, 2]]
    volatilities = {}
    volatilities['Maturities'] = volsdata_noint.loc[:, 0].map(MaturityFromStringToYear).values.astype(float)
    volatilities['Volatilities'] = np.divide(volsdata_noint.loc[:, 2].values.astype(float),100)

    curve = discount_curves[selected_disc_curve].loc[(discount_curves[selected_disc_curve].loc[:, 2] == 'Y'), [0, 1]]
    discount = {}
    discount['discount times'] = curve.loc[:, 0].values
    discount['discount factors'] = curve.loc[:, 1].values.astype(float)

    # converto le date dei fattori di sconto in intervalli in termini di giorni (intanto provo)
    discount['discount times'] = discount['discount times'] - discount['discount times'][0]
    discount['discount times'] = np.array([(d.days) / 365.2425 for d in discount['discount times']])

    bootstrapped_volatilities = Bootstrap_CapFloor_ATM(shift, discount, volatilities)

    writeBootstrapVolOnXls(xla, bootstrapped_volatilities, volsdata[selected_vols], discount_curves[selected_disc_curve])