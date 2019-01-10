



from sc_elab.core import funzioni_bond_fitting as bf
from sc_elab.core import funzioni_base as fb

import time
import datetime

def set_prms_for_fit(model_fit):
    
    

    
    if (model_fit == 'SVE'): #->SVE
    
        #bound_min_sve = [0.0001,  0.0001, -10.00, -10.050, -10.00, -10.00]

        bound_min = [0.0001,  0.0001, -10.00, -10.00, -10.00, -10.00]
        bound_max = [10.0,     50.0,  10.03,   10.0,   10.5,   10.0]
        x0  = [1.0,       10.0,   0.03,   0.03,   0.03,   0.03]
        n_par = 6


    elif (model_fit == 'CIR'): #->CIR

        bound_min = [  -0.1,  0.1, 0.001, 0.001]
        bound_max = [ 10.00, 10.0, 10.00, 1.000]
        x0  = [0.0001,  1.0, 0.015,  0.01]
        n_par = 4


    elif (model_fit == 'NS'): #->NS

        bound_min = [0.0001,  -10.0, -10.0, -10.0]
        bound_max = [100,     +10.0, +10.0, +10.0]
        x0        = [5.0001,  0.03, 0.03, 0.03]
        n_par = 4

    
    x_bnd = []
    for i in range(0, n_par):
    
        bndTmp = []
        
        b_min = bound_min[i]
        b_max = bound_max[i]
        
        bndTmp.append(b_min)
        bndTmp.append(b_max)
        
        x_bnd.append(bndTmp)


    return x0, x_bnd

    



if __name__ == "__main__":

    t_start = time.time()
    
    #---------------------------------------------------------------------------        
    #----------------------- Carico dati di mercato ----------------------------
    #---------------------------------------------------------------------------    

    """
    dataRef = datetime.datetime(2011,12,30)
    recovery_rate = 0.4

    #input_rf          = "input_test/bond_fitting_data/test_rf_31_12_11_01.txt"    
    #C:\Users\monacoa\workspace\scenario\sc_elab\core\input_test\dond_fitting_data
    
    input_rf          = "input_test/bond_fitting_data/test_rf_31_12_11_01.txt"    
    input_infl        = "input_test/bond_fitting_data/test_infl_v0.txt"
    inflation_ts_file = "input_test/bond_fitting_data/inflation_ts_v2.txt"    

    portfolio_file    = "input_test/bond_fitting_data/TEST_PTF_ITA_30_06_12_FIXED_ANAG_V2.txt"
    #portfolio_file    = "input_test/bond_fitting_data/TEST_PTF_ITA_30_06_12_FLOAT_ANAG_V2.txt"
    #portfolio_file    = "input_test/bond_fitting_data/TEST_PTF_ITA_30_06_12_ZC_ANAG_V2.txt"
    """
    
    
    
    # ----------------- TEST REF V2 ----------------------------------------------

    """
    dataRef = datetime.datetime(2017,03,07)
    recovery_rate = 0.4
    inflation_ts_file  = "input_test/bond_fitting_data/inflation_ts_v2.txt"    
    input_infl         = "input_test/bond_fitting_data/test_infl_v0.txt"
    input_rf           = "input_test/bond_fitting_data/test_rf_07_03_17_flat_0.01.txt"
    #portfolio_file     = "input_test/bond_fitting_data/TEST_BTP_ITA_07_03_17_FLAT_0.03_v2.txt"
    portfolio_file    = "input_test/bond_fitting_data/TEST_BTPi_ITA_07_03_17_FLAT_0.03_v2.txt" # errore
    """

    #-------------------- TEST REAL DATA BTP --------------------------------------

    """
    dataRef = datetime.datetime(2017,03,07)
    recovery_rate       = 0.4
    inflation_ts_file     = "input_test/bond_fitting_data/inflation_ts_v2.txt"    
    input_rf              = "input_test/bond_fitting_data/test_rf_07_03_17.txt"
    input_infl            = "input_test/bond_fitting_data/test_infl_22_12_16.txt"
    portfolio_file        = "input_test/bond_fitting_data/TEST_BTP_ITA_07_03_17_v0.txt"
    """
    

    #-------------------- TEST REAL DATA BTPi --------------------------------------
    
    """
    dataRef = datetime.datetime(2014,9,30)
    recovery_rate       = 0.3
    inflation_ts_file = "input_test/bond_fitting_data/inflation_ts_v2.txt"    
    input_rf          = "input_test/bond_fitting_data/test_rf_30_09_2014_v0.txt"
    input_infl        = "input_test/bond_fitting_data/test_infl_22_12_16.txt"
    portfolio_file    = "input_test/bond_fitting_data/TEST_BTPi_ITA_30_09_14_v0.txt"
    """
    
    
    # ---------------------TEST ATLANTIA V0 ----------------------------------
    """
    dataRef = datetime.datetime(2015,05,15)
    recovery_rate = 0.3

    input_rf          = "input_test/bond_fitting_data/test_rf_12_05_2015_v0.txt"    

    input_infl        = "input_test/bond_fitting_data/test_infl_v0.txt"
    inflation_ts_file = "input_test/bond_fitting_data/inflation_ts_v2.txt"    

    portfolio_file    = "input_test/bond_fitting_data/TEST_ATLANTIA_V0.txt"
    """
    #-----------------------------------------------------------------------------------------
    
    # ---------------------TEST UBS V0 ----------------------------------
    
    """
    dataRef = datetime.datetime(2012,12,31)
    recovery_rate = 0.3

    input_rf          = "input_test/bond_fitting_data/test_v1/test_rf_31_12_12_UBS.txt"    

    input_infl        = "input_test/bond_fitting_data/test_v1/test_infl_v0.txt"
    inflation_ts_file = "input_test/bond_fitting_data/test_v1/inflation_ts_v2.txt"    

    portfolio_file    = "input_test/bond_fitting_data/test_v1/TEST_UBS_V0.txt"
    """
    
    #-----------------------------------------------------------------------------------------

    dataRef = datetime.datetime(2017,11,21)

    recovery_rate = 0.331

    #input_rf          = "input_test/bond_fitting_data/test_v1/test_rf_21_11_17_XXX.txt"    
    #input_rf          = "input_test/bond_fitting_data/test_v1/test_rf_21_11_17_XXX.txt"    
    input_rf          = "input_test/bond_fitting_data/test_v1/test_rf_21_11_17_YYY.txt"    

    input_infl        = "input_test/bond_fitting_data/test_v1/test_infl_v0.txt"
    inflation_ts_file = "input_test/bond_fitting_data/test_v1/inflation_ts_v2.txt"    

    #portfolio_file    = "input_test/bond_fitting_data/test_v1/TEST_PTF_MEDIOBANCA_V0.txt"
    portfolio_file    = "input_test/bond_fitting_data/test_v1/TEST_PTF_MEDIOBANCA_S.txt"
    
    
    
    #-------------------- TEST REAL DATA BTPi --------------------------------------

    
    """
    dataRef = datetime.datetime(2016,12,16)
    recovery_rate       = 0.3
    inflation_ts_file = "input_test/bond_fitting_data/inflation_ts_v2.txt"    
    input_rf          = "input_test/bond_fitting_data/test_rf_22_12_16.txt"
    input_infl        = "input_test/bond_fitting_data/test_infl_22_12_16.txt"
    portfolio_file    = "input_test/bond_fitting_data/TEST_BTPi_ITA_22_12_16_v0.txt"
    """
    
    
    """
    # ---------------------TEST CAIXA V0 ----------------------------------
    dataRef = datetime.datetime(2014,9,30)
    recovery_rate = 0.3

    #input_rf          = "input_test/bond_fitting_data/test_rf_31_12_11_01.txt"    
    input_rf          = "input_test/bond_fitting_data/test_rf_caixa_30_09_2014_v0.txt"    

    input_infl        = "input_test/bond_fitting_data/test_infl_v0.txt"
    inflation_ts_file = "input_test/bond_fitting_data/inflation_ts_v2.txt"    

    portfolio_file    = "input_test/bond_fitting_data/TEST_CAIXA_V0.txt"

    #-----------------------------------------------------------------------------------------
    """


    #prms_file          = "input_test/bond_fitting_data/test_v1/model_params_sve_v0.txt"
    #prms_file          = "input_test/bond_fitting_data/model_params_cir_v0.txt"
    #prms_file          = "input_test/bond_fitting_data/model_params_ns_v0.txt"
    #prms_file          = "input_test/bond_fitting_data/model_params_sve_v0.txt"

    prms_file          = "input_test/bond_fitting_data/test_v1/model_params_sve_XXX.txt"
    #prms_file          = "input_test/bond_fitting_data/test_v1/model_params_sve_YYY.txt"
    
    #--------------------------------------------------------------------------------
    
    #out_file         = 'output_data/data_out_for_ITA.txt'
    #out_file_curve   = 'output_data/data_curve_for_ITA.txt'
    #out_file_prices   = 'output_test/bond_fitting_data/test_1/prices_out_atlantia_v0.txt'
    #out_file_curve    = 'output_test/bond_fitting_data/test_1/curves_out_atlantia_v0.txt'

    out_file_prices   = 'output_test/bond_fitting_data/test_1/prices_out_mediobanca_v1.txt'
    out_file_curve    = 'output_test/bond_fitting_data/test_1/curves_out_mediobanca_v1.txt'

    #-------------------- SETUP ------------------------------------------------------------
    model_bond_0 = 'RMV'
    rf_model_fit = 'LIN'
    hr_model_fit = '2'

    #model_bond_0 = 'RFV'


    """
    flag_plot = 1
    
    flag_plot_price = False
    flag_plot_price = True
    flag_plot_price = False

    model_bond_0 = 'RMV'
    #model_bond_0 = 'RFV'
    """

    #-----------------------------------------------------------------------------
    #---------------------------  LOAD DATA  -------------------------------------
    #-----------------------------------------------------------------------------

    zc_times, zc_rf, data_zc_rf          = bf.load_curve_fromFile(input_rf, dataRef)
    zc_infl_t, zc_infl_rf, data_zc_infl  = bf.load_curve_fromFile(input_infl, dataRef)
    
    data_zc_rf                    = bf.set_model_prms(data_zc_rf)
    data_zc_infl                  = bf.set_model_prms(data_zc_infl)
    
    dictPortfolio                  = bf.loadPortfolio_fromFile_v2(portfolio_file)
    ts_dates, ts_val, data_ts_infl = bf.loadTS_fromFile(inflation_ts_file)
    
    #print 'ts_val: ', ts_val
    #print 'data_ts_infl: ', data_ts_infl
    
    #bf.FQ(999)
    
    flag_make_graph = True
    flag_plot = False
    flag_plot2 = False

    flag_dump = True
    
    #opt_dict    = fb.set_data_opt_for_fit(hr_model_fit, flag_make_graph, dataRef)
    opt_dict_rf = fb.set_data_opt_for_fit(rf_model_fit, flag_make_graph, dataRef)
    opt_dict_hr = fb.set_data_opt_for_fit(hr_model_fit, flag_make_graph, dataRef)
    
    #opt_dict['interp'] = hr_model_fit #2 = SVE modello per curva risk free e curva inflation
    
    
    rf_prms   = fb.fitting_with_plot(data_zc_rf['MatDate'], data_zc_rf['ValoreNodo'], opt_dict_rf, flag_plot)
    infl_prms = fb.fitting_with_plot(data_zc_infl['MatDate'], data_zc_infl['ValoreNodo'], opt_dict_rf, flag_plot2)
    
    data_zc_rf['Model'] = rf_model_fit 
    data_zc_rf['prms']  = rf_prms

    #data_zc_rf['Model'] = opt_dict['interp']
    #data_zc_rf['prms']  = infl_prms
    
    dict_start_params, x0, x_bnd, h_model  = bf.loadModelParams(prms_file)
    
    #-----------------------  SET OPZIONI ELABORAZIONI -------------------------------
    
 
    opt_elab = bf.settingDefaultOptions(prms_file, model_bond_0, recovery_rate, dataRef, h_model, flag_make_graph, flag_dump, out_file_prices, out_file_curve)
    
    #x0, x_bnd = set_prms_for_fit(h_model)
    
    #print 'x_bnd: ', x_bnd
    #print 'x_bnd2: ', x_bnd2
    #print 'x02: ', x02
    #print 'x0: ', x0
    
    #print data_zc_rf.keys()
    ""
    
    """
    print 'data_zc_rf[Model]: ', data_zc_rf['Model']
    print 'data_zc_rf[ValoreNodo]: ', data_zc_rf['ValoreNodo']
    print 'data_zc_rf[MatDate]: ', data_zc_rf['MatDate']
    print 'data_zc_rf[prms]: ', data_zc_rf['prms']
    
    print 'opt_elab: ',opt_elab
    print 'dictPortfolio: ', dictPortfolio
    """
    
    
    #print 'data_zc_infl: ', data_zc_infl.keys()
    #print 'data_ts_infl: ', data_ts_infl.keys()


    data_zc_rf['DiscountFactors'] = []
    
    isinList =  dictPortfolio.keys()
    
    is0 = isinList[0]
    #print dictPortfolio[is0].keys()
    #print dictPortfolio[is0]['cash flow']
    
    
    #['end date', 'cash flow', 'index rate', 'day count', 'emission date', 'clean price', 'coupon', 'coupon times', 'fixed rate', 'bond type', 'BDay', 'ytm', 'weights', 'inflRatio', 'emission price', 'isin', 'freq', 'tenor rate', 'coupon dates']

    
    #print 'opt_elab[MakeGraph]: ', opt_elab['MakeGraph']
    #----------------------- SET PARAMETRI -------------------------------------------

    #bf.compute_bond_fitting(opt_elab, dictPortfolio, zc_times, zc_rf, zc_infl_t, zc_infl_rf, ts_dates, ts_val, x0, x_bnd, dict_start_params)
    
    #print 'x0: ', x0
    #print 'x_bnd: ', x_bnd
    
    flag_res, res_out = bf.compute_bond_fitting(opt_elab, dictPortfolio, data_zc_rf, data_zc_infl, data_ts_infl, x0, x_bnd)
    
    
    
    #print 'flag_res: ', flag_res
    #print 'res_out: ', res_out
    
    """
    print 'res_out: ', res_out['opt_clean_prices']
    print 'res[dict_opt_prms]: ', res_out['dict_opt_prms']
    print 'res[survProbCum]: ', res_out['survProbCum']
    print 'res[outputDates]: ', res_out['outputDates']
    """
     
    t_end = time.time()
    
    
    t_elapse = t_end - t_start
    
    print 'Tempo di elaborazione: ', t_elapse
    
    
    import sys
    sys.exit()

    