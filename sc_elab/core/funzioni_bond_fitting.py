
import time

from scipy import optimize
from scipy.optimize import minimize
from scipy.optimize import fmin

import matplotlib.pyplot as plt
from scipy import interpolate
from os import sys


from sc_elab.core.mdates import holidays
from sc_elab.core.mdates import daycount
from sc_elab.core.mdates import busdayrule
from sc_elab.core.mdates import dateutils


#from Tkinter import *
import Tkinter as Tk

import tkMessageBox

import datetime
import dateutil.relativedelta
from datetime import datetime as dtime

from dateutil.relativedelta import *
import numpy as np
from math import exp, log



def find_indx_bf(list_ref, value):

		
		indx = 0
		n_list = len(list_ref)
		last_indx = int(n_list - 1)
		
		if (float(value) <= float(list_ref[0])):

			indx = 0

		elif (float(value) >= float(list_ref[last_indx])):
		
			indx = last_indx
		else:
		
			for i in range(1, n_list):
			
				if (float(list_ref[i]) == float(value)):
					
					indx = i
					break
				elif (float(list_ref[i]) > float(value)):

					indx = i - 1
					break
					
		return indx


def interp_lin(x_ref, y_ref, x_target):


	ln = len(x_ref)

	if (float(x_target) >= float(x_ref[ln - 1])):
	
		y_new = y_ref[ln-1]

	elif (float(x_target) <= float(x_ref[0])):

		y_new = y_ref[0]
		
	else:

		indx_m = find_indx_bf(x_ref, x_target)
		indx_p = indx_m + 1	

		x_m = float(x_ref[indx_m])
		x_p = float(x_ref[indx_p])

		"""
		print 'y_ref: ', y_ref
		print '------------------------------------'		
		print 'x_ref: ', x_ref
		print '------------------------------------'
		print 'x_target: ', x_target
		print '------------------------------------'
		print 'indx_p: ', indx_p
		print '------------------------------------'
		print 'indx_m: ', indx_m
		print '------------------------------------'
		print '------------------------------------'
		"""
		
		y_p = float(y_ref[indx_p])
		y_m = float(y_ref[indx_m])

		x_target = float(x_target)

		dx = (x_p - x_m)
		dy = (y_p - y_m)
		
		y_new = y_m + dy/dx*(x_target - x_m)
		
		
	return y_new

def from_date_to_ordinal(date_dates):


	serial_dates = []

	for i in range(0, len(date_dates)):
		dateTmp = date_dates[i].toordinal()
		serial_dates.append(dateTmp)

	return  serial_dates  

def interp_exp(x_ref, y_ref, x_target):


	ln = len(x_ref)

	if (float(x_target) >= float(x_ref[ln - 1])):
	
		y_new = y_ref[ln-1]

	elif (float(x_target) <= float(x_ref[0])):

		y_new = y_ref[0]
	else:
	
		#print 'x_target: ', x_target
		#print 'x_ref: ', x_ref
		#print '---------------------------------'
		#print 'y_ref: ', y_ref
		#print 'y_ref: ', y_ref
		#-----------------------------------------'
		
		indx_m = find_indx_bf(x_ref, x_target)
		indx_p = indx_m + 1
		
		x_m = float(x_ref[indx_m])
		x_p = float(x_ref[indx_p])
		
		y_p = float(y_ref[indx_p])
		y_m = float(y_ref[indx_m])

		ln_y_p = log(y_p + 0.0000000001)
		ln_y_m = log(y_m + 0.0000000001)

		x_target = float(x_target)

		dx = (x_p - x_m)
		ln_dy = (ln_y_p - ln_y_m)
		
		ln_y_new = ln_y_m + ln_dy/dx*(x_target - x_m)
		
		y_new = exp(ln_y_new)
		
		
	return y_new



def print_curve_onFile(py_times, s_curve, output_file):


	fout = open(output_file, "w")

	for i in range(0, len(py_times)):

		fout.write(str(py_times[i]))
		fout.write("\t")
		fout.write(str(s_curve[i]))
		fout.write("\n")
		
		fout.close()

	

def dump_vec(dt, t_end, t_vec, sp_vec, pd_vec, zrf_vec, psum_vec, price, output_file):


	fout = open(output_file, "w")
	
	fout.write('\n\n')
	fout.write('Scadenza TITOLO: ')
	fout.write(str(t_end))
	fout.write('\n')

	fout.write('Dt utilizzato: ')
	fout.write(str(dt))
	fout.write('\n')
	fout.write('\n')
	fout.write('----------------------------------------------------------\n')
	fout.write('\n')
	fout.write('Times \t\t\t S(t) \t\t\t Pd \t\t\t  zc_rf \t\t\t   Price_sum\n\n')


	for i in range(0, len(t_vec)):

		fout.write(str(t_vec[i]))
		fout.write("\t\t")
		
		fout.write(str(sp_vec[i]))
		fout.write("\t\t")

		fout.write(str(pd_vec[i]))
		fout.write("\t\t")

		fout.write(str(zrf_vec[i]))
		fout.write("\t\t")

		fout.write(str(psum_vec[i]))
		fout.write("\n")
		
	fout.write("---------------------\n")
	fout.write("Final Price: ")
	fout.write(str(price))
		
	fout.close()



def dump_vec_rmv(h_model, model_params, date_vec, t_vec, dt_vec, rf_vec, zrf_vec, sp_vec, cf_vec, cf_rf_vec, cf_ry_vec, clean, dirty, rateo, output_file):

	fout = open(output_file, "w")

	fout.write('Modello: \n')
	fout.write(str(h_model))

	fout.write('Params: \n')
	fout.write(str(model_params))
	fout.write('\n')

	fout.write('Dates\t Times\t Dt\t  rf_rate\t   z_rf\t Sp\t cf\t  Sp^(LGD))C_i\t C_ry\n ')

	for i in range(0, len(t_vec)):
	
		fout.write(str(date_vec[i]))
		fout.write("\t")

		fout.write(str(t_vec[i]))
		fout.write("\t")
		
		fout.write(str(dt_vec[i]))
		fout.write("\t")

		fout.write(str(rf_vec[i]))
		fout.write("\t")

		fout.write(str(zrf_vec[i]))
		fout.write("\t")

		fout.write(str(sp_vec[i]))
		fout.write("\t")

		fout.write(str(cf_vec[i]))
		fout.write("\t")

		fout.write(str(cf_rf_vec[i]))
		fout.write("\t")

		fout.write(str(cf_ry_vec[i]))
		fout.write("\n")


		
	fout.write("---------------------\n")

	fout.write("Clean price: ")
	fout.write(str(clean))
	fout.write("\n")

	fout.write("Dirty price: ")
	fout.write(str(dirty))
	fout.write("\n")

	fout.write("Rateo: ")
	fout.write(str(rateo))
		
	fout.close()


def dump_vec_rfv(h_model, model_params, date_vec, t_vec, dt_vec, rf_vec, zrf_vec, sp_vec, cf_vec, cf_s_rf_vec, cf_def_vec, m_def_vec, clean, dirty, rateo, output_file):

	fout = open(output_file, "w")
	
	fout.write('Modello: \n')
	fout.write(str(h_model))

	fout.write('Params: \n')
	fout.write(str(model_params))
	fout.write('\n')

	fout.write('Dates\t Times\t Dt\t rf_rate\t z_rf\t Sp\t cf\t S*Z_rf*C_i\t dm*Z_rf*c_i*RR\tD m*Z_rf*1*RR \n ')


	
	for i in range(0, len(t_vec)):
	
		fout.write(str(date_vec[i]))
		fout.write("\t")

		fout.write(str(t_vec[i]))
		fout.write("\t")
		
		fout.write(str(dt_vec[i]))
		fout.write("\t")

		fout.write(str(rf_vec[i]))
		fout.write("\t")

		fout.write(str(zrf_vec[i]))
		fout.write("\t")

		fout.write(str(sp_vec[i]))
		fout.write("\t")

		fout.write(str(cf_vec[i]))
		fout.write("\t")

		fout.write(str(cf_s_rf_vec[i]))
		fout.write("\t")

		fout.write(str(cf_def_vec[i]))
		fout.write("\t")

		fout.write(str(m_def_vec[i]))
		fout.write("\n")


	fout.write("---------------------\n")

	fout.write("Clean price: ")
	fout.write(str(clean))
	fout.write("\n")

	fout.write("Dirty price: ")
	fout.write(str(dirty))
	fout.write("\n")

	fout.write("Rateo: ")
	fout.write(str(rateo))
		
	fout.close()




def dump_portfolio_out(ISIN_vec, t_end, price_vec, output_file, tipo_modello):


	fout = open(output_file, "w")
	
	fout.write('\n')
	fout.write('Tipo modello: %s' %tipo_modello)

	fout.write('\n')
	fout.write('----------------------------------------------------------\n')
	fout.write('\n')
	fout.write('ISIN \t\t\t Scadenza \t\t\t Price \n\n')


	for i in range(0, len(t_end)):

		fout.write(str(ISIN_vec[i]))
		fout.write("\t\t")
		
		fout.write(str(t_end[i]))
		fout.write("\t\t")

		fout.write(str(price_vec[i]))
		fout.write("\n")

		
	fout.write("-----------------------------------------------\n")
		
	fout.close()



def load_curve_fromFile(inputFile, refDate):



	fin = open(inputFile, 'r')
	listInput = fin.readlines()

	n_lines = len(listInput)
			
	
	py_values = []
	py_times = []
	py_dates = []
	data_out = {}

	
	for i in range(0, n_lines):
		
		line_splitted = listInput[i].split("\t")
		
		
		timeTmp = float(line_splitted[0])
		
		timeTmp_mnth = timeTmp*12.0
		timeTmp_mnth = int(timeTmp_mnth)
		
		if (timeTmp_mnth < 1.0):

			timeTmp_week = timeTmp*52.0
			timeTmp_week = int(timeTmp_week)
			dateTmp = refDate + relativedelta(weeks=timeTmp_week)

		else:
			
			dateTmp = refDate + relativedelta(months=timeTmp_mnth)
		
		#print 'line_splitted: ', line_splitted
		dateTmp = dateTmp.date()

		py_valTmp   = float(line_splitted[1])

		py_dates.append(dateTmp)
		py_times.append(timeTmp)
		py_values.append(py_valTmp)


	data_out['MatDate'] = py_dates
	data_out['ValoreNodo'] = py_values
	
	return py_times, py_values, data_out
	
def set_model_prms(data_raw_bench):
	
	data_raw_bench['Model'] = 2
	
	dict_params = {}
	
	dict_params['const1'] = [1.0]
	dict_params['const2'] = [1.0]
	dict_params['beta0']  = [0.03]
	dict_params['beta1']  = [0.0]
	dict_params['beta2']  = [0.0]
	dict_params['beta3']  = [0.0]
	
	
	data_raw_bench['prms'] = dict_params
	
	
	return data_raw_bench
	
def loadModelParams(inputFile):



	fin = open(inputFile, 'r')
	listInput = fin.readlines()

	n_lines = len(listInput)

	modelTmp = listInput[0].split("=")
	modelTmp = str(modelTmp[1].strip())
	
	dict_params = {}
	dict_limit_params = {}
	
	x0 = []
	#x0_min = []
	#x0_max = []
	x_bnd  = []
	
	
	
	
	for i in range(1, n_lines):
		
		prmsLineTmp  = listInput[i].split("=")
		
		#print  prmsLineTmp
		prmsValueTmp = float(prmsLineTmp[1].strip())
		prmsNameTmp  = str(prmsLineTmp[0].strip())
		
		
		if (prmsNameTmp[-3:] != 'min') and  (prmsNameTmp[-3:] != 'max') :

			dict_params[prmsNameTmp]    = prmsValueTmp
			x0.append(prmsValueTmp)
		else:
			dict_limit_params[prmsNameTmp]    = prmsValueTmp
		
		
	if (modelTmp == 'NS'):

		tau = dict_params['tau']; tau_min = dict_limit_params['tau_min']; tau_max = dict_limit_params['tau_max']

		b0 = dict_params['b0'];	b0_min = dict_limit_params['b0_min']; b0_max = dict_limit_params['b0_max']
		b1 = dict_params['b1'];	b1_min = dict_limit_params['b1_min']; b1_max = dict_limit_params['b1_max']
		b2 = dict_params['b2'];	b2_min = dict_limit_params['b2_min']; b2_max = dict_limit_params['b2_max']
				
		x0      = [tau, b0, b1, b2]	
		x_bnd   = [[tau_min, tau_max], [b0_min, b0_max], [b1_min, b1_max], [b2_min,b2_max]]

	elif (modelTmp == 'SVE'):

		tau1 = dict_params['tau1']; tau1_min = dict_limit_params['tau1_min']; tau1_max = dict_limit_params['tau1_max']
		tau2 = dict_params['tau2']; tau2_min = dict_limit_params['tau2_min']; tau2_max = dict_limit_params['tau2_max']

		b0 = dict_params['b0'];	b0_min = dict_limit_params['b0_min']; b0_max = dict_limit_params['b0_max']
		b1 = dict_params['b1'];	b1_min = dict_limit_params['b1_min']; b1_max = dict_limit_params['b1_max']
		b2 = dict_params['b2'];	b2_min = dict_limit_params['b2_min']; b2_max = dict_limit_params['b2_max']
		b3 = dict_params['b3'];	b3_min = dict_limit_params['b3_min']; b3_max = dict_limit_params['b3_max']
				
		x0      = [tau1, tau2, b0, b1, b2, b3]	
		x_bnd   = [[tau1_min, tau1_max], [tau2_min, tau2_max], [b0_min, b0_max], [b1_min, b1_max], [b2_min,b2_max], [b3_min,b3_max]]

		
		
	elif (modelTmp == 'CIR'):

		r0    = dict_params['r0'];        r0_min = dict_limit_params['r0_min'];        r0_max = dict_limit_params['r0_max'];
		kappa = dict_params['kappa']; kappa_min = dict_limit_params['kappa_min'];  kappa_max = dict_limit_params['kappa_max'];
		theta = dict_params['theta']; theta_min = dict_limit_params['theta_min']; theta_max = dict_limit_params['theta_max'];	
		sigma = dict_params['sigma']; sigma_min = dict_limit_params['sigma_min']; sigma_max = dict_limit_params['sigma_max'];

		x0      = [r0, kappa, theta, sigma]	
		x_bnd   = [[r0_min, r0_max], [kappa_min, kappa_max], [theta_min, theta_max], [sigma_min,sigma_max]]
		
	else:
	
		print 'modello non contrplato!!!'

		

	return dict_params, x0, x_bnd, modelTmp

def convertDate2Time(date_vec):
	
	time_vec = []
	ln 		 = len(date_vec)

	date_0 = date_vec[0]
	
	for i in range(0, ln):
		
		timeTmp = ((date_vec[i] - date_0).days)/365.2425
		time_vec.append(timeTmp)
	
	
	return time_vec 
	
	


	



def loadPortfolio_fromFile(inputFile):



	fin = open(inputFile, 'r')
	listInput = fin.readlines()

	n_lines = len(listInput)
			
	
	dictPortfolio = {}
	#py_T     = []
	#py_c     = []
	
	
	for i in range(1, n_lines):
		
		line_splitted = listInput[i].split("\t")
		#print 'line_splitted: ', line_splitted
		
		isin_Tmp    = (line_splitted[0])
		data_endTmp0 = (line_splitted[1])
		tipo_bondTmp = (line_splitted[2])
		freqTmp      = int(line_splitted[3])
		dayCountTmp  = (line_splitted[4])
		bdaysTmp     = (line_splitted[5])
		couponTmp    = float(line_splitted[6])/100.0
		fixedRateTmp    = float(line_splitted[7])/100.0

		priceTmp     = float(line_splitted[8])
		
		data_endTmp0 = data_endTmp0.split('/')
		data_endTmp  = datetime.datetime(int(data_endTmp0[2]),int(data_endTmp0[1]),int(data_endTmp0[0]))

		
		
		dictPortfolio[isin_Tmp] = {}
		
		dictPortfolio[isin_Tmp]['bond type']= tipo_bondTmp
		dictPortfolio[isin_Tmp]['day count']= dayCountTmp
		dictPortfolio[isin_Tmp]['end date'] = data_endTmp
		dictPortfolio[isin_Tmp]['coupon']   = couponTmp
		dictPortfolio[isin_Tmp]['isin']     = isin_Tmp
		dictPortfolio[isin_Tmp]['freq']     = freqTmp
		dictPortfolio[isin_Tmp]['BDay']     = bdaysTmp
		dictPortfolio[isin_Tmp]['clean price']      = priceTmp
		dictPortfolio[isin_Tmp]['fixed rate']       = fixedRateTmp
		dictPortfolio[isin_Tmp]['coupon dates']     = []
		dictPortfolio[isin_Tmp]['coupon times']     = []
		dictPortfolio[isin_Tmp]['cash flow']        = []
		

	
	return dictPortfolio


def loadTS_fromFile(inputFile):



	fin = open(inputFile, 'r')
	listInput = fin.readlines()

	n_lines = len(listInput)
			
	
	dates_ts_list  = []
	values_ts_list = []
	
	data_ts_raw = {} 
	
	
	for i in range(1, n_lines):
		
		line_splitted = listInput[i].split("\t")
		
		date_Tmp    = line_splitted[0]
		value_Tmp   = float(line_splitted[1])


		date_Tmp 	  = date_Tmp.split('/')
		dateDate_Tmp  = datetime.datetime(int(date_Tmp[2]),int(date_Tmp[1]),int(date_Tmp[0]))

		dates_ts_list.append(dateDate_Tmp)
		values_ts_list.append(value_Tmp)
		
	data_ts_raw['MatDate'] = dates_ts_list
	data_ts_raw['Values']  = values_ts_list
	
	return dates_ts_list, values_ts_list, data_ts_raw



def fromXLSToBondFittingPortfolio(dict_start):
	
	dict_out = {}
	id_list = dict_start.keys()
	
	for i in range(0, len(id_list)):
		
		id_Tmp = id_list[i]
		isin_Tmp = dict_start[id_Tmp]['Isin']


		if (isin_Tmp != None):

			dict_out[isin_Tmp] = {}


			dict_out[isin_Tmp]['isin']  	= isin_Tmp 			
			dict_out[isin_Tmp]['emission price']  	= dict_start[id_Tmp]['Prezzo emissione'] 
			dict_out[isin_Tmp]['emission date'] 	= dict_start[id_Tmp]['Data emissione']  
			dict_out[isin_Tmp]['inflRatio']			= float(dict_start[id_Tmp]['Prezzo rimborso / Inflation Ratio']/100.0) 
	
			
			dict_out[isin_Tmp]['end date'] 			= dict_start[id_Tmp]['Data scadenza']
			dict_out[isin_Tmp]['day count']			= dict_start[id_Tmp]['Basis']	
			dict_out[isin_Tmp]['BDay'] 				= (dict_start[id_Tmp]['Adjustment']).lower()
			
			dict_out[isin_Tmp]['freq'] 				= int(dict_start[id_Tmp]['Periodicita cedola (mesi)'])
			dict_out[isin_Tmp]['weights']  			= dict_start[id_Tmp]['Peso']	
			dict_out[isin_Tmp]['tenor rate']    	= dict_start[id_Tmp]['Tenor del tasso floater (anni)']	
			dict_out[isin_Tmp]['fixed rate'] 		= float(dict_start[id_Tmp]['Tasso cedolare annuo (Fisso/spread)']/100.0)

			dict_out[isin_Tmp]['bond type']			= dict_start[id_Tmp]['Tipo tasso'] 	
			
			#if (dict_out[isin_Tmp]['bond type'] == 'FLOATING' or dict_out[isin_Tmp]['bond type'] == 'FLOATER') and dict_start[id_Tmp]['Cedola in corso'] == None:
			#	# significa che ho intercettato un errore!
				
			if (dict_out[isin_Tmp]['emission date'] == None):

				root = Tk.Tk()
				root.withdraw()
				
				msg0 = "Data emissione non valorizzata correttamente!!" 
				tkMessageBox.showinfo("Attenzione!!", msg0)
		
				root.destroy()
				return
			
			try:
				dict_out[isin_Tmp]['coupon']   			= float(dict_start[id_Tmp]['Cedola in corso']/100.0)
			
			except:

				root = Tk.Tk()
				root.withdraw()
				
				msg0 = "Cedola in corso non valorizzata correttamente!!" 
				tkMessageBox.showinfo("Attenzione!!", msg0)
		
				root.destroy()
				return
				
			dict_out[isin_Tmp]['clean price']  		= dict_start[id_Tmp]['Prezzo-MID']
			dict_out[isin_Tmp]['ytm'] 				= dict_start[id_Tmp]['YTM/DM (MID)']	
			dict_out[isin_Tmp]['index rate']      	= dict_start[id_Tmp]['Indicizzazione']
	
	
			#dict_start[isin_Tmp]['Tasso di riferiemnto']
			#dict_start[isin_Tmp]['Tasso repo']	
			#dict_start[isin_Tmp]['Giorni di fixing']	
			#dict_start[isin_Tmp]['Tipo fixing']	
			#dict_start[id_Tmp]['Tipo rimborso'] 
	
			
			dict_out[isin_Tmp]['coupon dates']     = []
			dict_out[isin_Tmp]['coupon times']     = []
			dict_out[isin_Tmp]['cash flow']        = []
	
		else:
			
			continue
	
	return dict_out


def loadPortfolio_fromFile_v2(inputFile):



	fin = open(inputFile, 'r')
	listInput = fin.readlines()

	n_lines = len(listInput)
			
	
	dictPortfolio = {}
	
	
	for i in range(1, n_lines):
		
		line_splitted = listInput[i].split("\t")
		
		isin_Tmp    = line_splitted[0]

		try:
			e_price_Tmp    		= float(line_splitted[6])
		except:
			e_price_Tmp    		= 0.0
		
		inflRatio_Tmp    	= float(line_splitted[7])/(100.0)
		
		emissionDate_Tmp0	= line_splitted[8]

		
		try:
			tenor_rate_Tmp  = float(line_splitted[16])
		except:
			tenor_rate_Tmp  = float(line_splitted[15])
		
		try:
			ytm_Tmp    		= float(line_splitted[22])
		except:
			ytm_Tmp    		= 0.0
			
		
		
		#print 'line_splitted: ', line_splitted
		
		#FQ(999)
		weight_Tmp    	= float(line_splitted[29])
		indx_rate_Tmp 	= str(line_splitted[30])

		data_endTmp0 	= line_splitted[9]
		tipo_bondTmp 	= line_splitted[3]
		freqTmp      	= int(line_splitted[15])
		dayCountTmp  	= line_splitted[13]
		bdaysTmp     	= line_splitted[14]

		try: 
			couponTmp    	= float(line_splitted[17])/100.0
			fixedRateTmp    = float(line_splitted[18])/100.0
		except: 
			couponTmp    	= 0.0
			fixedRateTmp    = 0.0

		priceTmp     	= float(line_splitted[19])

		try:
			emissionDate_Tmp0 	= emissionDate_Tmp0.split('/')
			emissionDate_Tmp  	= datetime.datetime(int(emissionDate_Tmp0[2]),int(emissionDate_Tmp0[1]),int(emissionDate_Tmp0[0]))
		except:
			emissionDate_Tmp  	= datetime.datetime(int(2011),int(11),int(11))
			
		
		data_endTmp0 	= data_endTmp0.split('/')
		data_endTmp  	= datetime.datetime(int(data_endTmp0[2]),int(data_endTmp0[1]),int(data_endTmp0[0]))


		indx_rate_Tmp = indx_rate_Tmp.split('\n')[0]

		
		dictPortfolio[isin_Tmp] = {}
		
		dictPortfolio[isin_Tmp]['bond type']	= tipo_bondTmp
		dictPortfolio[isin_Tmp]['day count']	= dayCountTmp
		dictPortfolio[isin_Tmp]['end date'] 	= data_endTmp
		dictPortfolio[isin_Tmp]['emission date'] 	= emissionDate_Tmp
		dictPortfolio[isin_Tmp]['coupon']   	= couponTmp
		dictPortfolio[isin_Tmp]['isin']     	= isin_Tmp
		dictPortfolio[isin_Tmp]['freq']     	= freqTmp
		dictPortfolio[isin_Tmp]['BDay']     	= bdaysTmp
		dictPortfolio[isin_Tmp]['clean price']  = priceTmp
		dictPortfolio[isin_Tmp]['fixed rate']   = fixedRateTmp		
		dictPortfolio[isin_Tmp]['emission price'] = e_price_Tmp
		dictPortfolio[isin_Tmp]['inflRatio'] 	  = inflRatio_Tmp
		dictPortfolio[isin_Tmp]['tenor rate']     = tenor_rate_Tmp
		dictPortfolio[isin_Tmp]['ytm']     		  = ytm_Tmp
		dictPortfolio[isin_Tmp]['weights']     	  = weight_Tmp
		dictPortfolio[isin_Tmp]['index rate']      = indx_rate_Tmp
		
		dictPortfolio[isin_Tmp]['coupon dates']     = []
		dictPortfolio[isin_Tmp]['coupon times']     = []
		dictPortfolio[isin_Tmp]['cash flow']        = []




	
	return dictPortfolio






def compute_swap_curve_from_curve(time_zc_rate, value_zc_rate, freq, time_ref):

	t_out_vec  = []
	sw_out_vec = []
	
	for i in range(0, len(time_ref)):

		t_dd = time_ref[i]
		t_sw,  sw_out = compute_swap_rate_c(time_zc_rate, value_zc_rate, freq, t_dd)

		t_out_vec.append(t_sw)
		sw_out_vec.append(sw_out)

	return t_out_vec, sw_out_vec


def compute_swap_spread_curve_from_model(model_params, model, freq, time_ref, RR):

	t_out_vec  = []
	sw_out_vec = []
	
	#RR = 0.0


	for i in range(0, len(time_ref)):

		t_dd = time_ref[i]
		t_sw,  sw_out = compute_swap_rate_m(model_params, model, freq, t_dd)


		t_out_vec.append(t_sw)
		sw_out_vec.append(sw_out*(1.0 - RR))

	return t_out_vec, sw_out_vec


def compute_swap_rate_c(time_zc_rate, value_zc_rate, freq, T):

	dt_ref   = 1.0/freq 
	len_rate = int(T/dt_ref) 
	if (len_rate <= 1):
		len_rate = 2
	
	t_sw     = dt_ref*len_rate

	z_i_sum = 0.0

	for i in range(1, len_rate):

		t_i   = i*dt_ref
		rf_i  = interp_lin(time_zc_rate, value_zc_rate, t_i)
		z_i   = (1.0 + rf_i)**(-t_i)

		z_i_sum = z_i + z_i_sum
		
	
	z_end = z_i
	
	sw_out = (1.0 - z_end)/(z_i_sum*dt_ref)

	return t_sw,  sw_out

def fwd_rate(model, prms, t1, t2):

	z1 = z_model(model, prms, t1)
	z2 = z_model(model, prms, t2)
	
	fwd_out = -1.0/(t2 - t1)*log(z2/z1)

	return fwd_out 

def compute_swap_rate_m(model_params, model, freq, T):

	dt_ref   = 1.0/freq 
	len_rate = int(T/dt_ref) 
	t_sw     = dt_ref*len_rate
	
	z_i_sum = 0.0


	for i in range(1, len_rate):
	
		t_i = i*dt_ref

		if (model == 'NS'):
			z_i = z_ns(model_params, t_i)

		elif(model == 'SVE'):		
			z_i = z_sve(model_params, t_i)
		elif(model == 'CIR'):		
			z_i = z_cir(model_params, t_i)
		else:
			print 'Modello di curva non disponibile!!'
			raise Exception

		z_i_sum = z_i + z_i_sum
		
	z_end = z_i
	
	sw_out = (1.0 - z_end)/(z_i_sum*dt_ref)

	return t_sw,  sw_out




def z_model(model, parameters, t_ref):
	
	if (model == 'CIR'):
		z_out = z_cir(parameters, t_ref)

	elif (model == 'NS'):	
		z_out = z_ns(parameters, t_ref)

	elif (model == 'SVE'):	
		z_out = z_sve(parameters, t_ref)
	
	else:
		
		print 'Model %s non definito !' %(model)
		FQ(999)

	return z_out

def z_ns(parameters, T):

	

	tau = float(parameters['tau'])
	b0  = float(parameters['b0'])
	b1  = float(parameters['b1'])
	b2  = float(parameters['b2'])

	
	dummy = ((b0*T) + tau*(b1+b2)*(1 - exp(-T/(tau))) - b2*T* exp(-T/tau))
	
	z_ns_out = exp(-dummy)
	
	
	return z_ns_out

def z_sve(parameters, T):

	
	c1 = float(parameters['tau1'])
	c2 = float(parameters['tau2'])
	b0 = float(parameters['b0'])
	b1 = float(parameters['b1'])
	b2 = float(parameters['b2'])
	b3 = float(parameters['b3'])


	T = float(T)
	
	z_sve_out =exp(-((b0*T + b1*c1*(1.0-exp(-T/c1)) + b2*(c1-(T + c1)*exp(-T/c1)) + b3*(c2-(T + c2)*exp(-T/c2)))))
	
	return z_sve_out


	
def z_cir(parameters, T):

	
	#print 'T: ', T
	
	#print 'CCC'

	kappa = float(parameters['kappa']) 
	theta = float(parameters['theta']) 
	sigma = float(parameters['sigma']) 
	r0    = float(parameters['r0']) 
	
	g       = (kappa**2 + 2*sigma**2)**0.5
	
	alpha2  = (exp(g*(T)) - 1)

	num     = 2*g*exp((g + kappa)/2*(T))
	den     = (g + kappa)*alpha2 + 2*g

	B = (2*alpha2)/den;
	A = 2*kappa*theta/(sigma*sigma)*log(num/den)

	#-------------------------------------------------
	#-------------------------------------------------



	tmp = B*r0 - A
	tmp = max(-100, tmp)
	tmp = min(+100, tmp)

	
	z = exp (-tmp)
	
	#print 'z: ', z
	
	return z



def FQ(ref):

	print '----------FIN QUI (%s) TUTTO OK-------'%(ref)
	sys.exit()
	
	
def fromDf2Rates(df_dates, df_values):
	
	date_ref = df_dates[0]
	n_dates = len(df_dates)
	ref_daycount = 'ACT/ACT'

	rate_values = []
	time_values = []
	
	for i in range(1, n_dates):
		
		
		date_i = df_dates[i]	
		df_i   = df_values[i]	

		time_i = daycount.yearfrac(date_ref, date_i, ref_daycount)
		rate_i = -np.log(df_i)/time_i
		
		rate_values.append(rate_i)
		time_values.append(np.float64(time_i))

	#print 'rate_values: ', rate_values
	#print '99999999999999999999999999'

	time_values = np.array(time_values)
	rate_values = np.array(rate_values)
	
	
	time_values = np.insert(rate_values, 0, rate_values[0])
	rate_values = np.insert(time_values, 0, 0.0)
	
	
	return time_values, rate_values


def	computeBondCalendar(date_ref, date_end, frequency, dayCount, busDay, mkt_ref, bond_type):

	
	#try:


	ref_calendar = holidays.get_calendar(mkt_ref)

	#ref_daycount = 'ACT/ACT'
	#ref_daycount = 'ACTUAL/ACTUAL ISDA'
	ref_daycount_0 = 'ACT/365 FIXED'

	ref_daycount  = dayCount		
	ref_busscnv  = busDay

	bond_sc_dates = []
	bond_sc_times = []
	bond_sc_dt    = []
	#bond_sc_dt_n  = []

	dateTmp = date_end

	date_end_a  = busdayrule.rolldate(date_end, ref_calendar, ref_busscnv)
	time_end_a  = daycount.yearfrac(date_ref, dateTmp, ref_daycount_0)


	bond_sc_dates.append(date_end_a)
	bond_sc_times.append(time_end_a)

	k = 1

	dateTmp_o = date_end_a
	dateTmp_o_0 = date_end


	if (bond_type == 'ZC'):
		
		bond_sc_dates.append(date_ref)
		bond_sc_times.append(0.0)
		bond_sc_dt.append(time_end_a)
		bond_sc_dt.append(0.0)

		bond_sc_dates.reverse()
		bond_sc_times.reverse()
		bond_sc_dt.reverse()

		
	else:
		

		while (dateTmp_o >= date_ref):

			dateTmp_n_0 = date_end - relativedelta(months=k*frequency)

			dateTmp_n = busdayrule.rolldate(dateTmp_n_0, ref_calendar, ref_busscnv)
			timeTmp_n = daycount.yearfrac(date_ref, dateTmp_n, ref_daycount_0)
			dtTmp     = daycount.yearfrac(dateTmp_n_0, dateTmp_o_0, ref_daycount)


			bond_sc_dt.append(dtTmp)
			bond_sc_dates.append(dateTmp_n)
			bond_sc_times.append(timeTmp_n)

			dateTmp_o = dateTmp_n
			dateTmp_o_0 = dateTmp_n_0

			#timeTmp_o = timeTmp_n

			k = k + 1


		bond_sc_dt.append(dtTmp)

		bond_sc_dates.reverse()
		bond_sc_times.reverse()
		bond_sc_dt.reverse()

	return bond_sc_dates, bond_sc_times, bond_sc_dt

	#except:
			
	#	strE = 'computeBondCalendar'
	#	print strE

	

	
def computeBondYPrice(Y, data_portfolio, date_ref):

	try:

		#date_end_val   = data_portfolio['end date']    
		#freq_val       = data_portfolio['freq']       
		#dayCount_val   = data_portfolio['day count']   
		#bDay_val       = data_portfolio['BDay']       
		c_val          = data_portfolio['coupon']
		#bond_dates     = data_portfolio['coupon dates']
		coupon_times   = data_portfolio['coupon times']
		coupon_values  = data_portfolio['cash flow']
		bType  		   = data_portfolio['bond type']  

		#----------------------------------------------------------------------------------
		#------------------------------ START COMPUTATION ---------------------------------
		#----------------------------------------------------------------------------------


		c_ref      = c_val
		ln_b = len(coupon_times)


		price_tmp = 0.0


		price_sum = 0.0
		price_tmp = 0.0

		rateo =  (-coupon_times[0])*c_ref
		rateo = coupon_values[0]

		for i in range(1, ln_b):


			t_i   = coupon_times[i]		
			c_i   = coupon_values[i]

			z_i_rf = (1.0 + Y)**(-t_i)

			price_tmp = c_i*z_i_rf
			price_sum  = price_tmp + price_sum

			#-----------------------------------------------------------------------------
			#------------------------- CALCOLO RIMBORSO CAPITALE -------------------------
			#-----------------------------------------------------------------------------

		t_end = t_i


		z_rf_end  = (1.0 + Y)**(-t_end)
		priceTmp = z_rf_end*(1.0)
		price = price_sum +priceTmp

		dirty_price = 100*price
		#clean_price = 100*(price - rateo)
		
		if (bType == 'ZC'):
			dirty_price = dirty_price + rateo 			
		
		
		clean_price = (dirty_price - rateo)
		


		return 	clean_price, dirty_price 	

	except:
			
		strErr = 'computeBondYPrice'
		print strErr


def fromDataModelToYTM(zc_times, zc_rf, h_model, model_prms, LGD, t_i):
	
		if (h_model == 'CIR'):
			sp_p  = z_cir(model_prms, float(t_i))	 	
		elif (h_model == 'NS'):
			sp_p  = z_ns(model_prms, float(t_i))	 
		else:
			sp_p  = z_sve(model_prms, float(t_i))	 

		rf_i   = interp_lin(zc_times, zc_rf, float(t_i))
		z_i_rf = (1.0 + rf_i)**(-t_i)
		z_i_risky = z_i_rf*(sp_p**(LGD))
		
		y_from_model = z_i_risky**(-1.0/t_i)-1

		
		return y_from_model
	
	
	

def computeBondPriceFromCF(model_params, data_portfolio, opt_elab, zc_times, zc_rf, zc_infl_times, zc_infl, ts_infl_dates, ts_infl_values):





	#try:

	RR             = opt_elab['RR'] 
	tipo_modello   = opt_elab['BondModel'] 
	h_modello      = opt_elab['HRateModel'] 
	date_ref       = opt_elab['DataRef']
	
	

	ISIN_val       = data_portfolio['isin']
	#date_end_val   = data_portfolio['end date']    
	tipo_bond_val  = data_portfolio['bond type']  
	#freq_val       = data_portfolio['freq']       
	#dayCount_val   = data_portfolio['day count']   
	#bDay_val       = data_portfolio['BDay']       
	#c_val          = data_portfolio['coupon']
	coupon_dates   = data_portfolio['coupon dates']
	coupon_times   = data_portfolio['coupon times']
	coupon_dt      = data_portfolio['coupon dt']
	cf_val         = data_portfolio['cash flow']
	indx_security  = data_portfolio['index rate'] 

	startDate      = data_portfolio['emission date']  
	inflRatio_anag = data_portfolio['inflRatio'] 
	
	LGD       = 1.0 - RR
	
	#print 'ISIN_val: ', ISIN_val
	#print 'CTimes: ', coupon_times
	

	#----------------------------------------------------------------------------------
	#------------------------------ START COMPUTATION ---------------------------------
	#----------------------------------------------------------------------------------
	
	
	#if (ISIN_val == 'XS1689739347'):

	#	print 'AAAA'
	#	print 'model_params: ', model_params

	#ln_titoli = len(ISIN_val)

	#c_ref      = c_val
	#date_end   = date_end_val
	#frequency  = freq_val
	#dayCount   = dayCount_val
	#busDay     = bDay_val
	#s_old   = 1.0
	#P_i     = 1.0
	#price_tmp = 0.0


	ln_b = len(coupon_times)

	s_curve = []
	s_curve.append(1.0)


	sp_m = 1.0

	price_sum = 0.0
	price_tmp = 0.0

	rateo = cf_val[0]
	
	dump_flag = 0
	
	if (dump_flag == 1):

		date_vec     = []
		t_vec        = []
		dt_vec       = []
		rf_vec       = []
		zrf_vec      = []
		sp_vec       = []
		cf_vec       = []
		cf_rf_vec    = []
		cf_ry_vec    = []
		cf_def_vec   = []
		m_def_vec    = []
		cf_ry_vec    = []
		cf_s_rf_vec  = []


	
	#FQ(9988)
	for i in range(1, ln_b):


		t_i   = coupon_times[i]		
		dt    = coupon_dt[i]
		c_i   = cf_val[i] 

		
		if (h_modello == 'CIR'):
			
			#print 't_i: ', t_i

			sp_p  = z_cir(model_params, float(t_i))	 	
		elif (h_modello == 'NS'):
			sp_p  = z_ns(model_params, float(t_i))	 
		else:
			sp_p  = z_sve(model_params, float(t_i))	 

		rf_i   = interp_lin(zc_times, zc_rf, float(t_i))
		
		
		z_i_rf = (1.0 + rf_i)**(-t_i)

		#-----------------------------------------------------------------------------

		if (indx_security == 'CPTFEMU'):
			
			inflRatio_anag_n = float(inflRatio_anag/1.0)
			inflRatio = inflRatio_anag_n
			
			"""
			inflRatio_ts = inflationRatio(date_ref, startDate, ts_infl_dates, ts_infl_values)
			chk_inflRatio = np.abs(inflRatio_anag_n - inflRatio_ts)/inflRatio_anag_n
			
			if chk_inflRatio < 0.9:
				
				inflRatio = inflRatio_anag_n
			else:
	
				print 'infaltion ratio ricalcolato non in linea con quello presente nell anagrafica del titolo diff (IR_anag - IR_ts)/IR_anag >0.1'
				print 'inflRatio da anagrafica %s' %inflRatio_anag_n
				print 'inflRatio da ts %s' %inflRatio_ts
				
				raise Exception
	
				#str = 'infaltion ratio ricalcolato non in linea con quello presente nell anagrafica del titolo diff (IR_ts - IR_anag)/IR_anag >0.1'
				#print str
			"""
			
			#print 'zc_infl_times: ', zc_infl_times
			#print 'zc_infl: ', zc_infl
			#print '======================================'
	
			infl_rate = interp_lin(zc_infl_times, zc_infl, float(t_i))
			z_infl 	  = exp(-infl_rate*float(t_i))
			
			infl_r_i = (1.0/z_infl)*inflRatio			
	
		else:
			
			infl_r_i  = 1.0
			inflRatio = 1.0



		"""
		if (indx_security == 'CPTFEMU'):

			infl_rate = interp_lin(zc_infl_times, zc_infl, float(t_i))
			z_infl 	  = exp(-infl_rate*float(t_i))
			infl_r_i = 1.0/z_infl*inflRatio_anag
		else:
			infl_r_i = 1.0
		"""
		
		if (tipo_modello == 'RFV'):

			ds        = (sp_m - sp_p)

			price_tmp = sp_p*c_i*infl_r_i + RR*ds*(1.0 + c_i)*infl_r_i
			price_tmp = price_tmp*z_i_rf

		else:		 	
			ds        = 1.0
			z_i_risky = z_i_rf*(sp_p**(LGD))
			
			price_tmp = z_i_risky*c_i*infl_r_i

		price_sum  = price_tmp + price_sum
		

		"""
		print 't_i: ', t_i
		print 'z_i_rf: ', z_i_rf
		#print 'z_infl: ', z_infl
		print 'z_i_risky: ', z_i_risky
		print 'c_i: ', c_i
		print 'price_tmp: ', price_tmp		
		print 'price_sum: ', price_sum
		print '---------------------------------'
		"""

		sp_m = sp_p

		if (dump_flag == 1):

			date_vec.append(coupon_dates[i])
			t_vec.append(t_i)
			dt_vec.append(dt)
			rf_vec.append(rf_i)
			zrf_vec.append(z_i_rf)
			sp_vec.append(sp_p)
			cf_vec.append(c_i)
			
			if (tipo_modello == 'RMV'):

				cf_rf_vec.append((sp_p**(LGD))*c_i*infl_r_i)
				cf_ry_vec.append(price_tmp)

			else:
			
				cf_s_rf_vec.append(z_i_rf*sp_p*c_i*infl_r_i)
				cf_def_vec.append(z_i_rf*RR*ds*c_i*infl_r_i)
				m_def_vec.append(z_i_rf*RR*ds*1.0)
				cf_ry_vec.append(price_tmp)

			
			#-----------------------------------------------------------------------------

	t_end = t_i

	if (h_modello == 'CIR'):
		s_end  = z_cir(model_params, float(t_end))
	elif (h_modello == 'NS'):
		s_end  = z_ns(model_params, float(t_end))
	else:
		s_end  = z_sve(model_params, float(t_end))

	
	if (indx_security == 'CPTFEMU'):
		
		inflRatio = inflRatio_anag_n
		

		"""
		inflRatio_ts = inflationRatio(date_ref, startDate, ts_infl_dates, ts_infl_values)
		
		chk_inflRatio = np.abs(inflRatio_anag_n - inflRatio_ts)/inflRatio_anag
		
		if chk_inflRatio < 0.9:
			
			inflRatio = inflRatio_anag_n
			
		else:
		
			print 'infaltion ratio ricalcolato non in linea con quello presente nell anagrafica del titolo diff (IR_anag - IR_ts)/IR_anag >0.1'
			print 'inflRatio da anagrafica %s' %inflRatio_anag_n
			print 'inflRatio da ts %s' %inflRatio_ts
			
			raise Exception
		"""
		
		infl_rate = interp_lin(zc_infl_times, zc_infl, float(t_i))
		z_infl 	  = exp(-infl_rate*float(t_i))

		
		infl_r_i = (1.0/z_infl)*inflRatio

	else:
		
		infl_r_i  = 1.0
		inflRatio = 1.0


	rf_end = interp_lin(zc_times, zc_rf, t_end)
	z_rf_end  = (1.0 + rf_end)**(-t_end)

	if (tipo_modello == 'RFV'):

		zRisky   = s_end**(LGD)
		priceTmp = z_rf_end*s_end*(1.0)*infl_r_i
		price    = price_sum + priceTmp
	else: 

		zRisky   = z_rf_end*s_end**(LGD)
		

		priceTmp = zRisky*(1.0)*infl_r_i
		price = price_sum + priceTmp

	#pd_end = sp_m - s_end

	"""
	print 't_i: ', t_i
	print 'c_i: ', c_i
	print 'z_i_rf: ', z_i_rf
	print 'z_infl: ', z_infl
	print 'z_i_risky: ', z_i_risky
	print 'price_tmp: ', price_tmp		
	print 'price: ', price
	print '---------------------------------'
	"""
	
	dirty_price = 100*price
		
	if (tipo_bond_val == 'ZC'):

		dirty_price = dirty_price + rateo
	
	clean_price = dirty_price - rateo
	clean_price = clean_price/inflRatio
	
	#clean_price = 100*(price - rateo)
	
	#print 
	
	
	#dump_flag = 0
	

	if (dump_flag == 1):
	
		date_vec.append(coupon_dates[i])
		t_vec.append(t_i)
		dt_vec.append(dt)
		rf_vec.append(rf_end)
		zrf_vec.append(z_rf_end)
		sp_vec.append(s_end)
		cf_vec.append(c_i)


		if (tipo_modello == 'RMV'):

			cf_rf_vec.append(zRisky)
			cf_ry_vec.append(priceTmp)

			#output_file = 'dump/' + ISIN_val + '_' + h_modello + '_RMV_XX.txt'
			#output_file = 'output_test/bond_fitting_data/dump/' + ISIN_val + '_' + h_modello + '_RMV_XX.txt'
			output_file = 'output_test/bond_fitting_data/dump/' + ISIN_val + '_' + h_modello + '_RMV_XX.txt'

			
			dump_vec_rmv(h_modello, model_params, date_vec, t_vec, dt_vec, rf_vec, 
				    zrf_vec, sp_vec, cf_vec, cf_rf_vec,
				    cf_ry_vec, clean_price, dirty_price, rateo, output_file)

		else:	

			cf_s_rf_vec.append(z_rf_end*s_end*(1.0))
			cf_def_vec.append(0.0)
			m_def_vec.append(0.0)
			cf_ry_vec.append(price_tmp)

			output_file = 'dump/' + ISIN_val + '_' + h_modello + '_RFV.txt'

			dump_vec_rfv(h_modello, model_params, date_vec, t_vec, dt_vec, rf_vec, 
				    zrf_vec, sp_vec, cf_vec, cf_s_rf_vec, 
				    cf_def_vec, m_def_vec, clean_price, dirty_price, rateo, output_file)
	
	
	
	
	return 	clean_price, dirty_price 	

	#except:
			
	#	str = 'computeBondPriceFromCF'
	#	print str



#def blm_clean_price = ytm2blmPrice()

"""
def ytm2blmPrice(coupon_values, dt_vec, rateo, YTM):

	de        = (nextCdate - setDate)/deltaRef;
	c_tmp_new = 100*(CRate*h)/(1 + YTM)^de;
	df_tmp    = 1.0/(1 + YTM)^de;


	for i in range(0, coupon_values):

		df_tmp    = df_tmp/(1 + YTM)^dt_vec[i];
		c_tmp     = 100*(CRate*dt_vec[i])*df_tmp;
		c_tmp_new = c_tmp_new + c_tmp;


	c_end_m = 100*df_tmp;   
	blmPrice_secco   = (c_tmp_new + c_end_m) - rateo;

	return blmPrice_secco
"""



def computeBondCF(data_portfolio, zc_times, zc_rf, rf_prms, rf_model, mkt_ref, date_ref, zc_infl_t, zc_infl_val, infl_prms, infl_model):



	#try:


	#ISIN_val       = data_portfolio['isin']
	date_end_val   = data_portfolio['end date']    
	tipo_bond_val  = data_portfolio['bond type']  
	freq_val       = data_portfolio['freq']       
	dayCount_val   = data_portfolio['day count']   
	bDay_val       = data_portfolio['BDay']       
	c_val          = data_portfolio['coupon']


	#----------------------------------------------------------------------------------
	#------------------------------ START COMPUTATION ---------------------------------
	#----------------------------------------------------------------------------------

	c_ref      = c_val
	date_end   = date_end_val
	frequency  = freq_val
	dayCount   = dayCount_val
	busDay     = bDay_val

	coupon_values = []


	if (tipo_bond_val == 'ZC'):
		frequency = 12*100		
		coupon_dates, coupon_times, coupon_dt  = computeBondCalendar(date_ref, date_end, frequency, dayCount, busDay, mkt_ref, tipo_bond_val)
		coupon_dates[0] = date_ref

	else:

		coupon_dates, coupon_times, coupon_dt  = computeBondCalendar(date_ref, date_end, frequency, dayCount, busDay, mkt_ref, tipo_bond_val)


	ln_b = len(coupon_times)

	if (tipo_bond_val == 'FIXED'):
		rateo =  (-coupon_times[0])*c_ref*100.0

	elif(tipo_bond_val == 'FLOATER'):
		rateo =  (-coupon_times[0])*data_portfolio['fixed rate']*100.0  
	else:
		
		emissionPrice = data_portfolio['emission price']
		emiDate 	  = data_portfolio['emission date']
		
		#print 'date_end_val: ', date_end_val
		#print 'emiDate: ', emiDate
		
		time_zc_gg = (date_end_val - emiDate).days
		t_rateo_gg = (date_ref - emiDate).days

		time_zc = float(time_zc_gg)/365.2425
		t_rateo = float(t_rateo_gg)/365.2425
		rateo   = ((100.0 - emissionPrice)/time_zc)*t_rateo

	coupon_values.append(rateo)
	#print 'rateo: ', rateo	

	for i in range(1, ln_b):


		#t_i   = coupon_times[i]		
		t_im  = coupon_times[i - 1]
		dt    = coupon_dt[i]

		if (tipo_bond_val == 'FIXED'):

			c_i   = float(c_ref*dt)
		elif(tipo_bond_val == 'FLOATER'):

			if (i == 1):

				fixed_rate = data_portfolio['fixed rate'] 
				c_i        = float(fixed_rate*dt)

			else:
				
				dt_tenor = data_portfolio['tenor rate']
				dt_tenor = float(dt_tenor/12.0)
				
				t_i_n    = t_im + dt_tenor
				
				rf_im  	 = interp_lin(zc_times, zc_rf, float(t_im))
				rf_i   	 = interp_lin(zc_times, zc_rf, float(t_i_n))

				z_im   	 = exp(-t_im*rf_im)
				z_i    	 = exp(-t_i_n*rf_i)

				fwd 	 = -log(z_i/z_im)/dt_tenor

				spread 	 = c_ref
				#fwd     = (z_im/z_i - 1.0)/dt_tenor
				c_i      = float(fwd + spread)*dt

		else:
				c_i    = 0.0

		coupon_values.append(c_i)
		
	#print 'coupon_values: ', coupon_values


	return 	coupon_dates, coupon_times, coupon_values, coupon_dt  	

#except:

#	str = 'computeBondCF'
#	print str





	#except:
	
	#	strE = 'loss_bf_cf'
	#	print strE






"""
def computePYmodel_rate(dict_opt_params, time_ref, LGD, zc_times, zc_rf, model_type):

	
	try:
		dict_ry = {}

		for k in time_ref:

			t_ref     = dict_bond_times[k]

			if(model_type == 'CIR'):
				sp        = z_cir(dict_opt_params, t_ref)
			elif(model_type == 'NS'):
				sp        = z_ns(dict_opt_params, t_ref)
			else:
				sp        = z_sve(dict_opt_params, t_ref)

			zspread    = -log(sp**(LGD))/t_ref

			rf_i      = interp_lin(zc_times, zc_rf, t_ref)
			z_ry      = rf_i +  zspread

			dict_ry[k]= z_ry

		return dict_ry
	except:
	
		str = 'computePYmodel_rate'
		print str
"""


def set_prms_for_fit(model_fit):

	if (model_fit == 'SVE'): #->SVE

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


def settingDefaultOptions(prms_file, model_in, RR, date_ref, h_model, flag_make_graph, flag_dump, out_file_prices, out_file_curve):


	modello_h_rate = h_model
	model_bond  = model_in
	opt_elab = {}
	opt_elab['ElabFromYtm']  = False 

	opt_elab['RR']          = RR 
	opt_elab['BondModel']   = model_bond 
	opt_elab['HRateModel']  = modello_h_rate
	opt_elab['DataRef']     = date_ref
	opt_elab['MKTRef']      = 'de'
	opt_elab['MakeGraph']   = flag_make_graph
	opt_elab['MakeDump']    = flag_dump
	opt_elab['out_file_prices']    = out_file_prices
	opt_elab['out_file_curve']    = out_file_curve
	


	return opt_elab 

	
def setOptParams(ff, h_model):


	dict_opt_params = {}
	
	
	if (h_model == "CIR"):

		r0_opt    = ff.x[0];	kappa_opt = ff.x[1]
		theta_opt = ff.x[2];	sigma_opt = ff.x[3]

		dict_opt_params["r0"]    = r0_opt;	dict_opt_params["kappa"] = kappa_opt
		dict_opt_params["theta"] = theta_opt;	dict_opt_params["sigma"] = sigma_opt
	
	elif (h_model == "SVE"):

		tau1_opt= ff.x[0];	tau2_opt= ff.x[1];	b0_opt  = ff.x[2];
		b1_opt  = ff.x[3];	b2_opt  = ff.x[4];	b3_opt  = ff.x[5];

		dict_opt_params["tau1"] = tau1_opt;	dict_opt_params["tau2"] = tau2_opt;
		dict_opt_params["b0"]   = b0_opt;	dict_opt_params["b1"]   = b1_opt;
		dict_opt_params["b2"]   = b2_opt;	dict_opt_params["b3"]   = b3_opt;

	elif (h_model == "NS"):

		tau_opt = ff.x[0];	b0_opt  = ff.x[1];
		b1_opt  = ff.x[2];	b2_opt  = ff.x[3];

		dict_opt_params["tau"] = tau_opt;	dict_opt_params["b0"] = b0_opt;
		dict_opt_params["b1"] = b1_opt;		dict_opt_params["b2"] = b2_opt;
		

	return dict_opt_params


def	computeTimesDatesRef(opt_elab, dictPtf):
	
	
	
	
	date_eval = opt_elab['DataRef']
	date_0 = dtime.fromordinal(date_eval.toordinal())
	
	dateOut = []
	timeOut = []
	
	endDateRef = date_eval
	isinList = dictPtf.keys()
	
	for i in range(0, len(isinList)):
		
		isinTmp = isinList[i]
		endDate = dictPtf[isinTmp]['end date']
		
		if (endDate > endDateRef):
			endDateRef = endDate
			
			
	timeLast = 	((endDateRef - date_eval).days)/365.2425
	
	if (timeLast > 30):
		ln = 3 + int(timeLast)
	else:		
		ln = 33
		
		
	for i in range(0, ln):
		
		if (i<4):
			timeTmp = float(i*0.5)
			dateTmp = date_0 + relativedelta(months=int(12*i*0.5))
			
		else:
			timeTmp = i - 2
			dateTmp = date_0 + relativedelta(months=12*(i-2))
		#dateTmp = date_0 + relativedelta(months=12*timeTmp)
		
		dd = int(dateTmp.day)
		mm = int(dateTmp.month)
		yy = int(dateTmp.year)
		
		dateDateTmp = datetime.date(yy, mm, dd)
		
		dateOut.append(dateDateTmp)
		timeOut.append(float(timeTmp))

	return timeOut, dateOut


def	computeX2(list_opt_clean_prices, list_mkt_clean_prices):

	sum_diff = 0.0
	n_bond = len(list_opt_clean_prices)
	
	
	if (n_bond < 2):

		opt_clean_priceTmp = list_opt_clean_prices[0]
		mkt_clean_priceTmp = list_mkt_clean_prices[0]
	
		diff = (opt_clean_priceTmp - mkt_clean_priceTmp)
		sum_diff = sum_diff + diff*diff
		i = 1

	else:	
		for i in range(0, n_bond):
			
			opt_clean_priceTmp = list_opt_clean_prices[i]
			mkt_clean_priceTmp = list_mkt_clean_prices[i]
		
			diff = (opt_clean_priceTmp - mkt_clean_priceTmp)
			sum_diff = sum_diff + diff*diff
		
	x2 = sum_diff/float(i)
	
	return x2




def computePYmodel_rate_n(dict_opt_params, time_ref, freq, LGD, zc_times, zc_rf, model_type):


	
	freq = 2.0
	dt_gg = 1.0/365.2425
	dt_swp = 1.0/freq
	
	sw_ry_vec = []
	sw_sp_vec = []
	sw_rf_vec = []
	z_sp_vec  = []
	
	sp_vec = []
	h_vec  = []

	time_vec  = []
	
	z_rf_i = 1
	
	indx_ref = find_indx_bf(time_ref, dt_swp)
	
	time_ref_new = time_ref[0:indx_ref]
	n_t_ref      = len(time_ref)
	n_ref        = int(time_ref[n_t_ref-1]/dt_swp)
	
	for i in range(1, n_ref+1):
		
		
		time_ref_new.append(i*dt_swp)
	
	#for t_tmp in time_ref:
	for t_tmp in time_ref_new:
		
		if (t_tmp < 0.00001):
			t_tmp = 0.00001

		n_t = int(t_tmp/dt_swp) + 1
		
		sum_ry = 0.0
		sum_rf = 0.0
		sum_sp = 0.0
		
		h_tmp  = fwd_rate(model_type, dict_opt_params, t_tmp, t_tmp + dt_gg)
		
		if (t_tmp <= 0.00001):
			sp_tmp = 1.0
		else:
			sp_tmp = z_model(model_type, dict_opt_params, t_tmp)
			
		
		rf_tmp = interp_lin(zc_times, zc_rf, t_tmp)

		"""		
		print 't_tmp: ', t_tmp
		print 'rf_tmp: ', rf_tmp
		print 'h_tmp: ', h_tmp
		print 'sp_tmp: ', sp_tmp
		print '----------------------------'
		"""

		z_rf_tmp  = 1.0/(1.0 + rf_tmp*t_tmp)
		
		#z_rf_tmp  = exp(-rf_tmp*t_tmp)
		z_sp_tmp  = sp_tmp**(LGD)
		z_ry_tmp  = z_sp_tmp*z_rf_tmp

		z_spTmp = -log(z_ry_tmp/z_rf_tmp)/t_tmp
		
		"""
		print 't_tmp: ', t_tmp
		print 'n_t: ', n_t
		print '+++++++++++++++++++++++++'
		"""
		
		#print 'rf_tmp: ', rf_tmp
		#print '=============================='

		
		if (t_tmp > dt_swp):
		
			for i in range(0, n_t):

				t_sw_tmp     = dt_swp*(i+1)
				
				sp     = z_model(model_type, dict_opt_params, t_tmp)
				sp_i   = z_model(model_type, dict_opt_params, t_sw_tmp)

				rf_i    = interp_lin(zc_times, zc_rf, t_sw_tmp)

				z_rf_i  = 1/(1.0  + rf_i*dt_swp)**i

				z_sp_i  = sp_i**(LGD)
				z_ry_i  = z_sp_i*z_rf_i

				sum_ry = sum_ry + z_ry_i
				sum_rf = sum_rf + z_rf_i
				sum_sp = sum_sp + z_sp_i
				
				
				#print 't_tmp: ', t_tmp
				#print 'rf_i: ', rf_i
				#print '----------------------------'
		
		
		
		if (t_tmp <= dt_swp):
		
			sp   = z_model(model_type, dict_opt_params, t_tmp)
			z_sp_i  = sp**(LGD)

			rf_i    = interp_lin(zc_times, zc_rf, t_tmp)

			z_rf_i  = 1.0/(1.0  + rf_i*t_tmp)
			z_ry_i  = z_sp_i*z_rf_i
			#z_sp_i  = sp**(LGD)

			sw_ry = ((1.0/z_ry_i) - 1.0)/t_tmp
			sw_rf = ((1.0/z_rf_i) - 1.0)/t_tmp
			sw_sp = ((1.0/z_sp_i) - 1.0)/t_tmp


		else:
		
			sw_ry = (1.0 - z_ry_i)/(dt_swp*sum_ry)
			sw_rf = (1.0 - z_rf_i)/(dt_swp*sum_rf)	
			sw_sp = (1.0 - z_sp_i)/(dt_swp*sum_sp)

			"""
			print 't_tmp: ', t_tmp
			print 'sum_rf: ', sum_rf
			print 'z_rf_i: ', z_rf_i
			print 'sw_rf: ', sw_rf
			print '----------------------------'
			"""
			#FQ(999)


		z_sp_vec.append(z_spTmp)
		sp_vec.append(sp)
		h_vec.append(h_tmp)
		
		sw_ry_vec.append(sw_ry)
		sw_rf_vec.append(sw_rf)
		sw_sp_vec.append(sw_sp)
		time_vec.append(t_tmp)


	z_sp_vec[0] = z_sp_vec[1]
	sp_vec[0] 	= sp_vec[1]
	h_vec[0] 	= h_vec[1]
	
	sw_ry_vec[0] = sw_ry_vec[1]
	sw_rf_vec[0] = sw_rf_vec[1]
	sw_sp_vec[0] = sw_sp_vec[1]


	time_vec_n = np.interp(time_ref, time_ref_new, time_vec)
	z_sp_vec_n = np.interp(time_ref, time_ref_new, z_sp_vec)
	sp_vec_n = np.interp(time_ref, time_ref_new, sp_vec)
	h_vec_n = np.interp(time_ref, time_ref_new, h_vec)
	
	sw_ry_vec_n = np.interp(time_ref, time_ref_new, sw_ry_vec)
	sw_rf_vec_n = np.interp(time_ref, time_ref_new, sw_rf_vec)
	sw_sp_vec_n = np.interp(time_ref, time_ref_new, sw_sp_vec)
	
	if (time_vec_n[0] <= 0.001):
		sp_vec_n[0] = 1.0
		

	return sw_ry_vec_n, sw_rf_vec_n, sw_sp_vec_n, z_sp_vec_n, sp_vec_n, h_vec_n, time_vec_n
	#return sw_ry_vec, sw_rf_vec, sw_sp_vec, z_sp_vec, sp_vec, h_vec, time_vec

	#except:
	
	#	str = 'computePYmodel_rate'
	#	print str


def computeYTM(Y, dict_anag, refDate ):

	try:

		refPrice = dict_anag['clean price']
		clean_yprice, dirty_yprice = computeBondYPrice(Y, dict_anag, refDate)
		diff = abs (refPrice - clean_yprice)#/refPrice

		return diff

	except:
	
		strE = 'computeYTM'
		print strE




def computePortfolioYTM(dictPortfolio, refDate):


	try:
		isin_list = dictPortfolio.keys()
		ytm_out = {}
		Yguess  = 0.03 

		for isin_tmp in isin_list:

			ytm_tmp = fmin(computeYTM, Yguess, args=(dictPortfolio[isin_tmp],refDate), xtol=1e-12)

			ytm_out[isin_tmp] = ytm_tmp

		return ytm_out

	except:
	
		strErr = 'computePortfolioYTM'
		print strErr


def	inflationRatio(refDate, startDate, ts_dates, ts_values):


	INFL_LEG = 3

	IRef    = interpInflation(refDate, ts_dates, ts_values, INFL_LEG)
	IStart  = interpInflation(startDate, ts_dates, ts_values, INFL_LEG)
	
	
	#print 'IRef: ',IRef
	#print 'IStart: ',IStart

	iRatio  = IRef/IStart

	return iRatio

def seasonal_adjustment(infl_index_dates, infl_index_values, method):
	
	"""
	 Function to estimate the seasonal adjustments (annualized) given a vector 
	 of past index values according to two possible methods
	
	 INPUT:
	 infl_index_dates          Column vector of dates (increasing order)
	 infl_index_values         Column vector of index values
	 method                    (Avg = Averages of averages | Dum = Dummy regression model)
	 
	 OUTPUT:
	 additive_seasonal         Column vector (12 x 1) of seasonals starting
	                           from "month(infl_index_dates(end,1)) + 1"
	 reference_month           month of the first seasonal (e.g. reference_month = 1 means January) 
	
	 Notes: 
	 1- additive_seasonal will sum to 0
	 2- infl_index_values should be observed index values (not observed
	    "Reference Index" values for a specific product)
	 3- the Avg method potentially discards information; if there are exactly 
	    "n years + 1 month" of obs the two methods are equivalent 
	"""
	
	# Find and export reference month
	ln = len(infl_index_dates)
	m = infl_index_dates[ln-1].month
	if (m == 12):
		reference_month = 1
	else:
		reference_month = m + 1
	
	# Construct the vector of MoM inflation
	n             = len(infl_index_dates)
	MoM_inflation = log(infl_index_values[1:n-1]/ infl_index_values[0:n-2])
	n_years       = int((n - 1)/ 12)
	n_spare_month = (n - 1) - n_years*12.0
	
	flag_method = True
	
	
	ln_mom = len(MoM_inflation)
	if (flag_method == True):				
		# Discard observations of least recent spare months
		if (n_spare_month > 0):
			MoM_inflation = MoM_inflation[n_spare_month:ln_mom-1]
		
		# Detrend the MoM inflation for each year observed inflation
		
		inflation_demeaned = np.zeros(12*n_years)
		for i in range(1,n_years):
		
			targetData = MoM_inflation[12*(i-1):12*i-1]
			inflation_demeaned[12*(i-1):12*i-1] = MoM_inflation[12*(i-1):12*i-1] - np.mean(targetData)
		
		# Compute averages by month of the detrended inflation rates
		
		n_years = int(len(MoM_inflation)/12.0)
		
		additive_seasonals = np.zeros(12)
		for i in range(0, 12):
			
			sum_valTmp = 0.0
			for j in range(0, n_years):
				
				valTmp = inflation_demeaned[j*12 + i]
				sum_valTmp = sum_valTmp + valTmp
			
			sum_valTmp = sum_valTmp/n_years
			
			additive_seasonals[i] = sum_valTmp
			
	return additive_seasonals, reference_month 
	
	




def  interpInflation(date, ts_dates, ts_values, INFL_LEG):
	
	

	a = date + relativedelta(months=int(-INFL_LEG))
	b = date + relativedelta(months=int(-(INFL_LEG-1)))
	
	a_ordinal = a.toordinal()
	b_ordinal = b.toordinal()


	# recupero i valori
	
	dates_ordinal = from_date_to_ordinal(ts_dates)
	
	indx_a       = find_indx_bf(dates_ordinal, a_ordinal)
	indx_b       = find_indx_bf(dates_ordinal, b_ordinal)
	
	Ia = ts_values[indx_a]
	Ib = ts_values[indx_b]
	yy = date.year
	mm = date.month
	
	eom = dateutils.eom(yy, mm)
	eom_dd = eom.day 
	dd = date.day
	
	I0 = Ia + ( dd - 1.)/eom_dd*(Ib-Ia)
	

	return I0



def computePortfolioPricesFromCF(model_params, dictPortfolio, opt_elab, zc_times, zc_rf, rf_prms, rf_model, zc_infl_times, zc_infl, infl_prms, infl_model, ts_infl_dates, ts_infl_values):

	#try:
	isin_list = dictPortfolio.keys()

	clean_prices = {}
	dirty_prices = {}
	bond_times   = {}


	for i in range(0, len(isin_list)):


		isin_tmp = isin_list[i]

		dateEndTmp  = dictPortfolio[isin_tmp]['end date']
		
		#print 'end: ', type(dateEndTmp)
		#print 'ref: ', type(opt_elab['DataRef'])

		t_endTmp    = dateEndTmp -  opt_elab['DataRef']
		t_endTmp    = float(t_endTmp.days)/365.0

		bond_times[isin_tmp]  = t_endTmp

		clean_price, dirty_price = computeBondPriceFromCF(model_params, dictPortfolio[isin_tmp], opt_elab, zc_times, zc_rf, zc_infl_times, zc_infl, ts_infl_dates, ts_infl_values)

		clean_prices[isin_tmp] = clean_price
		dirty_prices[isin_tmp] = dirty_price


	return clean_prices, dirty_prices, bond_times  

	#except:

	#	strE = 'computePortfolioPricesFromCF'
	#	print strE




def computePortfolioCF(dictPortfolio, opt_elab, zc_times, zc_rf, rf_prms, rf_model, zc_infl_t, zc_infl_val, infl_prms, infl_model):


	#try:

	mkt_ref   = opt_elab['MKTRef']
	date_ref  = opt_elab['DataRef']
	isin_list = dictPortfolio.keys()


	for i in range(0, len(isin_list)):


		isin_tmp = isin_list[i]
		
		#print 'isin-tmp: ', isin_tmp
		coupon_dates, coupon_times, coupon_values, coupon_dt = computeBondCF(dictPortfolio[isin_tmp], zc_times, zc_rf, rf_prms, rf_model, mkt_ref, date_ref, zc_infl_t, zc_infl_val, infl_prms, infl_model)
		

		dictPortfolio[isin_tmp]['coupon dates']  = coupon_dates
		dictPortfolio[isin_tmp]['coupon times']  = coupon_times
		dictPortfolio[isin_tmp]['coupon dt']     = coupon_dt
		dictPortfolio[isin_tmp]['cash flow']     = coupon_values


	return dictPortfolio  

	
	#except:
	
	#	str = 'Error in computePortfolioCF'
	#	print str



def write_dump_out_v2(times, surv, h_rate, z_spread, sw_spread, sw_rf, out_file_curve):

	fout   = open(out_file_curve, 'w')
	
	fout.write('Times\t Survival\t HazardRate\t zSpread\t pySpread\t pyRiskFree\n')
	

	
	for i in range(0, len(times)):
	
		timeTmp 	=  times[i]
		survTmp 	=  surv[i]
		hTmp 		=  h_rate[i]
		zSpreadTmp  =  z_spread[i]
		swSpreadTmp =  sw_spread[i]
		swRFTmp 	=  sw_rf[i]
		
		fout.write('%2.5f\t %2.5f\t %2.5f\t %2.5f\t  %2.5f\t %2.5f\n'%(timeTmp, survTmp, hTmp, zSpreadTmp, swSpreadTmp, swRFTmp))
	
	fout.close()




def write_dump_out(dictPortfolio, sw_ry, sw_times, ytm_mkt, dict_bond_times, opt_clean_prices, out_file):
	
	list_k_ref = sorted(dict_bond_times, key = lambda key: dict_bond_times[key])

	list_opt_clean_prices   = []
	
	list_mkt_clean_prices   = []
	list_bond_times         = []

	list_ytm_model          = []
	list_ytm_mkt            = []
	list_isin               = []

	for k in list_k_ref:
	
		ytm_model_tmp   = interp_lin(sw_times, sw_ry, dict_bond_times[k])
		
		list_opt_clean_prices.append(opt_clean_prices[k])
		list_mkt_clean_prices.append(dictPortfolio[k]['clean price'])
		list_bond_times.append(dict_bond_times[k])
		list_isin.append(k)
		list_ytm_model.append(ytm_model_tmp)
		list_ytm_mkt.append(ytm_mkt[k])
	
	fout   = open(out_file, 'w')
	
	fout.write('ISIN\t bond_times\t clean_price\t mkt_price\t ytm_model\t ytm_mkt\n')
	

	
	for i in range(0, len(list_opt_clean_prices)):
	
		a0 =  list_isin[i]
		a1 =  list_bond_times[i]
		a2 =  list_opt_clean_prices[i]
		a3 =  list_mkt_clean_prices[i]
		a4 =  list_ytm_model[i]
		a5 =  list_ytm_mkt[i][0]
		
		fout.write('%s\t %2.5f\t %2.5f\t %2.5f\t %2.5f\t %2.5f\n'%(a0,a1,a2,a3,a4,a5))
	
	fout.close()
	

def set_var_out(dictPortfolio, zc_times, zc_rf, h_model, dict_opt_params, LGD, sw_ry, sw_times, ytm_mkt, dict_bond_times, opt_clean_prices):

	
	
	list_k_ref = sorted(dict_bond_times, key = lambda key: dict_bond_times[key])
		

	list_opt_clean_prices   = []
	list_mkt_clean_prices   = []
	list_bond_times         = []

	list_ytm_model          = []
	list_ytm_mkt            = []

	for k in list_k_ref:
	
		ytm_model_tmp   = interp_lin(sw_times, sw_ry, dict_bond_times[k])
		
		ytm_model_tmp2 = fromDataModelToYTM(zc_times, zc_rf, h_model, dict_opt_params, LGD, dict_bond_times[k])
		
		list_opt_clean_prices.append(opt_clean_prices[k])
		list_mkt_clean_prices.append(dictPortfolio[k]['clean price'])
		list_bond_times.append(dict_bond_times[k])
		list_ytm_model.append(ytm_model_tmp2)
		list_ytm_mkt.append(ytm_mkt[k])
	
	
	return list_bond_times, list_ytm_mkt, list_ytm_model, list_opt_clean_prices, list_mkt_clean_prices




def	plotResults(model_bond, time_ref, sw_rf, sw_spread, sw_ry, list_bond_times, list_ytm_mkt, list_ytm_model, list_opt_clean_prices_n, list_mkt_clean_prices_n, indx_to_plot, flag_plot_price):


	import matplotlib
	matplotlib.use('TkAgg')
	
	from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
	from matplotlib.backend_bases import key_press_handler
	
	
	from matplotlib.figure import Figure
	
	
	root = Tk.Tk()
	root.wm_title("Plot Bond fitting results")
	
	f = Figure(figsize=(5, 4), dpi=100)
	a = f.add_subplot(111)

	#flag_plot_price = False
	
	if (len(list_bond_times) < 2): 
		flag_plot_price = True
		flag_single_bond = True
	
	if(model_bond == 'RMV' and flag_plot_price != True):
	#if(model_bond == 'RMV'):

		a.plot(list_bond_times, list_ytm_mkt, 'go', label='YTM mkt')
		#a.plot(list_bond_times, list_ytm_model, 'kx', label='YTM model')
	
		a.plot(time_ref[1:indx_to_plot], sw_rf[1:indx_to_plot], '-.b', label='SW rf')
		a.plot(time_ref[1:indx_to_plot], sw_spread[1:indx_to_plot], '-.r', label = 'SW spread')	
		a.plot(time_ref[1:indx_to_plot], sw_ry[1:indx_to_plot], '-k', label = 'SW risky')
		a.set_title('Fitting clean bond prices using %s model'  %(model_bond))
		
		#a.legend(['YTM mkt', 'SW rf', 'SW spread', 'SW risky'], loc = 0)
		a.set_xlabel('Maturities [years]')
		a.set_ylabel('YTM')
		
		legend = a.legend(loc='upper left', shadow=False)
		
	else:

		a.plot(list_bond_times, list_opt_clean_prices_n, 'x', label = 'Fit prices')
		a.plot(list_bond_times, list_mkt_clean_prices_n, 'o', label = 'MKT prices')
		


		a.set_title('Fitting clean bond prices using %s model' %(model_bond) )
		a.set_xlabel('Maturities [years]')
		a.set_ylabel('Prices')
		legend = a.legend(loc='upper left', shadow=False)
		
		
		if (flag_single_bond == True): 
			a.axes.set_ylim([10,200])
	

	#plt.show()
	
	
	canvas = FigureCanvasTkAgg(f, master=root)
	
	
	canvas.show()
	
	
	canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
	
	toolbar = NavigationToolbar2TkAgg(canvas, root)
	toolbar.update()
	canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
	

	def on_key_event(event):
		print('you pressed %s' % event.key)
		key_press_handler(event, canvas, toolbar)

	canvas.mpl_connect('key_press_event', on_key_event)


	def _quit():
		root.quit()     # stops mainloop
		root.destroy()  # this is necessary on Windows to prevent
						# Fatal Python Error: PyEval_RestoreThread: NULL tstate
		
	button = Tk.Button(master=root, text='Quit', command=_quit)
	button.pack(side=Tk.BOTTOM)
	
	Tk.mainloop()



def chk_inflatio_ratio(dictPortfolio, date_ref, ts_infl_dates, ts_infl_values):
	
	list_isin = dictPortfolio.keys()

	for isinTmp in 	list_isin:

		inflRatio_anag_n = dictPortfolio[isinTmp]['inflRatio']
		startDate 		=  dictPortfolio[isinTmp]['emission date']

		
		inflRatio_ts = inflationRatio(date_ref, startDate, ts_infl_dates, ts_infl_values)
		
		chk_inflRatio = np.abs(inflRatio_anag_n - inflRatio_ts)/inflRatio_anag_n
		
		if (chk_inflRatio <0.9):
			
			
			continue
		
		else:
			
			#from Tkinter import *
			import tkMessageBox
			
			# significa che ho intercettato un errore!
			root = Tk.Tk()
			root.withdraw()
			msg0 = "Inflation ratio non coerente con i nostri dati storici sull'inflazione" 
			tkMessageBox.showinfo("Attenzione!!", msg0)
			
			root.destroy()
			return


	
	return 1



def loss_bf_cf(list_model_params, dictPortfolio, opt_elab, zc_times, zc_rf, rf_prms, rf_model, zc_infl_times, zc_infl, infl_prms, infl_model, ts_infl_dates, ts_infl_values):

	#try:
	
	#print 'i'

	dict_model_params = {}

	if (opt_elab['HRateModel'] == 'CIR'):
	
		
		r0    = list_model_params[0]
		kappa = list_model_params[1]
		theta = list_model_params[2]
		sigma = list_model_params[3]

		dict_model_params['kappa'] = kappa 
		dict_model_params['theta'] = theta
		dict_model_params['sigma'] = sigma
		dict_model_params['r0']    = r0


	elif (opt_elab['HRateModel'] == 'NS'):
	
		
		tau = list_model_params[0]
		b0  = list_model_params[1]
		b1  = list_model_params[2]
		b2  = list_model_params[3]

		dict_model_params['tau'] = tau 
		dict_model_params['b0']  = b0
		dict_model_params['b1']  = b1
		dict_model_params['b2']  = b2
		
	else:

		
		tau1 = list_model_params[0]
		tau2 = list_model_params[1]
		b0  = list_model_params[2]
		b1  = list_model_params[3]
		b2  = list_model_params[4]
		b3  = list_model_params[5]

		dict_model_params['tau1'] = tau1 
		dict_model_params['tau2'] = tau2 
		dict_model_params['b0']  = b0
		dict_model_params['b1']  = b1
		dict_model_params['b2']  = b2
		dict_model_params['b3']  = b3

	clean_prices, dirty_prices, bond_times = computePortfolioPricesFromCF(dict_model_params, dictPortfolio, opt_elab, zc_times, zc_rf, rf_prms, rf_model, zc_infl_times, zc_infl, infl_prms, infl_model, ts_infl_dates, ts_infl_values)
	diff_sum = 0.0
	k  = 0
	
	wTmpSum = 0.0

	for isin_tmp in clean_prices.keys():
		
		k = k + 1
		
		wTmp = float(dictPortfolio[isin_tmp]['weights'])
		
		
		wTmpSum = wTmpSum + wTmp

		bond_mkt 	= dictPortfolio[isin_tmp]['clean price']
		bond_mdl    = clean_prices[isin_tmp]

		diff 		= abs(bond_mdl - bond_mkt)
		#diff 		= abs(bond_mdl - bond_mkt)/bond_mkt

		diff 		= wTmp*diff*diff
		diff_sum 	= diff_sum  +  diff

	#diff_sum = float(diff_sum)/float(k)/wTmpSum
	diff_sum = float(diff_sum)/wTmpSum
	
	return diff_sum



def	compute_bond_fitting(opt_elab, dictPortfolio, data_zc_rf, data_zc_infl, data_ts_infl, x0, x_bnd):


	zc_rf_dates  = data_zc_rf['MatDate']
	zc_rf = data_zc_rf['ValoreNodo']
	zc_rf_df = data_zc_rf['DiscountFactors']
	
	rf_prms = data_zc_rf['prms']
	rf_model = data_zc_rf['Model']
	
	infl_prms = data_zc_infl['prms']
	infl_model = data_zc_infl['Model']
	
	zc_infl_dates  = data_zc_infl['MatDate']
	zc_infl_val    = data_zc_infl['ValoreNodo']

	ts_infl_dates  = data_ts_infl['MatDate']
	ts_infl_values = data_ts_infl['Values']
	
	is0 = dictPortfolio.keys()[0]
	indicizzazione = dictPortfolio[is0]['index rate']
	
	
	if (indicizzazione == 'CPTFEMU'):
		
		zc_infl_t = convertDate2Time(zc_infl_dates)
	else:
		zc_infl_t = []
		
	
	
	zc_times  = convertDate2Time(zc_rf_dates)
	

	if len(zc_rf_df) < 2:
		
		zc_rf_df = []
		for i in range(0, len(zc_times)):
			
			timeTmp = zc_times[i]
			rateTmp = zc_rf[i]
			
			dfTmp = np.exp(-timeTmp*rateTmp)
			zc_rf_df.append(dfTmp)
		
		data_zc_rf['DiscounntFactors'] = zc_rf_df


	if len(zc_rf) < 2:
		
		zc_rf = []
		for i in range(1, len(zc_times)):
			
			timeTmp = zc_times[i]
			zc_rf_dfTmp = zc_rf_df[i]
			
			rateTmp = -np.log(zc_rf_dfTmp)/timeTmp
			zc_rf.append(rateTmp)
		
		zc_rf.insert(0, zc_rf[0])

		data_zc_rf['ValoreNodo'] = zc_rf


	h_model        = opt_elab['HRateModel']
	model_bond     = opt_elab['BondModel']
	refDate        = opt_elab['DataRef']
	RR             = opt_elab['RR']
	
	LGD = 1.0 - RR
	
	method_opt = 'TNC'
	
	flag_dump = opt_elab['MakeDump']
	flag_plot = opt_elab['MakeGraph']
	
	isinList = dictPortfolio.keys()
	is0 	 = isinList[0]
	indx_bond=  dictPortfolio[is0]['index rate']
	
	if (indx_bond =='CPTFEMU'):
		
		date_ref =	opt_elab['DataRef']
		chk_inflatio_ratio(dictPortfolio, date_ref, ts_infl_dates, ts_infl_values)
	
	
	if (model_bond == 'RFV') or (indx_bond == 'CPTFEMU'):
		flag_plot_price = 1
	else:
		flag_plot_price = 0
	
	



	
	
	#-------------------------------------------------------------------------------------------
	#----------------------- Inizio elaborazione -----------------------------------------------
	#-------------------------------------------------------------------------------------------
	
	#print 'start compute portfolio CF'
	
	dictPortfolio_new =  computePortfolioCF(dictPortfolio, opt_elab, zc_times, zc_rf, rf_prms, rf_model, zc_infl_t, zc_infl_val, infl_prms, infl_model)


	#print 'start minimize'

	
	#start_clean_prices, start_dirty_prices, dict_bond_times = computePortfolioPricesFromCF(dict_start_params, dictPortfolio_new, opt_elab, zc_times, zc_rf, rf_prms, rf_model, zc_infl_t, zc_infl_val, infl_prms, infl_model, ts_infl_dates, ts_infl_values)
	
	ff  = optimize.minimize(loss_bf_cf, x0,args = (dictPortfolio_new, opt_elab, zc_times, zc_rf, rf_prms, rf_model, zc_infl_t, zc_infl_val, infl_prms, infl_model, ts_infl_dates, ts_infl_values), method = method_opt,  bounds = x_bnd)


	#-------------------------- unpack/set opt params -----------------------------------------

	dict_opt_params = setOptParams(ff, h_model)#
	#print 'dict_opt_params: ',dict_opt_params
	
	#------------------------ COMPUTE RESULTS -------------------------------------------------
	
	time_ref, dateOut = computeTimesDatesRef(opt_elab, dictPortfolio_new)
	freq = 2.0
	

	"""
	print '------- AAAAAAAAAAAAAAA------------'
	print 'dict_opt_params: ', dict_opt_params
	print 'opt_elab: ', opt_elab
	print 'zc_times: ', zc_times
	print 'zc_rf: ', zc_rf
	print 'rf_prms: ', rf_prms
	print 'rf_model: ', rf_model
	print '-------- BBBBBBBBBBBBBB-------------'
	"""
	
	
	opt_clean_prices, opt_dirty_prices, dict_bond_times     = computePortfolioPricesFromCF(dict_opt_params, dictPortfolio_new, opt_elab, zc_times, zc_rf, rf_prms, rf_model, zc_infl_t, zc_infl_val, infl_prms, infl_model, ts_infl_dates, ts_infl_values)
	
	ytm_mkt   = computePortfolioYTM(dictPortfolio, refDate)
	
	sw_ry, sw_rf, sw_spread, z_spread, surv, hr_values, sw_times = computePYmodel_rate_n(dict_opt_params, time_ref, freq, LGD, zc_times, zc_rf, h_model)
	
	
	list_bond_times, list_ytm_mkt, list_ytm_model, list_opt_clean_prices, list_mkt_clean_prices = set_var_out(dictPortfolio, zc_times, zc_rf, h_model, dict_opt_params, LGD, sw_ry, sw_times, ytm_mkt, dict_bond_times, opt_clean_prices)
	
	
	x2 = computeX2(list_opt_clean_prices, list_mkt_clean_prices)
	
	if (flag_dump == True):



		out_file_prices = opt_elab['out_file_prices']
		out_file_curve	= opt_elab['out_file_curve']

	
		write_dump_out(dictPortfolio, sw_ry, sw_times, ytm_mkt, dict_bond_times, opt_clean_prices, out_file_prices)
		write_dump_out_v2(sw_times, surv, hr_values, z_spread, sw_spread, sw_rf, out_file_curve)



	if (flag_plot == True):
		
		ln = len(list_bond_times)
		t_max = list_bond_times[ln-1]
		indx_max_to_plot = int(t_max) + 5
		
		
		
		plotResults(model_bond, time_ref, sw_rf, sw_spread, sw_ry, list_bond_times, list_ytm_mkt, list_ytm_model, list_opt_clean_prices, list_mkt_clean_prices, indx_max_to_plot, flag_plot_price)

	output_data ={}
	
	output_data['x2']    = x2

	output_data['pyRiskFree']       = sw_rf	
	output_data['zcSpread']         = z_spread
	output_data['pySpread']         = sw_spread
	output_data['hazardRate']       = hr_values
	output_data['survProbCum']      = surv
	output_data['outputDates']      = dateOut
	
	output_data['ytm_mkt']          = list_ytm_mkt
	output_data['ytm_model']        = list_ytm_model
	output_data['bondTimes']        = list_bond_times
	output_data['opt_clean_prices'] = list_opt_clean_prices
	output_data['mkt_clean_prices'] = list_mkt_clean_prices
	output_data['dict_opt_prms'] 	= dict_opt_params
	
	
	#--------------------------------------------------------------------------------------------------------------
	
	#elapse = time.time() - t_start
	#print 'Tempo di calcolo: ', elapse
	
	#---------------------------------------------------------------------------------------------------------------
	flag_res = 1
	
	return flag_res,  output_data

