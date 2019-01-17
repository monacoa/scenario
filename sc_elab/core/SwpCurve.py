
import sys
import os
from sc_elab.excel_hook.connection import *
import sc_elab.core.mdates.busdayrule as busD
import sc_elab.core.mdates.holidays as holy
import sc_elab.core.mdates.dateutils as dateutils
import datetime
from dateutil.relativedelta import relativedelta
import sc_elab.core.funzioni_base as fb
from DEF_core import dict_segm2,dict_segm
import numpy as np


def revDict(do):
    dn = {}
    for k in do.keys():
        dn[do[k]] = k
    return dn

def computeMaturityString(self, cn, c, type, code):
    if type == "Fut":
        return code
    qry = '''
        SELECT DESCRIZIONE FROM ZTerm WHERE CODICE_TERM = '%s'
    ''' % (code)

    cn.execute(qry)
    res = (cn.fetchall())[0][0]
    return res



class Segm:
    def __init__(self):

        self.name   = ""
        # --- dizionario dell'anagrafica del segmento
        self.anag   = {}

        self.tags   = []
        self.mats   = []
        self.dates  = []
        self.values = []
        self.usage  = []


    def show(self):
        print "ANAG:"
        print self.anag
        print "....."

        print "TAGS  :", self.tags
        print "MATS  :", self.mats
        print "dates :", self.dates
        print "values:", self.values
        print "usage :", self.usage

    def getDayCount(self):
        try : return self.anag['DAY COUNT CONV.']
        except: zzzzzzzzzzzz

    def getAdj(self):
        try: return self.anag['DATE ADJ.']
        except: zzzzzzzzzzzz

    def getCapitalization(self):
        if(self.name == "Dep") or (self.name == "Libor") or (self.name =="Fut"):
            return self.anag['REGIME CAPIT.']
        else: return ""


class Curve(object):
    def __init__(self):
        object.__init__(self)
        self.setDefaults()

    def setDefaults(self):
        # --------
        # anagrafica
        # --------
        self.description    = ""
        self.curr           = ""
        self.ref_date       = ""
        self.source         = "Bloomberg"
        self.quotation      = "MID"
        self.download_type  = "0"
        self.emittente      = '999'
        self.rating         = '999'
        self.settore        = '999'
        self.seniority      = '999'
        self.type           = 'Swap'
        self.floater_tenor  = ''
        self.cal = ''
        # ---------
        self.HWparms = {}
        # ----------
        # segmenti (dict of classes
        # ----------
        self.segms = {}



    def getCurveCode(self):
        codeL   = "C"
        codeL  += self.curr
        codeR   = "ZC"
        codeR  += "BLM0" if self.source=='Bloomberg'  else "OTH"
        dd      = "%02d"%self.ref_date.day
        mm      = "%02d"%self.ref_date.month
        yy      = (str(self.ref_date.year))
        codeR  += dd+mm+yy[-2:]
        return codeL, codeR

    def getStrSegms(self):
        sep = ""
        res = ""
        for k in dict_segm2.keys():
            try:
                tmp = self.segms[dict_segm2[k]]
                res += sep + k
                sep = ","
            except:
                pass

        return res


    def show(self):
        print "self.description"    , self.description
        print "self.curr"           , self.curr
        print "self.ref_date"       , self.ref_date
        print "self.type"           , self.type
        print "self.source"         , self.source
        print "self.quotation"      , self.quotation
        print "self.rating"         , self.rating
        print "self.settore"        , self.settore
        print "self.seniority"      , self.seniority
        print "print self.type"     , self.type
        print "self.floater_tenor"  , self.floater_tenor
        print "self.cal"            , self.cal
        print "self.download_type"  , self.download_type
        print "self.emittente"      , self.emittente

        for k in self.HWparms.keys():
            print "HWparms [",k,"]: ", self.HWparms[k]

        print "===== Begin segms ======"
        for k in self.segms.keys():
            print "------------>Segmento:", k
            self.segms[k].show()
            print "<------------Fine segmento"
        print "==== End  segms ===="

    def computeTags(self):
        con = Connection()
        cn = con.db_anag()
        for k in self.segms.keys():
            seg = self.segms[k]
            for mat,val in  zip(seg.mats, seg.values):

                seg.tags.append(computeMaturityString(self,cn, self, k, mat))
                val = val / 100.0

        cn.close()

    def fillAnagSegm(self):

        con = Connection()
        cn = con.db_anag()

        for k in self.segms.keys():

            seg = self.segms[k]
            seg.name = k
            if (k == 'GSwp1M') or (k == 'GSwp3M'):
                ts = 'LIBOR'

            else:
                ts = ((revDict(dict_segm)[k][1:])+" libor") if ((revDict(dict_segm)[k][1:]) == "Future") else  (revDict(dict_segm)[k][1:])

            qry = '''
                SELECT NOME_ATTRIBUTO, VALORE_ATTRIBUTO FROM MKT_Segmentazione_D_N WHERE CODICE_SEGMENTAZIONE IN
                (
                    SELECT CODICE_SEGMENTAZIONE FROM MKT_Curve
                    WHERE  TIPO_CURVA = 'SWAP'
                    AND    VALUTA     = '%s'
                    AND TIPO_SEGMENTO = '%s'
                )
            ''' % (self.curr, ts)

            cn.execute(qry)
            res = cn.fetchall()

            seg.anag = {}
            for record in res:
                name  = record[0]
                value = record[1]
                seg.anag[name]= value
        cn.close()


    def addWorkingDays(self, ds, days, adj):
        calendar = holy.get_calendar(self.cal)
        tmp = int(days)
        while tmp > 0:

            ds = ds + datetime.timedelta(days=1)
            ds = busD.rolldate_from_db(ds, calendar, adj)
            tmp = tmp - 1
        return ds



    def computeDates(self):
        #----
        #recupero il calendario di riferimento
        #----

        for s in self.segms.keys():
            #[(u'DATE ADJ.', u'Following'), (u'DAY COUNT CONV.', u'30/360'), (u'FREQ. PAYMENT', u'6M'), (u'LAG', u'2')]
            adj       = self.segms[s].anag['DATE ADJ.']
            day_count = self.segms[s].anag['DAY COUNT CONV.']
            fix_days  = self.segms[s].anag['LAG']


            #per i futures non faccio nulla
            if s[0] == 'F':
                i = 0
                for do in self.segms[s].mats:
                    i = i+1
                    print 'do: ', do
                    ds = self.addWorkingDays(do, fix_days, adj)
                    self.segms[s].dates.append(ds)
                    self.segms[s].tags.append(str(i))
            else:
                # case s[0] == 'D','L','S', 'G1', 'G3'
                rdate = dateutils.asdatetime(self.ref_date)
                for mat,tag in zip (self.segms[s].mats, self.segms[s].tags):
                    if   tag == "O/N" : date = self.addWorkingDays(rdate , 1, adj)
                    elif tag == "T/N" : date = self.addWorkingDays(rdate, 2, adj)
                    elif tag[-1]=="D" :
                        nday = int(tag[:-1])
                        date = self.addWorkingDays(rdate, nday, adj)
                    elif tag[-1] == "W":
                        date = self.addWorkingDays(rdate, fix_days, adj)
                        date = date + datetime.timedelta(weeks=(int(tag[:-1])))
                        calendar = holy.get_calendar(self.cal)
                        date = busD.rolldate_from_db(date, calendar, adj)
                    elif tag[-1] == "M":
                        calendar = holy.get_calendar(self.cal)
                        date = self.addWorkingDays(rdate, fix_days, adj)
                        date = date + relativedelta(months=(int(tag[:-1])))
                        date = busD.rolldate_from_db(date, calendar, adj)
                    elif tag[-1] == "Y":
                        calendar = holy.get_calendar(self.cal)
                        date = self.addWorkingDays(rdate, fix_days, adj)
                        date = date + relativedelta(years=(int(tag[:-1])))
                        date = busD.rolldate_from_db(date, calendar, adj)
                    else: qqqqqqqqqqqqqqqqqqqqqqqq
                    self.segms[s].dates.append(date)



    def init_finalize(self):
        print "SONO IN INIT FINALIZE!!!"
        #compute output strings
        self.computeTags()
        #compute anag segms
        self.fillAnagSegm()
        #compute output dates
        self.computeDates()




    def loadDataFromDB(self):
        con = Connection()

        c_a = con.db_data()

        #----------------------------------------
        # SELEZIONO IL GIUSTO TIPO E DESCRIZIONE
        #----------------------------------------


        qry = "SELECT distinct BloombergTicker FROM DProCurve where (TipoDato = 'CLibor' or  TipoDato = 'CSwap' or TipoDato = 'CDepositi' or TipoDato = 'CDeopositi') and  (Descrizione = '%s')"%self.description




        c_a.execute(qry)
        res = c_a.fetchall()

        blm_tckrs_lst_str_n = ""
        sep = ""
        for record in res:
            blm_tckrs_lst_str_n += sep + "'" + record[0] + "'"
            sep = ","


        blm_tckrs_lst_str = blm_tckrs_lst_str_n 


        #----------------------------------------
        # SELEZIONO IL GIUSTO CONTRIBUTOR ALLA DATA PREFISSATA
        #----------------------------------------

        c_d = con.db_data()


        qry = '''
             SELECT TipoDato, Contributor, BloombergTicker FROM DProCurve where BloombergTicker IN
             (
                SELECT DISTINCT BLOOMBERGTICKER FROM DProTS_master WHERE
                BLOOMBERGTICKER IN (%s) AND DATA = '%s'
                )
            ''' % (blm_tckrs_lst_str, str(self.ref_date).replace("-", ""))

        print 'QRA: ', qry
        c_d.execute(qry)
        res_n = c_d.fetchall()

        res_len = len(res_n)

        c_list   = []
        t_list   = []
        blm_list = []

        for i in range(0, res_len):
            c_list.append(res_n[i][1])
            t_list.append(res_n[i][0])
            blm_list.append(res_n[i][2])

        
        blm_tkr_list_str_o = ""
        sep = ""
        for record in res_n:
            blm_tkr_list_str_o += sep + "'" + record[2] + "'"
            sep = ","

        
        
        type_ctb_dict  =  fb.set_type_ctb_dict(c_list, t_list)                   
        
        

        list_type = type_ctb_dict.keys()
        n_type    = len(list_type)
        blm_tkr_list_str_n = "" 

        for i in range(0, n_type):
            
            typeTmp = list_type[i]

            ctb_list_tmp = type_ctb_dict[typeTmp]
            
            if (len(ctb_list_tmp) < 2):
                
                ctbTmp = ctb_list_tmp[0]
                
            else:
                
                if  ('THOMSON_REUTERS' in c_list):
                    
                    ctbTmp = 'THOMSON_REUTERS'
                    
                else: 
                    
                    ctbTmp = c_list[0]
                    
        

            
            qry = '''
                 (
                    SELECT DISTINCT BLOOMBERGTICKER FROM DProTS_master WHERE
                    BLOOMBERGTICKER IN (%s) AND Contributor = '%s' AND  TipoDato = '%s'
                    )
                ''' % (blm_tkr_list_str_o, ctbTmp, typeTmp)
            

        
            print 'qry xxx: ', qry
            
            c_d.execute(qry)
            res_a = c_d.fetchall()

            #blm_tkr_list_str_n = ""
            sep = ""
            for record in res_a:
                blm_tkr_list_str_n += sep + "'" + record[0] + "'"
                sep = ","

            

        #print 'blm_tkr_list_str_n: ', blm_tkr_list_str_n
        blm_tckrs_lst_str = blm_tkr_list_str_n

        # ---------
        # RECUPERO CURRENCY E NUMERO/TIPOLOGIA DI SEGMENTI
        # ---------

        c_d = con.db_data()

        qry = '''
             SELECT DISTINCT currency, TipoDato FROM DProCurve where BloombergTicker IN
             (
                SELECT DISTINCT BLOOMBERGTICKER FROM DProTS_master WHERE
                BLOOMBERGTICKER IN (%s) AND
                 DATA = '%s'
                )
            ''' % (blm_tkr_list_str_n, str(self.ref_date).replace("-", ""))

        c_d.execute(qry)
        res = c_d.fetchall()
        self.curr = res[0][0]



        
        if self.curr == "EUR":
            self.floater_tenor = "6M"
            self.cal = "de.eurex"
        elif self.curr == "USD":
            self.floater_tenor = "3M"
            self.cal = 'us'
        elif self.curr == 'GBP':
            self.cal = 'uk'
        elif self.curr=='CAD':
            self.cal = 'ca'
        else:
            self.cal = 'us'
            self.floater_tenor = "3M"
        #----

        if self.quotation == "MID":
            quotazione = "valoremid"
        elif self.quotation == "ASK":
            quotazione = "valoreask"
        elif self.quotation == "BID":
            quotazione = "valorebid"
        else:
            mmmmmmmmmmmmmmmm

        for record in res:
            segm = record[1]
            qry = ""
            if ((segm == "CDepositi") or (segm == "CSwap") or (segm == 'CLibor')):

                qry = '''
                        SELECT      DProCurve.maturityInt, DProTS_master.%s
                        FROM        DProCurve
                        INNER JOIN  DProTS_master
                        ON          DProCurve.BloombergTicker = DProTS_master.BloombergTicker
                        WHERE       DProTS_master.data='%s'
                        AND         DProCurve.TipoDato = '%s'
                        AND         DProCurve.BloombergTicker in (%s)
                        ORDER BY    DProCurve.maturityInt
                        ''' % (quotazione, str(self.ref_date).replace("-", ""), segm, blm_tckrs_lst_str)
            elif segm == "CFuture":
                qry = '''
                        SELECT      DProTS_master.ScadenzaFuture, DProTS_master.%s
                        FROM        DProCurve
                        INNER JOIN  DProTS_master
                        ON          DProCurve.BloombergTicker = DProTS_master.BloombergTicker
                        WHERE       DProTS_master.data= '%s'
                        AND         DProCurve.TipoDato ='%s'
                        AND         DProCurve.BloombergTicker in (%s)
                        ORDER BY    DProTS_master.ScadenzaFuture
                        ''' % (quotazione, str(self.ref_date).replace("-", ""), segm, blm_tckrs_lst_str)

            else:
                #CSwap1M, CSwap3M, CSwap6M
                #print "segmento:", segm
                self.floater_tenor = segm[-2:]
                #print "floater tenor:", self.floater_tenor

                qry = '''
                        SELECT      DProCurve.maturityInt, DProTS_master.%s
                        FROM        DProCurve
                        INNER JOIN  DProTS_master
                        ON          DProCurve.BloombergTicker = DProTS_master.BloombergTicker
                        WHERE       DProTS_master.data='%s'
                        AND         DProCurve.TipoDato = '%s'
                        AND         DProCurve.BloombergTicker in (%s)
                        ORDER BY    DProCurve.maturityInt
                        ''' % (quotazione, str(self.ref_date).replace("-", ""), segm, blm_tckrs_lst_str)


            c_d.execute(qry)
            res = c_d.fetchall()
            #print "*"*120
            #print qry
            #print res
            #print "*" * 120
            if ((segm == "CDepositi") or (segm == "CSwap") or (segm == 'CLibor')or (segm == 'CFuture')):
                s = Segm()
                s.name = dict_segm[segm]
                
                if (segm == 'CFuture'):
                    
                    print 'res: ', res
                else:
 
                    print 'res: ', res
                
                for record in res:
                    mat = record[0]
                    val = record[1]
                    s.mats.append(mat)
                    s.values.append(val)
                    self.segms[s.name] = s
            else:
                #spezzo i segmenti in due, gli swap a breve e gli altri
                s1= Segm()
                s1.name = dict_segm[segm]
                s2 = Segm()
                s2.name = 'Swp'
                for record in res:
                    mat = record[0]
                    val = record[1]
                    if (int(mat)<360):
                        s1.mats.append(mat)
                        s1.values.append(val)
                    else:
                        s2.mats.append(mat)
                        s2.values.append(val)
                    self.segms[s1.name] = s1
                    self.segms[s2.name] = s2

        print self.segms
        con.close()

    def bootstrap(self, data_opt):
        data_opt['Basis'  ]    = {}
        data_opt['BusConv']    = {}
        data_opt['RegimeRate'] = {}

        for sn in self.segms.keys():
            name = sn
            code = revDict(dict_segm2) [name]

            s = self.segms[name]
            data_opt['Basis'][code]   = s.getDayCount()
            data_opt['BusConv'][code] = s.getAdj()
            data_opt['RegimeRate'][code]=s.getCapitalization()
        #---
        data_opt['TenorSwap'] = self.floater_tenor
        data_opt['MKT']       = self.cal
        #---
        data_opt['ParConvexity'] = {}
        data_opt['ParConvexity']['A'] = self.HWparms['meanRS']
        data_opt['ParConvexity']['B'] = self.HWparms['sigma']

        data_opt['RefDate'] = self.ref_date
        #==========
        raw_data = {}
        raw_data ['UsaNodo']     = []
        raw_data['Nodo']         = []
        raw_data['ValoreNodo']   = []
        raw_data['TipoSegmento'] = []
        raw_data['MatDate']      = []

        for name in self.segms.keys():
            code = revDict(dict_segm2) [name]
            s = self.segms[name]
            for u,t,v,d in zip(s.usage, s.tags, s.values, s.dates):
                raw_data['UsaNodo'].append(u)
                raw_data['Nodo'].append(t)
                raw_data['ValoreNodo'].append(v)
                raw_data['TipoSegmento'].append(code)
                raw_data['MatDate'].append(d)
        try:
            res = fb.boot3s_elab_v2(data_opt, raw_data)
        except ValueError as ve:
            from Tkinter import *
            import tkMessageBox
            root = Tk()
            root.withdraw()
            msg = ve.message #"Missing input sheet for Swap Curves in your workbook... \nNothing to do for me!"
            tkMessageBox.showinfo("Warning!", msg)
            root.destroy()
            return None
        return res

    def fittingFromPY(self, optDict):

        raw_data = {}
        raw_data['ValoreNodo']   = []
        raw_data['MatDate']      = []

        for name in self.segms.keys():
            code = revDict(dict_segm2) [name]
            s = self.segms[name]
            for u,t,v,d in zip(s.usage, s.tags, s.values, s.dates):

                if (u == 'Y'):
                    raw_data['ValoreNodo'].append(v)
                    raw_data['MatDate'].append(d)
                else:
                    pass
        
        c_dates = raw_data['MatDate']
        c_rates = raw_data['ValoreNodo']
        
        c_rates = np.array(c_rates)/100.0
        
        print 'c_dates: ', c_dates
        print 'c_rates: ', c_rates
        print 'optDict: ', optDict
        
        return fb.fitting(c_dates, c_rates, optDict)


class BootstrappedCurve(Curve):
    def __init__(self):
        Curve.__init__(self)
        self.boot_dates = []
        self.boot_df    = []
        self.fit_usage  = []
        self.boot_rates = []
        self.code       = ""
        self.rendimento = "ZC"

    def show(self):
        Curve.show(self)
        print "------------------------------"
        print "Begin Show Bootstrap Vars:"
        print "Rendimneto:", self.rendimento
        print "Codice Segmenti Boots:", self.code

        for d,f,r,u in zip( self.boot_dates,self.boot_df,self.boot_rates, self.fit_usage):
            print "date:", d, "-- df:", f, "--rate:", r, "-- usage:", u
        print "End Show Bootstrap Vars"
        print "------------------------------"

    def fittingFromBoot(self, optDict):
        
        c_dates = self.boot_dates
        c_rates = self.boot_rates
        
        c_rates = np.array(c_rates)
        c_rates = c_rates/100.0
        
        res_out = fb.fitting(c_dates, c_rates, optDict)
        
        return fb.fitting(c_dates, c_rates, optDict)



class InterpolatedCurve(BootstrappedCurve):
    def __init__(self):
        BootstrappedCurve.__init__(self)
        self.code_interp         = ""
        self.interpolation_model = ""
        self.interp_dates        = []
        self.parms_list          = []
        self.interp_parms        = {}

    def show(self):
        BootstrappedCurve.show(self)
        print "------------------------------"
        print "Begin Show interpolation Vars:"
        print "Codice interpolazione:", self.code_interp
        print "Modello di interpolazione:", self.interpolation_model


        print "date:", self.interp_dates
        print "lista dei parametri del modello:", self.parms_list
        print "dizionario: ", self.interp_parms

        print "End Show interpolation Vars"
        print "------------------------------"

    def getCurveCode(self):
        return self.code_interp



class CdsCurve(Curve):
    def __init__(self):
        Curve.__init__(self)
        self.tags           = []
        self.mats           = []
        self.dates          = []
        self.values         = []
        self.bench_dates    = []
        self.bench_values   = []
        self.bench_df_val   = []
 
        self.infl_ts_val = []
        self.infl_ts_dates = []
 
        self.infl_curve_dates = []
        self.infl_curve_df_val = []
        
        self.bench_type     = ""
        self.risk_free_dates= []
        self.risk_free_df_val  = []
        self.bench_model = ""
        self.bench_prms = {}
        self.ratingProvider = ""
        self.sectorProvider = ""
        self.type           = "CDS"
        self.dayCount       = ""
        self.dayAdj         = ""
        self.lag            = 0
        self.capitalization = ""
        self.code = ""
        self.recovery = 0.0
        self.rendimento = ""
        self.node_type = ""
        self.emittente = ""
        self.cds_boot_method =  ""
        self.hr_model =  ""
        self.rf_interp_type = ""  
   

    def show(self):
        Curve.show(self)
        print "------------------------------"
        print "Begin Show CDS Vars:"
        print "Sector Provider:", self.sectorProvider
        print "Rating Provider:", self.ratingProvider
        print "Curve Type:", self.type
        print "Day count:", self.dayCount
        print "Day Adj:", self.dayAdj
        print "lag (days):", self.lag
        print "Capitalization:", self.capitalization
        print "Curve Code:", self.code
        print "Recovery rate:", self.recovery
        print "HR model:", self.hr_model
        print "Return type:", self.rendimento
        print "Node type:", self.node_type
        # ---
        for t,m,v  in zip( self.tags,self.mats,self.values):
            print "tag:", t, "-- mat:", m, "--value:", v
        # ---
        print "End Show CDS Vars"
        print "------------------------------"


    def mapBootCDS(self, code_m):

        dict_bcode      = {}
        dict_bcode[0] = 'LCS'
        dict_bcode[1] = 'CHR'
        try:        
            code_bf = dict_bcode[code_m]
        except:
            code_bf = 'CHR'
        
        return code_bf

    def mapCodeModelInv(self, code_0):
        
        dict_codeInv = {}
        dict_codeInv['0'] = 'LIN'
        dict_codeInv['1'] = 'AVD'
        dict_codeInv['2'] = 'SVE'
        dict_codeInv['3'] = 'CIR'
        dict_codeInv['4']  = 'NS'
        
        try:        
            code_f = dict_codeInv[code_0]
        except:
            code_f = 'LIN'
        
        return code_f


    def mapCodeModel(self, code_0):
        
        dict_code = {}
        dict_code['LIN'] = '0'
        dict_code['AVD'] = '1'
        dict_code['SVE'] = '2'
        dict_code['CIR'] = '3'
        dict_code['NS']  = '4'
        
        try:        
            code_f = dict_code[code_0]
        except:
            code_f = '0'
        
        return code_f

    def loadDataFromDB(self):
        self.description.strip()
        con = Connection()
        c_a = con.db_data()

        # ---
        #recupero elenco tickers
        # ---
        

        blm_tckrs_lst_str = ""
        sep = ""

        qry = "SELECT DISTINCT  BloombergTicker FROM DProCDS WHERE Descrizione = '%s' " % (
            self.description)
        c_a.execute(qry)
        res = c_a.fetchall()
        
        
        for record in res:
            blm_tckrs_lst_str += sep + "'" + record[0] + "'"
            sep = ","

        
        # ---------
        # RECUPERO DATI DA DB
        # ---------


        c_d = con.db_data()
        qry = '''
                     SELECT DISTINCT currency, Emittente, Seniority  FROM DProCDS where BloombergTicker IN
                     (
                        SELECT DISTINCT BLOOMBERGTICKER FROM DProTS_master WHERE
                        BLOOMBERGTICKER IN (%s) AND
                        DATA = '%s'
                     )
              ''' % (blm_tckrs_lst_str, str(self.ref_date).replace("-",""))
        #print qry
        c_d.execute(qry)
        res = c_d.fetchall()
        if len(res) <> 1: 'AAAA'
        self.curr      = res[0][0]
        code_sen    = res[0][2]
        code_issuer = res[0][1]
        
        
        dict_seniority = {}
        dict_seniority[0] = 'pippo'
        dict_seniority[1] = 'pippo'
        dict_seniority[2] = 'Senior Secured'
        dict_seniority[3] = 'Senior Subordinated'
        dict_seniority[4] = 'Senior Unsecured'
        dict_seniority[5] = 'Covered'
        dict_seniority[6] = 'Junior Subordinated'
        dict_seniority[999] = 'nd'
        
        code_sen_n = dict_seniority[code_sen]
        
        # ---
        # traduzione dei codici per seniority e emittente
        # ---
        """
        qry = '''
                 SELECT descrizione from ZEmittente where codice_emittente = '%s'
             '''%(code_issuer)
        #print qry
        c_a.execute(qry)
        self.emittente = c_a.fetchall()[0][0]
        """
        #print "self.emi", self.emittente

        # ---
        """
        qry = '''
                         SELECT descrizione from ZSeniority where codice_seniority = '%s'
                     ''' % (code_sen)
        #print qry
        c_a.execute(qry)
        self.seniority = c_a.fetchall()[0][0]
        """

        # ---
        # recupero le altre info di anagrafica dalla tabella MKT_Curve, dove sono censite le curve CDS (hanno data_rif == Null)
        # ---
        
        """
        qry = '''
                SELECT codice_curva, tipo_curva, rendimento, tipo_nodo, lag, recovery_rate, codice_segmentazione, codice_settore, codice_rating, periodicita
                FROM MKT_Curve WHERE
                  descrizione = '%s' 
            '''%self.description
        

        #print qry
        
        c_a.execute(qry)
        """
        """
        res             = c_a.fetchall()
        self.code       = res[0][0]
        self.type       = res[0][1]
        self.rendimento = res[0][2]
        self.node_type  = res[0][3]
        self.lag        = res[0][4]
        self.recovery   =  res[0][5]
        cod_segm        = res[0][6] #codice_segmentazione
        self.settore    = res[0][7]
        self.rating     = res[0][8]
        self.frequency  = res[0][9]
        """


        #self.code       = res[0][0]
        #self.code       = res[0][0]

        """
        qry = '''
                SELECT daycount, dayadj, capitalizzazione FROM MKT_Segmentazione_D where codice_segmentazione = '%s'
                '''%cod_segm

        c_a.execute(qry)
        res = c_a.fetchall()

        self.dayCount = res[0][0]
        self.dayAdj   = res[0][1]
        self.capitalization = res[0][2]
        """

        """
        if self.settore == 999:
            qry = '''
                    SELECT   codice_settore from MKT_EmittenteSettore where codice_emittente = '%s'
                    AND date <= '%s'
                    AND FONTE = '%s'
                    order by date desc
                  '''%(code_issuer, str(self.ref_date), self.sectorProvider)
            #print qry
            c_a.execute(qry)
            res = c_a.fetchall()
            self.settore = res[0][0]
        #ora converto il settore da codice a descrizione
        qry = ''' SELECT DESCRIZIONE FROM ZSettore WHERE CODICE_SETTORE = '%s' '''%(self.settore)
        #print qry
        c_a.execute(qry)
        res = c_a.fetchall()
        self.settore = res[0][0]
        """

        """
        if self.rating == 999:
            qry = '''
                         SELECT codice_rating from MKT_EmittenteRating where
                         codice_emittente = '%s' AND
                          DATE <= '%s' and
                          fonte = '%s'
                          order by date desc
                           LIMIT 1
                       ''' % (code_issuer,str(self.ref_date), self.ratingProvider)
            #print qry
            c_a.execute(qry)
            res = c_a.fetchall()
            if len(res) > 1: print 'AAA'
            self.rating = res[0][0]

        qry = ''' SELECT DESCRIZIONE FROM ZRating WHERE CODICE_RATING = '%s' ''' % (self.rating)
        #print qry
        c_a.execute(qry)
        res = c_a.fetchall()
        self.rating = res[0][0]
        """

        # --------------IMPOSTAZIONI DEFAULT 
        self.type       = 'CDS'
        self.rendimento = 'Par Yield'
        self.node_type  = 'CDS'
        self.lag       = 2
        self.recovery  =  40
        cod_segm       = 'SCDS' #codice_segmentazione
        self.settore   = 999
        self.rating    = 999
        self.frequency = 4
        self.seniority = code_sen_n
        
        #----
        if self.curr == "EUR":
            self.cal = "de.eurex"
        elif self.curr == "USD":
            self.cal = 'us'
        elif self.curr == 'GBP':
            self.cal = 'uk'
        elif self.curr == 'CAD':
            self.cal = 'ca'
        else:
            self.cal = 'us'
        # ----
        if self.quotation == "MID":
            quotazione = "valoremid"
        elif self.quotation == "ASK":
            quotazione = "valoreask"
        elif self.quotation == "BID":
            quotazione = "valorebid"
        else:
            print 'AAAA'

        #print res

        qry = '''
                     SELECT      DProCDS.maturityInt, DProTS_master.%s
                     FROM        DProCDS
                     INNER JOIN  DProTS_master
                     ON          DProCDS.BloombergTicker = DProTS_master.BloombergTicker
                     WHERE       DProTS_master.data='%s'
                     AND         DProCDS.BloombergTicker in (%s)
                     ORDER BY    DProCDS.maturityInt
              ''' % (quotazione, str(self.ref_date).replace("-", ""), blm_tckrs_lst_str)
        #print qry
        c_d.execute(qry)
        res = c_d.fetchall()
        #print "*" * 120
        #print qry

        #print res

        for record in res:
            #print "record:", record
            mat = float(record[0])/365.0
            #if (mat-int(mat)) < 1.e-2: mat = int(mat)
            mat = round(mat, 1)
            tag = str(mat)+"Y"
            val = float(record[1])
            #print mat, tag, val
            self.mats.append(mat)
            self.tags.append(tag)
            self.values.append(val)
            #print self.mats, self.tags, self.values

        self.show()



    def loadDataFromDB_old(self):
        self.description.strip()
        con = Connection()
        c_a = con.db_anag()

        # ---
        #recupero elenco tickers
        # ---
        qry = "SELECT DISTINCT  ID_object, Tipo FROM AnagMktObject WHERE description = '%s' " % (
            self.description)
        c_a.execute(qry)
        res = c_a.fetchall()
        if len(res) != 1:
            print 'AAA'

        self.type = res[0][1]
        # ----
        qry = '''
                 SELECT ticker FROM AnagMktObject_D where ID_object ='%s'
              ''' % res[0][0]
        #print qry
        c_a.execute(qry)
        res = c_a.fetchall()
        #print res
        blm_tckrs_lst_str = ""
        sep = ""
        for record in res:
            blm_tckrs_lst_str += sep + "'" + record[0] + "'"
            sep = ","

        # ---------
        # RECUPERO DATI DA DB RRM
        # ---------


        c_d = con.db_data()
        qry = '''
                     SELECT DISTINCT currency, Emittente, Seniority  FROM DProCDS where BloombergTicker IN
                     (
                        SELECT DISTINCT BLOOMBERGTICKER FROM DProTS_master WHERE
                        BLOOMBERGTICKER IN (%s) AND
                        DATA = '%s'
                     )
              ''' % (blm_tckrs_lst_str, str(self.ref_date).replace("-",""))
        #print qry
        c_d.execute(qry)
        res = c_d.fetchall()
        if len(res) <> 1: 'XXX'
        self.curr      = res[0][0]
        code_sen    = res[0][2]
        code_issuer = res[0][1]
        # ---
        # traduzione dei codici per seniority e emittente
        # ---
        qry = '''
                 SELECT descrizione from ZEmittente where codice_emittente = '%s'
             '''%(code_issuer)
        #print qry
        c_a.execute(qry)
        self.emittente = c_a.fetchall()[0][0]
        #print "self.emi", self.emittente

        # ---
        qry = '''
                         SELECT descrizione from ZSeniority where codice_seniority = '%s'
                     ''' % (code_sen)
        #print qry
        c_a.execute(qry)
        self.seniority = c_a.fetchall()[0][0]

        # ---
        # recupero le altre info di anagrafica dalla tabella MKT_Curve, dove sono censite le curve CDS (hanno data_rif == Null)
        # ---
        qry = '''
                SELECT codice_curva, tipo_curva, rendimento, tipo_nodo, lag, recovery_rate, codice_segmentazione, codice_settore, codice_rating, periodicita
                FROM MKT_Curve WHERE
                  descrizione = '%s' 
            '''%self.description

        #print qry
        c_a.execute(qry)
        res             = c_a.fetchall()
        self.code       = res[0][0]
        self.type    = res[0][1]
        self.rendimento = res[0][2]
        self.node_type = res[0][3]
        self.lag     = res[0][4]
        self.recovery =  res[0][5]
        cod_segm     = res[0][6] #codice_segmentazione
        self.settore = res[0][7]
        self.rating  = res[0][8]
        self.frequency = res[0][9]
        qry = '''
                SELECT daycount, dayadj, capitalizzazione FROM MKT_Segmentazione_D where codice_segmentazione = '%s'
                '''%cod_segm

        c_a.execute(qry)
        res = c_a.fetchall()

        self.dayCount = res[0][0]
        self.dayAdj   = res[0][1]
        self.capitalization = res[0][2]


        if self.settore == 999:
            qry = '''
                    SELECT   codice_settore from MKT_EmittenteSettore where codice_emittente = '%s'
                    AND date <= '%s'
                    AND FONTE = '%s'
                    order by date desc
                  '''%(code_issuer, str(self.ref_date), self.sectorProvider)
            #print qry
            c_a.execute(qry)
            res = c_a.fetchall()
            self.settore = res[0][0]
        #ora converto il settore da codice a descrizione
        qry = ''' SELECT DESCRIZIONE FROM ZSettore WHERE CODICE_SETTORE = '%s' '''%(self.settore)
        #print qry
        c_a.execute(qry)
        res = c_a.fetchall()
        self.settore = res[0][0]

        if self.rating == 999:
            qry = '''
                         SELECT codice_rating from MKT_EmittenteRating where
                         codice_emittente = '%s' AND
                          DATE <= '%s' and
                          fonte = '%s'
                          order by date desc
                           LIMIT 1
                       ''' % (code_issuer,str(self.ref_date), self.ratingProvider)
            #print qry
            c_a.execute(qry)
            res = c_a.fetchall()
            if len(res) > 1: zzzzzzzzz
            self.rating = res[0][0]

         # ora converto il settore da codice a descrizione
        qry = ''' SELECT DESCRIZIONE FROM ZRating WHERE CODICE_RATING = '%s' ''' % (self.rating)
        #print qry
        c_a.execute(qry)
        res = c_a.fetchall()
        self.rating = res[0][0]

        #----
        if self.curr == "EUR":
            self.cal = "de.eurex"
        elif self.curr == "USD":
            self.cal = 'us'
        elif self.curr == 'GBP':
            self.cal = 'uk'
        elif self.curr == 'CAD':
            self.cal = 'ca'
        else:
            self.cal = 'us'
        # ----
        if self.quotation == "MID":
            quotazione = "valoremid"
        elif self.quotation == "ASK":
            quotazione = "valoreask"
        elif self.quotation == "BID":
            quotazione = "valorebid"
        else:
            print 'AAAA'

        #print res

        qry = '''
                     SELECT      DProCDS.maturityInt, DProTS_master.%s
                     FROM        DProCDS
                     INNER JOIN  DProTS_master
                     ON          DProCDS.BloombergTicker = DProTS_master.BloombergTicker
                     WHERE       DProTS_master.data='%s'
                     AND         DProCDS.BloombergTicker in (%s)
                     ORDER BY    DProCDS.maturityInt
              ''' % (quotazione, str(self.ref_date).replace("-", ""), blm_tckrs_lst_str)
        #print qry
        c_d.execute(qry)
        res = c_d.fetchall()
        #print "*" * 120
        #print qry

        #print res

        for record in res:
            #print "record:", record
            mat = float(record[0])/360.
            if (mat-int(mat)) < 1.e-8: mat = int(mat)
            tag = str(mat)+"Y"
            val = float(record[1])
            #print mat, tag, val
            self.mats.append(mat)
            self.tags.append(tag)
            self.values.append(val)
            #print self.mats, self.tags, self.values

        self.show()




    def loadBenchDataFromDB(self, opt_download):
        
        self.description.strip()
        con = Connection()
        c_a = con.db_anag()

        ref_date   = opt_download['refDate']
        codeSeg      = opt_download['codeSeg']
        tipo_modello = opt_download['tipo_modello']
        valuta = opt_download['valuta']

        tipo_scarico = 0
        quotazioni = 'MID'
        fonte_dato = 'Bloomberg'
        
        dummy_curva = " (CODICE_CURVA NOT LIKE 'PDELTA%' AND CODICE_CURVA NOT LIKE 'MDELTA%')"
        dummy_model = " (CODICE_MODELLO NOT LIKE 'PDELTA%' AND CODICE_MODELLO NOT LIKE 'MDELTA%')"


        SS1 = 'PDELTA%' 
        SS2 = 'MDELTA%' 

        # RECUPERO CODICE CURVA BENCHMARK            

        qry0 = """
                SELECT CODICE_CURVA FROM MKT_Curve WHERE
                CODICE_EMITTENTE = 999 
                AND FONTE_DATO = '%s'
                AND QUOTAZIONE = '%s'
                AND TIPO_CURVA = 'SWAP' 
                AND TIPO_SCARICO = '%s' 
                AND DATA_RIF = '%s' 
                AND VALUTA = '%s' 
                AND TIPO_NODO = 'Discount Factor'
                AND CODICE_CURVA like '%s'      
                AND CODICE_CURVA not like '%s'
                AND CODICE_CURVA not like '%s'
                
            """  %(fonte_dato, quotazioni, tipo_scarico, ref_date, valuta, codeSeg, SS1, SS2)

        c_a.execute(qry0)
        res = c_a.fetchall()
        
        
        if len(res) != 0:
            code_curve = res[0][0]
        else:
            return 0

        # RECUPERO CURVA BENCHMARK            
        qry1 = "SELECT TERM, VALORE FROM MKT_Curve_D WHERE CODICE_CURVA = '%s' ORDER BY TERM" %code_curve

        c_a.execute(qry1)
        res = c_a.fetchall()
        
        
        
        #print 'res: ', res
        
        valueList = []
        matDateList = []
        
        for i in range(0, len(res)):
            
            resDateTmp = res[i][0]
            resValueTmp = res[i][1]
            
            matDateList.append(resDateTmp)
            valueList.append(resValueTmp)
        
        self.bench_dates = matDateList
        self.bench_df_val = valueList
        
        
        #print 'self.bench_dates: ', self.bench_dates
        #print 'ln: ', len(self.bench_dates)

        #fb.FQ(999)
        
        #------------------------------------------
        
        
        qry2 = """
            SELECT CODICE_MODELLO, DESCRIZIONE FROM MKT_ParametriInterp
            WHERE CODICE_EMITTENTE = 999 
            AND FONTE_DATO = '%s'
            AND QUOTAZIONE = '%s'
            AND TIPO_CURVA = 'SWAP' 
            AND TIPO_SCARICO = '%s'
            AND DATA_RIF = '%s'
            AND VALUTA = '%s'
            AND TIPO_MODELLO = '%s'
            AND CODICE_MODELLO like '%s'
            AND CODICE_MODELLO not like '%s'
            AND CODICE_MODELLO not like '%s'
            
        
        """  %(fonte_dato, quotazioni, tipo_scarico, ref_date, valuta, tipo_modello, codeSeg, SS1, SS2)
        c_a.execute(qry2)
        res = c_a.fetchall()



        if len(res) != 0:
            codice_modello = res[0][0]
            des_modello  = res[0][1]
        else:
            return 0
        #print 'code_modello: ', code_modello
        #print 'des_modello: ', des_modello
        
        qry3 = """
                SELECT VALORE FROM MKT_ParametriInterp_D WHERE 
                CODICE_MODELLO = '%s'
                ORDER BY CODICE_PARAMETRO
            """  %(codice_modello)

        c_a.execute(qry3)
        res = c_a.fetchall()

        if len(res) != 0:
            pass
        else:
            return 0        
        
        mdl_prms_list = []
        for i in range(0, len(res)):
            
            parTmp = res[i][0]
            mdl_prms_list.append(parTmp)
        
        
        prms_dict = fb.packModelPrms(tipo_modello, self.bench_dates, mdl_prms_list)
        
        self.bench_prms = prms_dict       
        self.bench_model = tipo_modello
        
        return 1

    def loadBenchDataFromDB_Bond(self, opt_download):

        self.description.strip()
        con = Connection()
        c_a = con.db_anag()

        ref_date = opt_download['refDate']
        codeSeg = opt_download['codeSeg']
        tipo_modello = opt_download['tipo_modello']
        valuta = opt_download['valuta']

        tipo_scarico = 0
        quotazioni = 'MID'
        fonte_dato = 'Bloomberg'

        dummy_curva = " (CODICE_CURVA NOT LIKE 'PDELTA%' AND CODICE_CURVA NOT LIKE 'MDELTA%')"
        dummy_model = " (CODICE_MODELLO NOT LIKE 'PDELTA%' AND CODICE_MODELLO NOT LIKE 'MDELTA%')"

        SS1 = 'PDELTA%'
        SS2 = 'MDELTA%'

        # RECUPERO CODICE CURVA BENCHMARK

        qry0 = """
                SELECT CODICE_CURVA FROM MKT_Curve WHERE
                CODICE_EMITTENTE = 999 
                AND FONTE_DATO = '%s'
                AND QUOTAZIONE = '%s'
                AND TIPO_CURVA = 'SWAP' 
                AND TIPO_SCARICO = '%s' 
                AND DATA_RIF = '%s' 
                AND VALUTA = '%s' 
                AND TIPO_NODO = 'Discount Factor'
                AND CODICE_CURVA like '%s'      
                AND CODICE_CURVA not like '%s'
                AND CODICE_CURVA not like '%s'

            """ % (fonte_dato, quotazioni, tipo_scarico, ref_date, valuta, codeSeg, SS1, SS2)

        c_a.execute(qry0)
        res = c_a.fetchall()

        if len(res) != 0:
            code_curve = res[0][0]
        else:
            return 0

        # RECUPERO CURVA BENCHMARK
        qry1 = "SELECT TERM, VALORE FROM MKT_Curve_D WHERE CODICE_CURVA = '%s' ORDER BY TERM" % code_curve

        c_a.execute(qry1)
        res = c_a.fetchall()

        # print 'res: ', res

        valueList = []
        matDateList = []

        for i in range(0, len(res)):
            resDateTmp = res[i][0]
            resValueTmp = res[i][1]

            matDateList.append(resDateTmp)
            valueList.append(resValueTmp)

        self.bench_dates = matDateList
        self.bench_df_val = valueList

        # ------------------------------------------
        # NON RECUPERO I PARAMETRI DEL MODELLO  PER I BOND
        # IL BOND FITTING UTILIZZA SEMPRE UN MODELLO DI TIPO LINEARE
        # CHE VIENE RICALCOLATO DALLA ROUTINE
        # ------------------------------------------

        if tipo_modello == 'LIN':
            pass
        else:
            return 2

        '''
        qry2 = """
            SELECT CODICE_MODELLO, DESCRIZIONE FROM MKT_ParametriInterp
            WHERE CODICE_EMITTENTE = 999 
            AND FONTE_DATO = '%s'
            AND QUOTAZIONE = '%s'
            AND TIPO_CURVA = 'SWAP' 
            AND TIPO_SCARICO = '%s'
            AND DATA_RIF = '%s'
            AND VALUTA = '%s'
            AND TIPO_MODELLO = '%s'
            AND CODICE_MODELLO like '%s'
            AND CODICE_MODELLO not like '%s'
            AND CODICE_MODELLO not like '%s'


        """ % (fonte_dato, quotazioni, tipo_scarico, ref_date, valuta, tipo_modello, codeSeg, SS1, SS2)
        c_a.execute(qry2)
        res = c_a.fetchall()

        if len(res) != 0:
            codice_modello = res[0][0]
            des_modello = res[0][1]
        else:
            return 0
        # print 'code_modello: ', code_modello
        # print 'des_modello: ', des_modello

        qry3 = """
                SELECT VALORE FROM MKT_ParametriInterp_D WHERE 
                CODICE_MODELLO = '%s'
                ORDER BY CODICE_PARAMETRO
            """ % (codice_modello)

        c_a.execute(qry3)
        res = c_a.fetchall()

        if len(res) != 0:
            pass
        else:
            return 0

        mdl_prms_list = []
        for i in range(0, len(res)):
            parTmp = res[i][0]
            mdl_prms_list.append(parTmp)
        prms_dict = fb.packModelPrms(tipo_modello, self.bench_dates, mdl_prms_list)
        '''

        self.bench_prms = {} #prms_dict
        self.bench_model = tipo_modello

        return 1

    def loadInflCurveFromDB(self, opt_download):
        
        self.description.strip()
        con = Connection()
        c_a = con.db_anag()

        ref_date   = opt_download['refDate']
        valuta = opt_download['valuta']

        tipo_scarico = 0
        quotazioni = 'MID'
        fonte_dato = 'Bloomberg'


        # RECUPERO CODICE CURVA BENCHMARK            

        qry0 = """
                SELECT CODICE_CURVA FROM MKT_Curve WHERE
                CODICE_EMITTENTE = 999 
                AND FONTE_DATO = '%s'
                AND QUOTAZIONE = '%s'
                AND TIPO_CURVA = 'Inflazione' 
                AND TIPO_SCARICO = '%s' 
                AND DATA_RIF = '%s' 
                AND VALUTA = '%s' 
                AND TIPO_NODO = 'Discount Factor'                
            """  %(fonte_dato, quotazioni, tipo_scarico, ref_date, valuta)

        c_a.execute(qry0)
        res = c_a.fetchall()
        
        if len(res) != 0:
            code_curve = res[0][0]
        else:
            print 'Curva di inflazione non trovata nel database!!!!'
            return 0

        # RECUPERO CURVA BENCHMARK            
        qry1 = "SELECT TERM, VALORE FROM MKT_Curve_D WHERE CODICE_CURVA = '%s' ORDER BY TERM" %code_curve

        c_a.execute(qry1)
        res = c_a.fetchall()
        #print 'res: ', res
        
        valueList = []
        matDateList = []
        
        for i in range(0, len(res)):
            
            resDateTmp = res[i][0]
            resValueTmp = res[i][1]
            
            matDateList.append(resDateTmp)
            valueList.append(resValueTmp)
        
        self.infl_curve_dates = matDateList
        self.infl_curve_df_val = valueList
        
        return 1
        
    def loadInflTSFromDB(self, opt_download):
        
        self.description.strip()
        con = Connection()
        #c_a = con.db_anag()
        c_a = con.db_data()

        valuta = opt_download['valuta']

        # RECUPERO SERIE STORICA TS            

        qry0 = """
                SELECT DATA, VALORE FROM TS_master WHERE
                AssetClass = 'CPTFEMU' 
                AND CodicePaese = '%s'
            """  %(valuta)

        c_a.execute(qry0)
        res = c_a.fetchall()
        
        if len(res) != 0:

            valueList = []
            matDateList = []
            
            for i in range(0, len(res)):
                
                resDateTmp = res[i][0]
                resValueTmp = res[i][1]
                
                matDateList.append(resDateTmp)
                valueList.append(resValueTmp)
            
            
        else:
            return 0

        self.infl_ts_dates = matDateList
        self.infl_ts_val = valueList

        
        return 1
        
        



    def bootstrap(self, data_opt):
        
        data_opt['Basis'  ]    = {}
        data_opt['BusConv']    = {}
        data_opt['RegimeRate'] = {}

        #================================

        data_opt['MKT']     = self.cal
        data_opt['RefDate'] = self.ref_date
        
        self.hr_model
        
        print 'AAAAAAAAAAAAAAAAA'
        print 'self.recovery: ', self.recovery
        print 'self.capitalization: ', self.capitalization
        print 'self.hr_model: ', self.hr_model
        print self.cds_boot_method
        

        print 'AAAAAAAAAAAAAAAAA'

        if (self.curr == 'EUR'):
            tenor = 3 # 3 pagamenti ogni 3Mesi
        else:
            tenor = 3

        data_opt['currency']       = self.curr # 
        data_opt['tenor']          = tenor # 
        data_opt['Basis']          = 'ACT/365'
        data_opt['interp']         = self.hr_model # 
        data_opt['BusConv']        = 'modfollow'
        data_opt['fixingDays']     = 2
        data_opt['compounding']    = 0  #0 = semplice, 1 = composto, 2 = continuo
        data_opt['ReocveryRate']   = self.recovery
        data_opt['hr_bootMethod']  = self.hr_model #0 = LCS, 1 = CHR
        
        
         
         
        data_raw_cds = {}
        data_raw_cds['ValoreNodo']  = self.values 
        data_raw_cds['MatTimes']    = self.mats
        data_raw_cds['Nodo']        = self.tags

        data_raw_bench = {}
        data_raw_swp = {}


        
        data_raw_bench['MatDate']    =self.bench_dates
        data_raw_bench['ValoreNodo'] =self.bench_values
        
        data_raw_bench['DiscountFactor'] =self.bench_df_val
        data_raw_bench['Model'] = self.mapCodeModel(self.bench_model)
        data_raw_bench['prms'] =self.bench_prms
        data_raw_bench['Type'] = self.bench_type
        
        data_raw_swp['MatDate']    =self.risk_free_dates
        data_raw_swp['DiscountFactor'] =self.risk_free_df_val

        data_raw_swp = data_raw_bench

        

        try:
            
            import funzioni_boot_cds as f_cds
            
            
            res = f_cds.boot_cds(data_opt, data_raw_cds, data_raw_bench, data_raw_swp)
            #res = fb.boot3s_elab_v2(data_opt, raw_data)

        except ValueError as ve:
            from Tkinter import *
            import tkMessageBox
            root = Tk()
            root.withdraw()
            msg = ve.message #"Missing input sheet for Swap Curves in your workbook... \nNothing to do for me!"
            tkMessageBox.showinfo("Warning!", msg)
            root.destroy()
            return None
        return res


    def getCurveCode(self):
        codeL   = "C"
        codeL  += self.curr
        codeR   = "ZC"
        codeR  += "BLM0" if self.source=='Bloomberg'  else "OTH"
        dd      = "%02d"%self.ref_date.day
        mm      = "%02d"%self.ref_date.month
        yy      = (str(self.ref_date.year))
        codeR  += dd+mm+yy[-2:]
        return codeL, codeR


