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





class BondPortfolio(object):
    def __init__(self):
        object.__init__(self)
        self.setDefaults()

    def setDefaults(self):
        # --------
        # anagrafica
        # --------
        
        self.bond_portfolio = []
        
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
        self.recoveryRate = ''
        self.ref_date = ''
        # ---------
        # ----------
        # segmenti (dict of classes
        # ----------



    """
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

    """
    
    """
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
    """

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


    def getCodiceSettore(self, m_Emittente, bond_date):
        
        c_date_qry = str(bond_date).replace("-", "")
    
        qry2 = '''
            SELECT CODICE_SETTORE FROM MKT_EmittenteSettore
            WHERE CODICE_EMITTENTE ='%s'  AND DATE <= '%s'
            AND FONTE = 'Bloomberg'  ORDER BY DATE DESC
            '''%(m_Emittente, c_date_qry)
    
        con2 = Connection()
        c_anag = con2.db_anag()
        
        c_anag.execute(qry2)
        res2 = c_anag.fetchone()
        con2.close()
        
        return  res2[0]


    def getDescrzioneSettore(self, cod_sector):
    
        qry = ''' SELECT DESCRIZIONE FROM ZSettore WHERE CODICE_SETTORE = '%s' '''%(cod_sector)
    
        con3 = Connection()
        c_a = con3.db_anag()
        
        c_a.execute(qry)
        res = c_a.fetchone()
        return res[0]





    def getDescrzioneEmittente(self, cod_emittente):
        
            qry = ''' SELECT DESCRIZIONE FROM ZEmittente WHERE CODICE_EMITTENTE = '%s' '''%(cod_emittente)
        
            con3 = Connection()
            c_a = con3.db_anag()
            
            c_a.execute(qry)
            res = c_a.fetchone()
            con3.close()
    
            return res[0]


    def getDescrzioneRating(self, cod_rating):
        
            qry = ''' SELECT DESCRIZIONE FROM ZRating WHERE CODICE_RATING = '%s' '''%(cod_rating)
        
            con3 = Connection()
            c_a = con3.db_anag()
            
            c_a.execute(qry)
            res = c_a.fetchone()
            con3.close()
            
            return res[0]

    
    def getCodeSeniority(self, des_seniority):
        
            qry = ''' SELECT CODICE_SENIORITY FROM ZSeniority WHERE DESCRIZIONE = '%s' '''%(des_seniority)
        
            con3 = Connection()
            c_a = con3.db_anag()
            
            c_a.execute(qry)
            res = c_a.fetchone()
            con3.close()
            
            return res[0]


    
    def getRecoveryValue(self, cod_recovery):
        
            #print 'cod_rating: ', cod_recovery
            qry = ''' SELECT VALORE FROM MKT_RecoveryRate_D WHERE CODICE_RR = '%s' '''%(cod_recovery)
        
            con3 = Connection()
            c_a = con3.db_anag()
            
            c_a.execute(qry)
            res = c_a.fetchone()
            con3.close()
            
            
            try:
                rv_out = res[0]
            except:
                rv_out = 99999

            
            return rv_out



    def getCodiceRating(self, bond_date, cod_emittente):
    
        c_date_qry = str(bond_date).replace("-", "")
        
        fonteRating = 'StandardPoors18'
        
        qry3 = '''
            SELECT CODICE_RATING FROM MKT_EmittenteRating WHERE CODICE_EMITTENTE = '%s'
              AND DATE <= '%s' ORDER BY DATE DESC''' %(cod_emittente, c_date_qry)
    
        con = Connection()
        c_anag = con.db_anag()
    
            
        c_anag.execute(qry3)
        res = c_anag.fetchall()
        con.close()
        
        return res[0][0]

    def getCodiceRecovery(self, cod_seniority, cod_sector, bond_date):
        
        c_date_qry = str(bond_date).replace("-", "")
    
        cod_settore = cod_sector
        c_date_qry = c_date_qry
        fonte_rr = 'JPMorgan'
        
    
        sql = '''
                SELECT  CODICE_RR from MKT_RecoveryRate 
                WHERE CODICE_SENIORITY = '%s'
                AND CODICE_SETTORE = '%s'
                AND CODICE_RATING = '999' 
                AND DATA_RIF <= '%s'
                AND FONTE_DATO = '%s'
                ORDER BY DATA_RIF DESC
            ''' %(cod_seniority, cod_settore, c_date_qry, fonte_rr)
    
    
        con = Connection()
        c_anag = con.db_anag()
    
        c_anag.execute(sql)
        res = c_anag.fetchall()
        c_anag.close()

        
        if len(res) > 1:
            
            r_out = res[0][0]
        else:
            try:
                r_out = res[0]
            except:
                r_out = 999
                
                
        return r_out
        
    


    
    def loadDataFromDB(self, bond_date):
        
        
        con = Connection()
        c_a = con.db_data()

        c_date_qry = str(bond_date).replace("-", "")

    
        qry = ''' 
                SELECT * FROM Bond_master WHERE Data = '%s'  AND descrizione = '%s'  
                ORDER BY Scadenza ASC
            ''' %(c_date_qry, self.description)

        c_a.execute(qry)
        res = c_a.fetchall()
        
        
        
        
        num_fields = len(c_a.description)
        field_names = [i[0] for i in c_a.description]
        
        bond_portfolio = []
        
        for row in res:
            bond_portfolio.append(dict(zip(field_names, row)))

        self.bond_portfolio = bond_portfolio
        con.close()
        
        self.seniority = self.bond_portfolio[0]['PUCollatType']
        code_emittente = self.bond_portfolio[0]['Emittente']
        self.curr = self.bond_portfolio[0]['Divisa']
        
        
        cod_seniority = self.getCodeSeniority(self.seniority)

        
        cod_sector = self.getCodiceSettore(code_emittente, bond_date)
        des_cod_sector = self.getDescrzioneSettore(cod_sector)
        des_emittente = self.getDescrzioneEmittente(code_emittente)

        cod_rating    = self.getCodiceRating(bond_date, code_emittente)
        des_cod_rating = self.getDescrzioneRating(cod_rating)    

        cod_rr   = self.getCodiceRecovery(cod_seniority, cod_sector, bond_date)    
        value_rr = self.getRecoveryValue(cod_rr)
        
        
        
        self.settore = des_cod_sector
        self.emittente = des_emittente
        self.rating = des_cod_rating
        self.recoveryRate = float(value_rr)/100.0
        
        
        
        
        


    """
    sql = "SELECT CODICE_SETTORE FROM MKT_EmittenteSettore" & _
          " WHERE CODICE_EMITTENTE =" & m_Emittente & _
          " AND DATE <= " & Apici(PWDate(CDate(DataLettura))) & _
          " AND FONTE = " & Apici(fonte_settore) & _
          " ORDER BY DATE DESC"

    """
          
    """
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

    """
    
    """
    def fittingFromPY(self, optDict):

        raw_data = {}
        raw_data['ValoreNodo']   = []
        raw_data['MatDate']      = []

        for name in self.segms.keys():
            code = revDict(dict_segm2) [name]
            s = self.segms[name]
            for u,t,v,d in zip(s.usage, s.tags, s.values, s.dates):

                if (u == 'y'):
                    raw_data['ValoreNodo'].append(v)
                    raw_data['MatDate'].append(d)
                else:
                    pass
        
        c_dates = raw_data['MatDate']
        c_rates = raw_data['ValoreNodo']
        
        c_rates = np.array(c_rates)/100.0
        
        return fb.fitting(c_dates, c_rates, optDict)
    """
