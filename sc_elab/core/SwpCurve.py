import sys
import os
from sc_elab.excel_hook.connection import *
import sc_elab.core.mdates.busdayrule as busD
import sc_elab.core.mdates.holidays as holy
import sc_elab.core.mdates.dateutils as dateutils
import datetime
from dateutil.relativedelta import relativedelta

dict_segm =\
    {\
      "CDepositi": "Dep" \
    , "CLibor"  : "Libor"\
    , "CFuture" : "Fut" \
    , "CSwap"   : "Swp" \
    }

dict_segm2 =\
    {\
      "D"   : "Dep"     \
    , "L"   : "Libor"   \
    , "F"   : "Fut"     \
    , "S"   : "Swp"     \
    , "G"   : "ST_Swp"  \
    }



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



class Curve:

    def __init__(self):
        #--------
        #anagrafica
        #--------
        self.description    = ""
        self.curr           = ""
        self.ref_date       = ""
        self.type           = ""
        self.source         = "Bloomberg"
        self.quotation      = "MID"
        self.download_type  = ""
        self.emittente      = '999'
        self.rating         = 'NR'
        self.settore        = '999'
        self.seniority      = '999'
        self.type           = 'Swap'
        self.floater_tenor  = ''
        self.cal            = ''
        #---------
        self.HWparms          = {}
        #----------
        #segmenti (dict of classes
        #----------
        self.segms          = {}


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
        #data = self.raw_data
        #dd = []
        con = Connection()
        cn = con.db_anag()

        for k in self.segms.keys():

            seg = self.segms[k]
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
                # case s[0] == 'D','L','S', 'G'
                rdate = dateutils.asdatetime(self.ref_date)
                for mat,tag in zip (self.segms[s].mats, self.segms[s].tags):
                    if   tag == "O/N" : date = self.addWorkingDays(rdate , 1, adj)
                    elif tag == "T/N" : date = self.addWorkingDays(rdate, 2, adj)
                    elif tag[-1]=="D" :
                        nday = int(tag[:-1])
                        date = self.addWorkingDays(rdate, 2, adj)
                    elif tag[-1] == "W":
                        date = rdate + datetime.timedelta(weeks=(int(tag[:-1])))
                        calendar = holy.get_calendar(self.cal)
                        date = busD.rolldate_from_db(date, calendar, adj)
                    elif tag[-1] == "M":
                        calendar = holy.get_calendar(self.cal)
                        date = rdate + relativedelta(months=(int(tag[:-1])))
                        date = busD.rolldate_from_db(date, calendar, adj)
                    elif tag[-1] == "Y":
                        calendar = holy.get_calendar(self.cal)
                        date = rdate + relativedelta(years=(int(tag[:-1])))
                        date = busD.rolldate_from_db(date, calendar, adj)
                    else: qqqqqqqqqqqqqqqqqqqqqqqq
                    self.segms[s].dates.append(date)



    def init_finalize(self):
        #compute output strings
        self.computeTags()
        #compute anag segms
        self.fillAnagSegm()
        #compute output dates
        self.computeDates()


    def loadDataFromDB(self):
        con = Connection()
        c_a = con.db_anag()
        qry = "SELECT DISTINCT  ID_object, Tipo, rating FROM AnagMktObject WHERE description = '%s' " % (
        self.description)

        c_a.execute(qry)
        res = c_a.fetchall()

        if len(res) != 1:
            aaaaaaaaaaaaaaaaaaaaaaa
        self.rating = res[0][2]
        self.type = res[0][1]
        # ----
        qry = '''
                SELECT ticker FROM AnagMktObject_D where ID_object ='%s'
                ''' % res[0][0]

        c_a.execute(qry)
        res = c_a.fetchall()
        blm_tckrs_lst_str = ""
        sep = ""
        for record in res:
            blm_tckrs_lst_str += sep + "'" + record[0] + "'"
            sep = ","

        # ---------
        # RECUPERO CURRENCY E NUMERO/TIPOLOGIA DI SEGMENTI
        # ---------

        c_d = con.db_data()
        qry = '''
             SELECT DISTINCT currency, TipoDato FROM alm.DProCurve where BloombergTicker IN
             (
                SELECT DISTINCT BLOOMBERGTICKER FROM alm.DProTS_master WHERE
                BLOOMBERGTICKER IN (%s) AND
                 DATA = '%s'
                )
            ''' % (blm_tckrs_lst_str, str(self.ref_date).replace("-", ""))

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
                        SELECT      alm.DProCurve.maturityInt, alm.DProTS_master.%s
                        FROM        alm.DProCurve
                        INNER JOIN  alm.DProTS_master
                        ON          alm.DProCurve.BloombergTicker = alm.DProTS_master.BloombergTicker
                        WHERE       alm.DProTS_master.data='%s'
                        AND         alm.DProCurve.TipoDato = '%s'
                        AND         alm.DProCurve.BloombergTicker in (%s)
                        ORDER BY    alm.DProCurve.maturityInt
                        ''' % (quotazione, str(self.ref_date).replace("-", ""), segm, blm_tckrs_lst_str)
            elif segm == "CFuture":
                qry = '''
                        SELECT      alm.DProTS_master.ScadenzaFuture, alm.DProTS_master.%s
                        FROM        alm.DProCurve
                        INNER JOIN  alm.DProTS_master
                        ON          alm.DProCurve.BloombergTicker = alm.DProTS_master.BloombergTicker
                        WHERE       alm.DProTS_master.data= '%s'
                        AND         alm.DProCurve.TipoDato ='%s'
                        AND         alm.DProCurve.BloombergTicker in (%s)
                        ORDER BY    alm.DProTS_master.ScadenzaFuture
                        ''' % (quotazione, str(self.ref_date).replace("-", ""), segm, blm_tckrs_lst_str)

            else:
                xxxxxxxxxxxxxxxxxxxxxxxx


            c_d.execute(qry)
            res = c_d.fetchall()

            s = Segm()
            s.name = dict_segm[segm]
            for record in res:
                mat = record[0]
                val = record[1]
                s.mats.append(mat)
                s.values.append(val)

            self.segms[s.name] = s

        con.close()
