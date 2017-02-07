import sys
import os



dict_segm =\
    {\
     "CDepositi": "Dep" \
    , "CLibor"  : "Libor"\
    , "CFuture" : "Fut" \
    , "CSwap"   : "Swp" \
}

class Curve:

    def __init__(self):
        #anagrafica
        self.description = ""
        self.curr = ""
        self.ref_date = ""
        self.type = ""
        self.source = "Bloomberg"
        self.quotation = "MID"
        self.download_type = ""
        self.emittente = '999'
        self.rating = 'NR'
        self.settore = '999'
        self.seniority = '999'
        self.type = 'Swap'
        self.floater_tenor = ''
        #dati grezzi
        self.raw_data = {"Dep":[], "Libor":[], "Fut":[], "Swp":[]}

    #def load_data_from_db (self):
    #    return True


