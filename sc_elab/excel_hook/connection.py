import sys
import os
import pyodbc


class Connection:
    def __init__(self):
        self.db_dta = None
        self.db_ang = None

    def db_data(self):
        if not self.db_dta:
            self.db_dta = (pyodbc.connect(r'DSN=dati_closing;UID=ricky;PWD=balboa')).cursor()
        return self.db_dta

    def db_anag(self):

        if not self.db_ang:
            self.db_ang  = (pyodbc.connect(r'DSN=penelope;UID=pricingunit;PWD=pr1cing4')).cursor()

        return self.db_ang

    def db_save(self):
        return self.db_anag()

    def close(self):
        if self.db_dta: self.db_dta.close()
        del(self.db_dta)

        if self.db_ang: self.db_ang.close()
        del (self.db_ang)
