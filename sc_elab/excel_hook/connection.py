import sys
import os
import pyodbc

import pymysql

class Connection:
    def __init__(self):
        self.db = None
        self.db_dta = None
        self.db_ang = None

    """
    def db_data(self):
        if not self.db_dta:
            self.db = (pyodbc.connect(r'DSN=dati_closing;UID=ricky;PWD=balboa'))
            self.db_dta = self.db.cursor()
        return self.db_dta
    """

    
    def db_data(self):

        if not self.db_dta:
            self.db =  pymysql.connect(host='cor-dbserver-mysql.prometeia',
                                         port = 3306,
                                         user='pricingunit',
                                         password='pr1cing4',
                                         db='db_mkt_data')

            #self.db  = (pyodbc.connect(r'DSN=db_mercato;UID=root;PWD=lucap'))
            #self.db  = (pyodbc.connect(r'DSN=db_mkt_data;UID=pricingunit;PWD=pr1cing4'))

            ### gestire l'errore nel caso in cui non riesca a connettersi

            self.db_dta = self.db.cursor()
        return self.db_dta


    def db_anag(self):

        if not self.db_ang:
            self.db =  pymysql.connect(host='cor-dbserver-mysql.prometeia',
                                         port = 3306,
                                         user='pricingunit',
                                         password='pr1cing4',
                                         db='penelope')

            #self.db  = (pyodbc.connect(r'DSN=penelope;UID=pricingunit;PWD=pr1cing4'))
            self.db_ang = self.db.cursor()

        return self.db_ang

    def db_save(self):
        return self.db_anag()

    def close(self):
        if self.db_dta: self.db_dta.close()
        del(self.db_dta)

        if self.db_ang: self.db_ang.close()
        del (self.db_ang)

    def commit(self):
        self.db.commit()
