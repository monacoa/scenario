import sys
import os
import pyodbc

import pymysql
import mysql.connector



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

            db_credential = {}
            db_credential['host'] = "localhost"                         
            db_credential['user'] = "xxx3"         
            db_credential['pwd']  = "yyy2"                        
            db_name = 'test_db_mkt_data'

    
            self.db = mysql.connector.connect(host     = db_credential['host'],
                                         user     = db_credential['user'],
                                         passwd   = db_credential['pwd'],
                                         database = db_name)

            self.db_dta = self.db.cursor()


        return self.db_dta
    




    def db_anag(self):

        if not self.db_ang:
            self.db  = (pyodbc.connect(r'DSN=penelope;UID=pricingunit;PWD=pr1cing4'))
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
