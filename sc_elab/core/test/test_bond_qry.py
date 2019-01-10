import sys
from sc_elab.excel_hook.connection import *


def FQ(label):
    print ('------------- FIN QUI TUTTO OK  %s ----------' %(label))
    sys.exit()



def getBondListFromDb(bond_date, type):


    c_date_qry = str(bond_date).replace("-", "")
    con = Connection()
    
    
    #c_anag = con.db_anag()
    c_anag = con.db_data()
    

    
    """    
    qry = '''
          SELECT DISTINCT description, ID_object FROM AnagMktObject WHERE ID_object in
          ( SELECT DISTINCT ID_object FROM AnagMktObject_D WHERE ticker IN
            (%s
            )
          )ORDER BY description
          ''' % tkrs

    """
    qry = '''
          SELECT DISTINCT Descrizione FROM view_Bond_master  WHERE Data = '%s'" 
          ''' %(c_date_qry)
          
          
    print "qry:", qry
    c_anag.execute(qry)
    res = c_anag.fetchall()
    print 'risultato qry: ', res
    ll = []
    for c in res:
        print "c", c
        curva = str(c[0])
        id    = int(c[1])
        a =[curva,id]
        ll.append(a)
        print "ll", ll
    con.close()
    return ll
    # -----------

def getCodiceSettore(m_Emittente, bond_date):
    
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


def getdescrzioneSettore(cod_sector):

    qry = ''' SELECT DESCRIZIONE FROM ZSettore WHERE CODICE_SETTORE = '%s' '''%(cod_sector)

    con3 = Connection()
    c_a = con3.db_anag()
    
    c_a.execute(qry)
    res = c_a.fetchone()
    con3.close()


    return res[0]

def getDescrzioneEmittente(cod_emittente):
    
        qry = ''' SELECT DESCRIZIONE FROM ZEmittente WHERE CODICE_EMITTENTE = '%s' '''%(cod_emittente)
    
        con3 = Connection()
        c_a = con3.db_anag()
        
        c_a.execute(qry)
        res = c_a.fetchone()
        con3.close()

        return res[0]


def getDescrzioneRating(cod_rating):
    
        qry = ''' SELECT DESCRIZIONE FROM ZRating WHERE CODICE_RATING = '%s' '''%(cod_rating)
    
        con3 = Connection()
        c_a = con3.db_anag()
        
        c_a.execute(qry)
        res = c_a.fetchone()
        con3.close()
        
        return res[0]


def getCodeSeniority(des_seniority):
    
        qry = ''' SELECT CODICE_SENIORITY FROM ZSeniority WHERE DESCRIZIONE = '%s' '''%(des_seniority)
    
        con3 = Connection()
        c_a = con3.db_anag()
        
        c_a.execute(qry)
        res = c_a.fetchone()
        con3.close()
        
        return res[0]



def getRecoveryValue(cod_recovery):
    
        #print 'cod_rating: ', cod_recovery
        qry = ''' SELECT VALORE FROM MKT_RecoveryRate_D WHERE CODICE_RR = '%s' '''%(cod_recovery)
    
        con3 = Connection()
        c_a = con3.db_anag()
        
        c_a.execute(qry)
        res = c_a.fetchone()
        con3.close()
        
        return res[0]



def retriveCodiceRating(bond_date, cod_emittente):

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

def getCodiceRecovery(cod_seniority, cod_sector, bond_date):
    
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
    
    return res[0][0]
    


if __name__ == "__main__":
    




    # ----------- DATA ETST ----------------------    
    bond_date = '2015-12-30'
    des = 'Bond Italia GOV Inflation EUR'
    
    c_date_qry = str(bond_date).replace("-", "")
    
    

    qry = ''' 
            SELECT * FROM ricky.view_Bond_master WHERE Data = '%s' and Descrizione = '%s'
            ORDER BY Scadenza ASC
        ''' %(c_date_qry, des)
          
          
    print "qry:", qry

    con = Connection()
    c_anag = con.db_data()

    c_anag.execute(qry)
    res = c_anag.fetchall()
    
    num_fields = len(c_anag.description)
    field_names = [i[0] for i in c_anag.description]
    
    print 'field_names: ', field_names
    
    bond_portfolio = []
    
    for row in res:
        bond_portfolio.append(dict(zip(field_names, row)))
    
    con.close()
    
    
    fieldRef = 'Divisa'


    print 'CHKK value: ', bond_portfolio[0][fieldRef]
    print '---------------------------------------------'
    
    FQ(999)
    
    
    
    
    #self.bond_portfolio = bond_portfolio
    #columns = [column[0] for column in cursor.description]
    
    
    cod_emittente = bond_portfolio[0]['Emittente']
    cod_sector = getCodiceSettore(cod_emittente, bond_date)
    des_sector = getdescrzioneSettore(cod_sector)
    
    print 'cod_sector: ', cod_sector
    print 'des_sector: ', des_sector
    print '-------------------------------'
    
    des_emittente = getDescrzioneEmittente(cod_emittente)

    print 'cod_emittente: ', cod_emittente
    print 'des_emittente: ', des_emittente
    print '-------------------------------'
    
    
    cod_rating = retriveCodiceRating(bond_date, cod_emittente)

    des_cod_rating = getDescrzioneRating(cod_rating)    

    print 'cod_rating: ',cod_rating
    print 'des_cod_rating: ',des_cod_rating
    print '-------------------------------'
    
    des_seniority = bond_portfolio[0]['PUCollatType']
    cod_seniority = getCodeSeniority(des_seniority)

    print 'des_seniority: ',des_seniority
    print 'cod_seniority: ',cod_seniority
    print '-------------------------------'
    

    cod_rr   = getCodiceRecovery(cod_seniority, cod_sector, bond_date)    
    value_rr = getRecoveryValue(cod_rr)

    print 'cod_recovery: ', cod_rr
    print 'value_rr: ', value_rr

    


    cod_emittente = bond_portfolio[0]['Emittente']
    cod_sector = getCodiceSettore(cod_emittente, bond_date)
    des_sector = getdescrzioneSettore(cod_sector)
    
    des_emittente = getDescrzioneEmittente(cod_emittente)
    cod_rating = retriveCodiceRating(bond_date, cod_emittente)

    des_cod_rating = getDescrzioneRating(cod_rating)    

    cod_seniority = 4
    
    cod_rr   = getCodiceRecovery(cod_seniority, cod_sector, bond_date)    
    value_rr = getRecoveryValue(cod_rr)



