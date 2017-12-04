
from connection import *

def getCurvesListFromDb(curve_date, type):

    if type == "SWP":
        table        = "DProCurve"
        where_clause = " where TipoDato in ('CDepositi', 'CSwap', 'CLibor', 'CFuture') "
    elif type == "CDS":
        table = "DProCDS"
        where_clause = ""

    c_date_qry = str(curve_date).replace("-", "")
    con = Connection()
    qry = '''
          select distinct BloombergTicker FROM alm.DProTS_master where BloombergTicker in
          ( select distinct BloombergTicker from alm.%s %s
          )and Data = '%s'
          ''' %(table, where_clause, c_date_qry)
    #print qry
    c_data = con.db_data()
    c_data.execute(qry)
    tickers = c_data.fetchall()
    tkrs = ""
    sep  = ""
    for tic in tickers:
        tkrs += sep + "'" + tic[0] + "'"
        sep = ","
    #================
    c_anag = con.db_anag()
    qry = '''
          SELECT DISTINCT description, ID_object FROM AnagMktObject WHERE ID_object in
          ( SELECT DISTINCT ID_object FROM AnagMktObject_D WHERE ticker IN
            (%s
            )
          )ORDER BY description
          ''' % tkrs
    #print "qry:", qry
    c_anag.execute(qry)
    res = c_anag.fetchall()
    #print res
    ll = []
    for c in res:
        #print "c", c
        curva = str(c[0])
        id    = int(c[1])
        a =[curva,id]
        ll.append(a)
        #print "ll", ll
    con.close()
    return ll
    # -----------

def getBondListFromDb(bond_date, type):

    c_date_qry = str(bond_date).replace("-", "")
    con = Connection()
    
    
    c_anag = con.db_data()
    

    
    qry = '''
          SELECT DISTINCT Descrizione FROM view_Bond_master  WHERE Data = '%s' 
          ''' %(c_date_qry)
          
          
    #print "qry:", qry
    c_anag.execute(qry)
    res = c_anag.fetchall()
    #print 'risultato qry: ', res
    ll = []
    k = 0
    for c in res:
        #print "c", c
        bond_des = str(c[0])
        id    = k + 1
        a =[bond_des,id]
        
        ll.append(a)
        #print "ll", ll
    con.close()
    return ll
    # -----------



def getDatesListFromDb(type):
    cn  = Connection()
    cur = cn.db_data()
    if type == "SWP":

        qry = '''
                     select distinct Data from alm.DProTS_master where BloombergTicker in
                     ( select distinct BloombergTicker from alm.DProCurve where TipoDato in
                         ('CDepositi', 'CSwap', 'CLibor', 'CFuture')
                      )
                     order by Data desc
                    '''
    elif type == "CDS":
        qry = '''
                    select distinct Data from alm.DProTS_master where BloombergTicker in
                    (
                        select distinct BloombergTicker from alm.DProCDS
                    )
                    order by Data desc
                    '''

    elif type == "BOND":
        qry = '''
                    select distinct Data from view_Bond_master where EMITTENTE 
                    is not Null order by Data desc
                    '''



    cur.execute(qry)
    res         = cur.fetchall()
    date_list   = []
    for record in res:
        date_list.append(record[0])
    cn.close()
    return date_list

def getProvidersFromDb (table):
    Con = Connection()
    c   = Con.db_anag()
    qry = '''
            select distinct FONTE from %s
            '''%table
    c.execute(qry)
    res=c.fetchall()
    #print res
    sects = []
    for record in res:
        sects.append(record[0])
    Con.close()
    return sects

def getCdsCodeFromDb(cds_curve):
    con = Connection()
    c = con.db_ang()
    qry ='''
            select codice_curva from MKT_Curve where
            tipo_curva = 'CDS' and
            emittente in (select codice_emittente from ZEmittente where descrizione = '%s') and
            seniority in (select codice_seniority from ZSeniority where descrizione = '%s') and
    '''