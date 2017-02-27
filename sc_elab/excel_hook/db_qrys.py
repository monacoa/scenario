
from connection import *

def getCurvesListFromDb(curve_date):
    c_date_qry = str(curve_date).replace("-", "")
    con = Connection()
    qry = '''
          select distinct BloombergTicker FROM alm.DProTS_master where BloombergTicker in
          ( select distinct BloombergTicker from alm.DProCurve
            where TipoDato in ('CDepositi', 'CSwap', 'CLibor', 'CFuture')
          )and Data = '%s'
          ''' % c_date_qry
    print qry
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
    print "qry:", qry
    c_anag.execute(qry)
    res = c_anag.fetchall()
    print res
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


def getDatesListFromDb():
    cn  = Connection()
    cur = cn.db_data()
    qry = '''
                     select distinct Data from alm.DProTS_master where BloombergTicker in
                     ( select distinct BloombergTicker from alm.DProCurve where TipoDato in
                         ('CDepositi', 'CSwap', 'CLibor', 'CFuture')
                      )
                     order by Data desc
                    '''
    cur.execute(qry)
    res         = cur.fetchall()
    date_list   = []
    for record in res:
        date_list.append(record[0])
    return date_list