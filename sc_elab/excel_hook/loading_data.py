import pandas as pd
from Tkinter import *

import tkMessageBox
#import mysql.connector
import pyodbc

from sc_elab.core.db_tipodato_censimento import tipologia_dict
from sc_elab.core.table_traslation import from_bond_table_to_datastream
from sc_elab.core.db_data_structure_v0 import table_dict
from sc_elab.core.table_traslation import min_field_set
from sc_elab.core.Tipologia_curva_dizionario import corrisp_tabella

#def FQ(label):
#    print ('------------- FIN QUI TUTTO OK  %s ----------' %(label))
#    sys.exit()


def TEST_insert_bond_data_record(connection_status, table_data, data_values):

    cursor = connection_status['cursor']
    db = connection_status['db_connection']
    db.autocommit = False

    field_ref = 'ISIN'  # pk della tabella
    isin_list = data_values[field_ref]
    n_records = len(isin_list)

    msg = 'Inserimento andato a buon fine'
    result_val = True
    verboseFlag = False

    list_f_table = table_dict['Bond_master'].keys()

    print 'data_values: ', data_values.keys()

    for i in range(0, n_records):

        data_dict_to_insert_tmp = {}

        for f_tab in list_f_table:

            f_dt = from_bond_table_to_datastream[f_tab]
            data_dict_to_insert_tmp[f_tab] = {}

            if (f_dt == 'NULL'):
                data_dict_to_insert_tmp[f_tab][0] = 'nan'
                pass
            else:
                data_dict_to_insert_tmp[f_tab][0] = data_values[f_dt][i]
#                print f_dt + 'value is: ' + str(data_values[f_dt][i])

        isinTmp = data_dict_to_insert_tmp['isin'][0]

        try:
            result_data = insert_data_on_table(cursor, table_data, data_dict_to_insert_tmp,verboseFlag)

        except Exception as e:

            msg = 'Inserimento ticker %s in %s non riuscito: %s!!' % (isinTmp, table_data, e)
            result_val = False
            db.close()

            return result_val, msg

    db.commit()
    db.close()

    return result_val, msg


def TEST_insert_cds_data_record(connection_status, table_data, table_anag, data_values):

    cursor = connection_status['cursor']
    db = connection_status['db_connection']
    db.autocommit = False

    field_ref = 'BloombergTicker'  # pk della tabella
    ticker_list = data_values['Ticker']
    n_records = len(ticker_list)

    msg = 'Inserimento andato a buon fine'
    result_val = True
    verboseFlag = False

    for i in xrange(0, n_records):

        data_dict_to_insert_tmp = {}

        for f_tab in min_field_set:

            data_dict_to_insert_tmp[f_tab] = {}

        tickerTmp = data_values['Ticker'][i]

        value_ref = tickerTmp

        print 'record n. %s' % i
        result_anag, anag_dict_res = retrive_record_from_table(cursor, table_anag, field_ref, value_ref)

        if (result_anag == False):
            msg = 'Ticker %s non censito in %s!!' % (tickerTmp, table_anag)
            result_val = False

            return result_val, msg

        # print 'anag_dict_res[TipoTicker]: ', anag_dict_res['TipoTicker']
        # print 'data_values[Ticker][i]', data_values['Ticker'][i]

        data_dict_to_insert_tmp['Contributor'][0] = anag_dict_res['Contributor']
        data_dict_to_insert_tmp['Datatype'][0] = 'Livello'
        data_dict_to_insert_tmp['ValoreBid'][0] = 0.0
        data_dict_to_insert_tmp['ValoreAsk'][0] = 0.0
        data_dict_to_insert_tmp['TipoDato'][0] = 'CDS'
        data_dict_to_insert_tmp['id'][0] = data_values['Ticker'][i] + anag_dict_res['TipoTicker']

        data_dict_to_insert_tmp['ValoreMid'][0] = data_values['Valore'][i]
        data_dict_to_insert_tmp['LastUpdate'][0] = data_values['Data'][i]
        data_dict_to_insert_tmp['DataScarico'][0] = data_values['Data'][i]
        data_dict_to_insert_tmp['Data'][0] = data_values['Data'][i]
        data_dict_to_insert_tmp['BloombergTicker'][0] = data_values['Ticker'][i]

        try:
            result_data = insert_data_on_table(cursor, table_data, data_dict_to_insert_tmp, verboseFlag)

        except Exception as e:

            msg = 'Inserimento ticker %s in %s non riuscito: %s!!' % (tickerTmp, table_data, e)
            result_val = False
            db.close()

            return result_val, msg

    db.commit()
    #db.close()

    return result_val, msg


def insert_anag_record(connection_status, table_anag, data_anag):

    cursor = connection_status['cursor']
    db = connection_status['db_connection']
    db.autocommit = False

    field_ref = 'BloombergTicker'
    ticker_list = data_anag[field_ref]
    n_records = len(ticker_list)

    msg = 'Inserimento andato a buon fine'
    result_val = True

    verboseFlag = False
    field_list = data_anag.keys()
    data_dict_to_insert_tmp = {}
    for f in field_list:
        data_dict_to_insert_tmp[f] = {}


    for i in range(0, n_records):

        for f in field_list:
            #data_dict_to_insert_tmp[f][0] = data_anag[f][i]
            data_dict_to_insert_tmp[f][0] = data_anag[f][i]
        #tickerTmp = data_anag[field_ref][i]
        tickerTmp = data_anag[field_ref][i]

        value_ref = tickerTmp

        print 'Anagrafica inserimento - record n. %s' %tickerTmp
        result_anag, anag_dict_res = retrive_record_from_table(cursor, table_anag, field_ref, value_ref)

        if (result_anag == True):
            msg = 'Ticker %s gia censito in %s!!' % (tickerTmp, table_anag)
            result_val = False

            return result_val, msg

        try:
            result_data = insert_data_on_table(cursor, table_anag, data_dict_to_insert_tmp,verboseFlag)

        # except:
        except Exception as e:

            msg = 'Inserimento ticker %s in %s non riuscito: %s!!' % (tickerTmp, table_anag, e)
            result_val = False

            return result_val, msg

    #db.commit()
    #db.close()

    return result_val, msg


def chk_data_fields(field_ref_list, data_dict1):
    filed_list = data_dict1.keys()

    val_out = True
    msg = 'Nome dei campi del file di caricamento corretti!'

    for f in filed_list:

        if (f not in field_ref_list):
            msg = 'Il campo %s non e presente nel file di input!!' % f
            val_out = False

            return val_out, msg

    return val_out, msg


def insert_data_on_table(cursor, table_name, data_dict, verboseFlag):
    # data_dict [NOME CAMPO TABELLA][N. RECORD] = VALORE

    field_list = data_dict.keys()
    f0 = field_list[0]
    n_fields = len(field_list)

    c_str = "'"
    v_str = " , "
    n_records = len(data_dict[f0].keys())

    if (verboseFlag == True): print 'n_records: ', n_records

    if (n_records > 50000):
        unit_ref = 10000
    elif (n_records > 30000):
        unit_ref = 5000
    elif (n_records > 10000):
        unit_ref = 2000
    elif (n_records > 5000):
        unit_ref = 1000
    elif (n_records > 1000):
        unit_ref = 500
    elif (n_records > 500):
        unit_ref = 200
    else:
        unit_ref = 100

    for i in range(0, n_records):

        valList = '('
        fList = '('

        if (verboseFlag == True):
            c = divmod(i, unit_ref)

            if (c[1] == 0) and (i != 0):

                print 'Caricati i primi %s records' % i

            elif (c[1] == 0) and (i == 0):

                print 'Caricamento dei primi %s records' % unit_ref

        first_time = True
        for j in range(0, n_fields):

            fTmp = field_list[j]
            if j==10:
                app3='here'
            fTmp_str = str(fTmp)

            valTmp = data_dict[fTmp][i]

            try:
                valTmp_str = str(valTmp)
            except:
                valTmp_str = valTmp.encode('utf-8')

            if (valTmp_str != 'nan') and (valTmp_str != 'NaT'):

                if (first_time == True):

                    valList = valList + c_str + valTmp_str + c_str
                    fList = fList + fTmp_str

                    first_time = False
                else:

                    valList = valList + v_str + c_str + valTmp_str + c_str
                    fList = fList + v_str + fTmp_str
            else:

                continue

        valList = valList + ')'
        fList = fList + ')'

        qry_to_execute = 'insert into ' + table_name + ' ' + fList + ' ' + 'values' + valList

        cursor.execute(qry_to_execute)

    return True

def retrive_record_from_table(cursor, table_name, field_ref, ticker_tmp):

    qry_to_execute = "SELECT * FROM %s WHERE %s = '%s' " % (table_name, field_ref, ticker_tmp)

    cursor.execute(qry_to_execute)
    result = cursor.fetchall()

    cursor.execute("SHOW columns FROM %s" % (table_name))
    result_names = cursor.fetchall()

    dict_results = {}

    if len(result) > 0:

        result_flag = True

        for i in range(0, len(result[0])):
            fTmp = result_names[i][0]
            vTmp = result[0][i]

            dict_results[fTmp] = vTmp

    else:
        result_flag = False

    return result_flag, dict_results

def retrive_record_from_table_check_anag(cursor, table_name, field_ref, ticker_tmp):

    qry_to_execute = "SELECT * FROM %s WHERE %s = '%s' " % (table_name, field_ref, ticker_tmp)

    cursor.execute(qry_to_execute)
    result = cursor.fetchall()

    cursor.execute("SHOW columns FROM %s" % (table_name))
    result_names = cursor.fetchall()

    dict_results = {}

    if len(result) > 0:

        result_flag = True

        for i in range(0, len(result[0])):
            fTmp = result_names[i][0]
            vTmp = result[0][i]

            dict_results[fTmp] = [vTmp]

    else:
        result_flag = False

    return result_flag, dict_results



def retrive_table_from_ticker(cursor, table_name, tickerTmp, field_ref):

    #QUESTA QUERY CERCA PER OGNI TABELLA

    table_name_anag =""
    result_flag = False

    qry_to_execute = "SELECT * FROM %s WHERE %s = '%s' " % (table_name, field_ref, tickerTmp)

    cursor.execute(qry_to_execute)
    result = cursor.fetchall()
    if len(result) > 0:
        result_flag = True
    else:
        result_flag = False

    return result_flag, result


def Insert_data_record(connection_status, table_data, anag_dict_res, data_values):

    cursor = connection_status['cursor']
    db = connection_status['db_connection']
    db.autocommit = False

    field_ref = 'BloombergTicker'  # pk della tabella
    ticker_list = data_values['Ticker']
    n_records = len(ticker_list)

    msg = 'Inserimento andato a buon fine'
    result_val = True
    verboseFlag = False

    data_dict_to_insert_tmp = {}
    data_dict_to_insert_tmp['BloombergTicker'] = {}
    data_dict_to_insert_tmp['Contributor'] = {}
    data_dict_to_insert_tmp['Data'] = {}
    data_dict_to_insert_tmp['Datatype'] = {}
    data_dict_to_insert_tmp['ValoreBid'] = {}
    data_dict_to_insert_tmp['ValoreAsk'] = {}
    data_dict_to_insert_tmp['ValoreMid'] = {}
    data_dict_to_insert_tmp['LastUpdate'] = {}
    data_dict_to_insert_tmp['DataScarico'] = {}
    data_dict_to_insert_tmp['TipoDato'] = {}
    data_dict_to_insert_tmp['id'] = {}

    tickerTmp = data_values['Ticker'][0]

    data_dict_to_insert_tmp['Contributor'][0] = anag_dict_res['Contributor'][0]
    data_dict_to_insert_tmp['Datatype'][0] = 'Livello'
    data_dict_to_insert_tmp['ValoreBid'][0] = 0.0
    data_dict_to_insert_tmp['ValoreAsk'][0] = 0.0

    if anag_dict_res['tableDB'] == 'DProCDS':
        data_dict_to_insert_tmp['TipoDato'][0] = 'CDS'
    else:
        data_dict_to_insert_tmp['TipoDato'][0] = anag_dict_res['TipoDato'][0]

    data_dict_to_insert_tmp['id'][0] = data_values['Ticker'][0] + anag_dict_res['TipoTicker'][0]

    data_dict_to_insert_tmp['ValoreMid'][0] = data_values['Valore'][0]
    data_dict_to_insert_tmp['LastUpdate'][0] = data_values['Data'][0]
    data_dict_to_insert_tmp['DataScarico'][0] = data_values['Data'][0]
    data_dict_to_insert_tmp['Data'][0] = data_values['Data'][0]
    data_dict_to_insert_tmp['BloombergTicker'][0] = data_values['Ticker'][0]

    try:
        result_data = insert_data_on_table(cursor, table_data, data_dict_to_insert_tmp, verboseFlag)

    except Exception as e:

        msg = 'Inserimento ticker %s in %s non riuscito: %s!!' % (tickerTmp, table_data, e)
        result_val = False
        db.close()

        return result_val, msg

    #db.commit()
    #db.close()

    return result_val, msg





from sc_elab.excel_hook.connection import Connection


def load_all_new_anag(connection_status,dict):

    for tkt in dict.keys():
        data_anag = dict[tkt].copy()
        table_name_anag = data_anag.pop('tableDB')
        val_results, msg_anag = insert_anag_record(connection_status, table_name_anag, data_anag)

        if val_results == False:
            root = Tk()
            tkMessageBox.showinfo("Errore", msg_anag)
            root.destroy()
            return

def take_anagrafica(ticker, presenti,assenti):

    if ticker in presenti.keys():
        out = presenti[ticker].copy()
        tmp_out = out.pop(u'insertdate')
    else:
        out = assenti[ticker].copy()

    return out

def close_loading(db,msg,header):
    #db.close()
    root = Tk()
    tkMessageBox.showinfo(header, msg)
    root.destroy()


#def test_load_nuovi_dati(file_new_data='C:/Users/scalambrinm/workspace/scenario/sc_elab/core/input/files_caricamento_datastream/test_gennaio.xlsx'):
def test_load_nuovi_dati(file_new_data):


        con = Connection()
        cursor = con.db_data()

        connection_status = {}
        connection_status['type'] = True
        connection_status['cursor'] = cursor
        connection_status['db_connection'] = con.db

        # CARICAMENTO DATI
        #file_new_data = 'input/files_caricamento_datastream/Curva_Depositi_EUR_ICAP.xlsx'

        elenco_tabelle = tipologia_dict
        tabelle_values = elenco_tabelle
        num_tabelle = len(tabelle_values)

        #PRIMA COSA: CONTROLLO SE ESISTE UN FOGLIO CHIAMATO "FOGLIO"

        dftest = pd.read_excel(file_new_data, None)
        lista_fogli = dftest.keys()

        if not('Dati' in lista_fogli):
            msg_errore = 'Non esiste alcun foglio di inserimento dati'
            close_loading(con.db,msg_errore,"Chk")
            return

        anag_wb = {}
        for check_foglio_dati in lista_fogli:
            tmp = check_foglio_dati.split('Anagrafica_')
            if len(tmp) > 1:
                anag_wb[tmp[1]] = dftest[check_foglio_dati]
                ##### VERIFICA BONTA ANAGRAFICA
                ##### INSERIRE QUI
                ##### corrisp_tabella

        #COMINCIO A LEGGERE I VALORI PRESENTI ALL'INTERNO DELLO SCARICO
        df1 = dftest[u'Dati']
        data_values = dftest[u'Dati'].to_dict()
        #controllo l'intestazione
        #field_ref_list = ['Data', 'Valore', 'Ticker']
        #controllo che esistano effettivamente dei valori all'interno del foglio 'Dati'
        uscita0 = False
        uscita1 = False
        uscita2 = False

        num_isin   = 0
        num_valori = 0

        #### controllo il foglio Dati
        if 'Data' in data_values.keys():
            lista_date= data_values['Data'] # IL CAMPO DATA E' L'UNICO CHE ACCOMUNA IL TEMPLATE DEI BOND CON GLI ALTRI TEMPLATE
            num_date = len(lista_date)
        else:
            uscita0 =True

        if 'Ticker' in data_values.keys():
            lista_valori = data_values ['Ticker'] # CHIAVE DI TUTTE LE TABELLE TRANNE I BOND
            num_valori = len (lista_valori)
        else:
            uscita1= True

        if 'ISIN' in data_values.keys():
            lista_valori = data_values[u'ISIN'] #CHIAVE DEI BOND
            num_isin = len (lista_valori)
        else:
            uscita2= True

        if ((uscita1 == True and uscita2 == True) or uscita0==True) :
            msg_alert= 'Attenzione! Manca il codice identificativo dello strumento'
            close_loading(con.db,msg_alert,"Chk")
            return

        # ------- NOTA ==> DA PERSONALIZZARE PER CARICAMENTO MASSIMO
        if num_isin == 0 :
            elenco_escluso_bond = elenco_tabelle.copy()
            del elenco_escluso_bond['bond_master']

            anag_assenti = {}
            anag_presenti = {}

            # ciclo la tabella Dati per verificare la presenza o meno dell'anagrafica
            for j in xrange(0, num_valori):
                find_table = False
                field_ref = 'BloombergTicker'
                tickerTmp = data_values[u'Ticker'][j]

                # verifico che sia gia' presente all'interno del DB
                if not (tickerTmp in anag_presenti.keys()):
                    for table_name in (elenco_escluso_bond):
                        result_ticker, table_name_test = retrive_record_from_table_check_anag(cursor , table_name, field_ref, tickerTmp)
                        if result_ticker == True:
                            table_name_anag = table_name
                            find_table = True
                            ### VERIFICARE COME FATTO table_name_test
                            anag_presenti[tickerTmp] = table_name_test
                            anag_presenti[tickerTmp][u'tableDB'] = table_name

                    # se non e' presente verifico che sia presente all'interno del foglio Excel
                    if find_table == False:
                        if not (tickerTmp in anag_assenti.keys()):
                            tickerMissing = tickerTmp
                            flag_no_anag = False
                            for sheet in anag_wb.keys():
                                df2 = anag_wb[sheet]
                                # controllo se tutti i valori della colonna sono nulli
                                if df2[field_ref].isnull().all():
                                    index = [False]
                                else:
                                    index = df2[field_ref] == tickerMissing

                                if sum(index) == 1:
                                    flag_no_anag = True
                                    anag_assenti[tickerMissing] = df2.loc[index, :].to_dict('list')

                                    if sheet in corrisp_tabella.keys():
                                        anag_assenti[tickerMissing][u'tableDB'] = corrisp_tabella[sheet]
                                    else:
                                        msg_alert = 'Il foglio %s non ha alcuna associazione con una tabella del database' % (sheet)
                                        close_loading(con.db, msg_alert, "Chk")
                                        return
                                elif sum(index) > 1:
                                    msg_errore = "Il Ticker %s ha un'anagrafica duplicata. Ricontrolla" % (
                                        tickerMissing)
                                    close_loading(con.db, msg_errore, "Chk")
                                    return

                            if flag_no_anag == False:
                                msg_errore = 'Il Ticker %s nel foglio dati non esiste in alcuna anagrafica. Ricontrolla'%(tickerMissing)
                                close_loading(con.db, msg_errore, "Chk")
                                return

                        else:
                            msg_errore = 'Il Ticker %s duplicato nel foglio Dati' %(tickerTmp)
                            close_loading(con.db, msg_errore, "Chk")
                            return

            # carico tutte le anagrafiche mancanti
            load_all_new_anag(connection_status,anag_assenti)

            for j in xrange(0, num_valori):
                tickerTmp  = data_values[u'Ticker'][j]
                tickerAnag = take_anagrafica(tickerTmp,anag_presenti,anag_assenti)

                table_name_data = 'DProTS_master'
                data_tick = df1.loc[df1['Ticker'] == tickerTmp,:].to_dict('list')

                print 'try insert data: '+ str(tickerTmp)
                val_results, msg_insert = Insert_data_record(connection_status, table_name_data, tickerAnag,
                                                             data_tick)
                if val_results == False:
                    close_loading(con.db, msg_insert, "Chk")
                    return

            con.db.commit()
            close_loading(con.db, "Inserimento andato a buon fine", "Work done")
            return



        else:
            find_table = True
            table_name_data = 'Bond_master'

            print 'try insert data'
            val_results, msg_insert = TEST_insert_bond_data_record(connection_status, table_name_data, data_values)
            if val_results == False:
                close_loading(con.db, msg_insert, "Chk")
                return

                #root = Tk()
                #tkMessageBox.showwarning("Chk", msg_insert)
                #root.destroy()
                #return
            else:
                close_loading(con.db, "Inserimento andato a buon fine", "Work done")
                return


#if __name__ == "__main__":
#    file_new_data = 'C:/Users/scalambrinm/workspace/scenario/sc_elab/core/input/files_caricamento_datastream/test_gennaio.xlsx'
#    test_load_nuovi_dati(file_new_data)