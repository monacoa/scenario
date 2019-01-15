import pandas as pd
from Tkinter import *

import tkMessageBox
#import mysql.connector
import pyodbc

from sc_elab.core.db_tipodato_censimento import tipologia_dict
from sc_elab.core.table_traslation import from_bond_table_to_datastream
from sc_elab.core.db_data_structure_v0 import table_dict
from sc_elab.core.table_traslation import min_field_set

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

    list_f_table = table_dict['bond_master'].keys()

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

    for i in range(0, n_records):

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
            data_dict_to_insert_tmp[f][0] = data_anag[f][i]

        tickerTmp = data_anag[field_ref][i]

        value_ref = tickerTmp

        print 'record n. %s' % i
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

    db.commit()
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


def TEST_insert_data_record(connection_status, table_data, table_anag, data_values):

    cursor = connection_status['cursor']
    db = connection_status['db_connection']
    db.autocommit = False

    field_ref = 'BloombergTicker'  # pk della tabella
    ticker_list = data_values['Ticker']
    n_records = len(ticker_list)

    msg = 'Inserimento andato a buon fine'
    result_val = True
    verboseFlag = False

    for i in range(0, n_records):

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

        tickerTmp = data_values['Ticker'][i]

        value_ref = tickerTmp

        print 'record n. %s' % i
        result_anag, anag_dict_res = retrive_record_from_table(cursor, table_anag, field_ref, value_ref)

        if (result_anag == False):
            msg = 'Ticker %s non censito in %s!!' % (tickerTmp, table_anag)
            result_val = False

            return result_val, msg


        data_dict_to_insert_tmp['Contributor'][0] = anag_dict_res['Contributor']
        data_dict_to_insert_tmp['Datatype'][0] = 'Livello'
        data_dict_to_insert_tmp['ValoreBid'][0] = 0.0
        data_dict_to_insert_tmp['ValoreAsk'][0] = 0.0
        data_dict_to_insert_tmp['TipoDato'][0] = anag_dict_res['TipoDato']
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
    db.close()

    return result_val, msg




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

from sc_elab.excel_hook.connection import Connection

#def test_load_nuovi_dati(file_new_data='C:/Users/scalambrinm/workspace/scenario/sc_elab/core/input/files_caricamento_datastream/Curva_Depositi_EUR_ICAP.xlsx'):
def test_load_nuovi_dati(file_new_data):


        # SETUP CREDENZIALI
     #   db_credential = {}
     #   #db_credential['host'] = "10.103.65.195"
     #   #db_credential['user'] = "xxx3"
     #   #db_credential['pwd'] = "yyy2"
#
#
     #   db_credential['host'] = "localhost"
     #   db_credential['user'] = "root"
     #   db_credential['pwd'] = "DatabaseRepl1ca"
     #   db_name = "db_mercato_matteo"
#
     #   db = pyodbc.connect(r'DSN=db_mercato_matteo;UID=lucap;PWD=lucap')
     #   cursor = db.cursor()

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
        dati_exist = False

        for check_foglio_dati in lista_fogli:
            test3 = check_foglio_dati.find(u'Dati')
            if not(test3==-1):
                dati_exist= True
                break

        if dati_exist == False:
            root = Tk()
            msg_errore = 'Non esiste alcun foglio di anagrafica'
            tkMessageBox.showinfo("Chk", msg_errore)
            root.destroy()
            return

        #COMINCIO A LEGGERE I VALORI PRESENTI ALL'INTERNO DELLO SCARICO
        df1 = pd.read_excel(file_new_data, sheet_name='Dati')
        #lista_fogli= dftest.keys()
        # leggo il contenuto del foglio Dati
        data_dict1 = df1.to_dict()
        data_values = data_dict1
        #controllo l'intestazione
        field_ref_list = ['Data', 'Valore', 'Ticker']
        #controllo che esistano effettivamente dei valori all'interno del foglio 'Dati'
        uscita0 = False
        uscita1= False
        uscita2 = False
        num_isin = 0
        num_valori = 0
        try:
            lista_date= data_values ['Data'] # IL CAMPO DATA E' L'UNICO CHE ACCOMUNA IL TEMPLATE DEI BOND CON GLI ALTRI TEMPLATE
            num_date = len(lista_date)
        except:
            uscita0 =True

        try:
            lista_valori = data_values ['Ticker'] # CHIAVE DI TUTTE LE TABELLE TRANNE I BOND
            num_valori = len (lista_valori)
        except:
            uscita1= True

        try:
            lista_valori = data_values[u'ISIN'] #CHIAVE DEI BOND
            num_isin = len (lista_valori)
        except:
            uscita2= True

        if ((uscita1 == True and uscita2 == True) or uscita0==True) :
            root = Tk()
            msg_alert= 'Attenzione! Manca il codice identificativo dello strumento'
            tkMessageBox.showinfo("Chk", msg_alert)
            root.destroy()
        #elif num_valori == 0:
        #        msg_alert ='Attenzione! Nel foglio Dati non esiste alcun valore'
        #elif num_date ==0:
        #    msg_alert = 'Attenzione! Nel foglio Dati non vengono riportate le date'


        #DOPPIO CHECK: CONTROLLO CHE IL BLOOMBERG TICKER SIA CENSITO ED IN CASO DOVE SIA CENSITO
        find_table = False # BOOLEANA CHE MI INDICA SE IL FILE E' CENSITO O MENO
        # ------- NOTA ==> DA PERSONALIZZARE PER CARICAMENTO MASSIMO
        if num_isin != 0 :
            table_name_anag = 'bond_master'
            find_table = True
        else:
            elenco_escluso_bond = elenco_tabelle.copy()

            del elenco_escluso_bond['bond_master']

            # QUESTO CICLO E' CREATO PER TESTARE LA PRESENZA DI TUTTI I TICKER ALL'INTERNO DEI DATABASE DI ANAGRAFICA
            # AL MOMENTO CICLA SU TUTTI I BLOOMBERG TICKER ANCHE SE IN REALTA' RESTITUISCE UNA SOLA STRINGA,
            # PERO' QUANDO CI SARA' UN CARICAMENTO MASSIMO IN CUI I TICKER POSSONO ARRIVARE DA TABELLE DIVERSE
            # IL TABLE_NAME_ANAG DIVENTERA' UN VETTORE

            for j in range(0, num_valori):
                find_table = False
                for table_name in (elenco_escluso_bond):
                        field_ref = 'BloombergTicker'
                        tickerTmp = data_values[u'Ticker'][j]
                        result_ticker, table_name_test = retrive_table_from_ticker(cursor , table_name, tickerTmp, field_ref)
                        if result_ticker == True:
                            table_name_anag = table_name
                            find_table =True

                if find_table == False:
                    tickerMissing = tickerTmp


        # SE IL DATO NON E' CENSITO ALLORA CERCO DI CENSIRLO
        if find_table == False:

            test_campi = False
            #test5=lista_fogli.find(u'Anagrafica')
            #CERCO LA PRESENZA DEL FOGLIO ANAGRAFICA
            for foglio in lista_fogli:
                test2 = foglio.find('Anagrafica')
                if not (test2==-1):
                    nome_tabella= foglio
                    #lunghezza_nome_foglio = len(foglio)
                    #if lunghezza_nome_foglio > 10:
                    #    nome_tabella= foglio.replace('Anagrafica'+'_','')
                    #else:
                    #   nome_tabella = 'Anagrafica'
                    test4 =True
                    break
            if test2 == -1:
                root = Tk()
                msg_alert = 'Non esiste alcun foglio di anagrafica'
                tkMessageBox.showinfo("Chk", msg_alert)
                root.destroy()
                return


            #QUI CERCO DI CAPIRE DOVE VA IL DATO DA ANAGRAFARE IN FUNZIONE DELLA PRIMA RIGA DEL FOGLIO ANAGRAFICA.
            #IN PARTICOLARE CONFRONTO IL FOGLIO ANAGRAFICA CON LA LISTA DI TUTTI I CAMPI DELLE VARIE TABELLE
            #SE NESSUNO COINCIDE, ALLORA IL DATO NON PUO' EFFETTIVAMENTE ESSERE INSERITO IN NESSUNA TABELLA


            anagrafica_tabelle = table_dict.copy()
            df2 = pd.read_excel(file_new_data, sheet_name=nome_tabella)
            data_dict2 = df2.to_dict() # LA PRIMA RIGA DEL FOGLIO ANAGRAFICA VIENE LETTA COME UN INSIEME DI CHIAVI
            nome_campi = data_dict2.keys()
            for tabella_appoggio in anagrafica_tabelle:
                lista_campi=anagrafica_tabelle[tabella_appoggio].keys()
                lista_campi.sort()
                nome_campi.sort()
                if lista_campi == nome_campi:    #SE SODDISFATTA VIENE SODDISFATTA, SIGNIFICA CHE HO TROVATO IL NOME DELLA TABELLA
                    test_campi = True
                    nome_tabella =  tabella_appoggio
                    break

            if test_campi == False:
                root = Tk()
                msg_alert = 'Impossibile censire. I campi da anagrafare non corrispondono a quelli del database'
                tkMessageBox.showinfo("Chk", msg_alert)
                root.destroy()
                return

            #ORA POSSIEDO TUTTE LE INFORMAZIONI PER CENSIRE L'ANAGRAFICA
            table_name_anag = nome_tabella
            df2 = pd.read_excel(file_new_data, sheet_name='Anagrafica')
            data_dict2 = df2.to_dict()
            data_anag = data_dict2
            val_results, msg_anag = insert_anag_record(connection_status, table_name_anag, data_anag)
            ################################if

            if val_results == False:
                root = Tk()
                tkMessageBox.showinfo("Errore", msg_anag)
                root.destroy()
                return

        #SE NON SONO USCITO MAI FINORA SIGNIFICA CHE POSSO PROCEDERE ALL'INSERIMENTO DEL DATO

        if table_name_anag== 'DProCDS':
            table_name_data = 'DProTS_master'
            field_ref_list = ['Data', 'Valore', 'Ticker']
            val_results0, msg_field = chk_data_fields(field_ref_list, data_dict1)
            # tkMessageBox.showinfo("Chk", msg_field)
            print 'try insert data'
            df2 = pd.read_excel(file_new_data, sheetname='Anagrafica')
            data_dict2 = df2.to_dict()
            data_anag = data_dict2
            val_results, msg_insert = TEST_insert_cds_data_record(connection_status, table_name_data,table_name_anag, data_values)
            if val_results == False:
                root = Tk()
                tkMessageBox.showwarning("Chk", msg_insert)
                root.destroy()
                return
            else:
                root = Tk()
                tkMessageBox.showinfo("Work done", "inserimento andato a buon fine")
                root.destroy()
                return



        elif (table_name_anag == 'bond_master'):
            """
            print 'chk field'
            field_ref_list     = ['Data', 'Valore', 'Ticker']
            val_results0, msg_field = dbt.chk_data_fields(field_ref_list, data_dict1)
            tkMessageBox.showinfo("Chk", msg_field)    
            """

            table_name_data = 'bond_master'

            print 'try insert data'
            val_results, msg_insert = TEST_insert_bond_data_record(connection_status, table_name_data, data_values)
            if val_results == False:
                root = Tk()
                tkMessageBox.showwarning("Chk", msg_insert)
                root.destroy()
                return
            else:
                root = Tk()
                tkMessageBox.showinfo("Work done", "inserimento andato a buon fine")
                root.destroy()
                return

        else:
            print 'try insert data'
            table_name_data = 'DProTS_master'
            val_results, msg_insert = TEST_insert_data_record(connection_status, table_name_data, table_name_anag,
                                                             data_values)
            if val_results == False:
                root=Tk()
                tkMessageBox.showwarning("Chk", msg_insert)
                root.destroy()
                return
            else:
                root = Tk()
                tkMessageBox.showinfo("Work done","inserimento andato a buon fine")
                root.destroy()
                return



#if __name__ == "__main__":
#    file_new_data = 'C:/Users/scalambrinm/workspace/scenario/sc_elab/core/input/files_caricamento_datastream/Curva_Depositi_EUR_ICAP.xlsx'
#    test_load_nuovi_dati(file_new_data)