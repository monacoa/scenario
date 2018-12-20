



# ------------ DProTSmaster ---------------------------------
table_dict_1 = {}
#table_dict_1['InsertDate']       = ['TYPE', 'PRIMARY KEY', FOREIGN KEY, REFERENCES]

table_dict_1['BloombergTicker']   = ['varchar(50)',    'True', 'False', 'False']
table_dict_1['Contributor']       = ['varchar(20)',   'False', 'False', 'False']
table_dict_1['Data']              = ['DATE',           'True', 'False', 'False']
table_dict_1['Datatype']          = ['varchar(20)',   'False', 'False', 'False']
table_dict_1['ValoreBid']         = ['double',        'False', 'False', 'False']
table_dict_1['ValoreMid']         = ['double',        'False', 'False', 'False']
table_dict_1['ValoreAsk']         = ['double',        'False', 'False', 'False']
table_dict_1['ScadenzaFuture']    = ['DATE',          'False', 'False', 'False']
table_dict_1['DaysToExpiration']  = ['INT(4)',        'False', 'False', 'False']
table_dict_1['LastUpdate']        = ['DATE',          'False', 'False', 'False']
table_dict_1['DataScarico']       = ['DATE',          'False', 'False', 'False']
table_dict_1['id']                = ['varchar(50)',   'False', 'False', 'False']
table_dict_1['TipoDato']          = ['varchar(50)',   'False', 'False', 'False']
table_dict_1['Volume']            = ['double',        'False', 'False', 'False']


# ------------ DProCurve ---------------------------------------------------------
table_dict_2 = {}
table_dict_2['Rating']              = ['INT(4)',        'False', 'False', 'False']
table_dict_2['Currency']            = ['varchar(20)',   'False', 'False', 'False']
table_dict_2['MaturityInt']         = ['INT(4)',        'False', 'False', 'False']
table_dict_2['TipoTicker']          = ['varchar(10)',   'False', 'False', 'False']
table_dict_2['SettorePaese']        = ['INT(4)',        'False', 'False', 'False']
table_dict_2['BloombergTicker']     = ['varchar(50)',    'True', 'False', 'False']
table_dict_2['TipoDato']            = ['varchar(20)',   'False', 'False', 'False']
table_dict_2['Descrizione']         = ['varchar(100)',  'False', 'False', 'False']
table_dict_2['Seniority']           = ['INT(4)',        'False', 'False', 'False']
table_dict_2['Emittente']           = ['INT(4)',        'False', 'False', 'False']
table_dict_2['Contributor']         = ['varchar(20)',   'False', 'False', 'False']
table_dict_2['FareScarico']         = ['INT(4)',        'False', 'False', 'False']


# ------------ DProCDS --------------------------------------------------------------
table_dict_3 = {}
table_dict_3['Rating']              = ['INT(4)',        'False', 'False', 'False']
table_dict_3['Currency']            = ['varchar(20)',   'False', 'False', 'False']
table_dict_3['MaturityInt']         = ['INT(4)',        'False', 'False', 'False']
table_dict_3['TipoTicker']          = ['varchar(10)',   'False', 'False', 'False']
table_dict_3['SettorePaese']        = ['INT(4)',        'False', 'False', 'False']
table_dict_3['BloombergTicker']     = ['varchar(50)',    'True', 'False', 'False']
table_dict_3['Descrizione']         = ['varchar(100)',  'False', 'False', 'False']
table_dict_3['Seniority']           = ['INT(4)',        'False', 'False', 'False']
table_dict_3['Emittente']           = ['INT(4)',        'False', 'False', 'False']
table_dict_3['Contributor']         = ['varchar(20)',   'False', 'False', 'False']
table_dict_3['FareScarico']         = ['INT(4)',        'False', 'False', 'False']


# ------------view_bond_master --------------------------------------------------------
table_dict_4 = {}
table_dict_4['id']                 = [ 'VARCHAR(30)'  ,  'False', 'False', 'False' ]
table_dict_4['isin']               = [ 'VARCHAR(30)'  ,  'True', 'False', 'False' ]
table_dict_4['Descrizione']        = [ 'VARCHAR(100)' ,  'False', 'False', 'False' ]
table_dict_4['Divisa']             = [ 'VARCHAR(30)'  ,  'False', 'False', 'False' ]
table_dict_4['Scadenza']           = [ 'DATE'  ,         'False', 'False', 'False' ]
table_dict_4['MaturityType']       = [ 'VARCHAR(30)'  ,  'False', 'False', 'False' ]
table_dict_4['ValoreRimborso']     = [ 'double'  ,       'False', 'False', 'False' ]
table_dict_4['DayCountConvention'] = [ 'VARCHAR(30)'  ,  'False', 'False', 'False' ]
table_dict_4['CouponType']         = [ 'VARCHAR(30)'  ,  'False', 'False', 'False' ]
table_dict_4['FrequenzaCoupon']    = [ 'INT'  ,          'False', 'False', 'False' ]
table_dict_4['Indicizzazione']     = [ 'VARCHAR(30)'  ,  'False', 'False', 'False' ]
table_dict_4['Spread']             = [ 'double'  ,       'False', 'False', 'False' ]
table_dict_4['PUAdjustmentRule']   = [ 'VARCHAR(30)'  ,  'False', 'False', 'False' ]
table_dict_4['tipoStrumento']      = [ 'VARCHAR(30)'  ,  'False', 'False', 'False' ]
table_dict_4['PUCollatType']       = [ 'VARCHAR(30)'  ,  'False', 'False', 'False' ]
table_dict_4['FixingMethod']       = [ 'VARCHAR(30)'  ,  'False', 'False', 'False' ]
table_dict_4['Contributor']        = [ 'VARCHAR(30)'  ,  'False', 'False', 'False' ]
table_dict_4['Emittente']          = [ 'INT'  ,          'False', 'False', 'False' ]
table_dict_4['data']               = [ 'DATE'  ,          'True', 'False', 'False' ]
table_dict_4['CedolaInCorso']      = [ 'double'  ,       'False', 'False', 'False' ]
table_dict_4['Rating']             = [ 'INT'  ,          'False', 'False', 'False' ]
table_dict_4['Settore']            = [ 'INT'  ,          'False', 'False', 'False' ]
table_dict_4['CleanDirty']         = [ 'VARCHAR(30)'  ,  'False', 'False', 'False' ]
table_dict_4['prezzoBid']          = [ 'double'  ,       'False', 'False', 'False' ]
table_dict_4['prezzoMid']          = [ 'double'  ,       'False', 'False', 'False' ]
table_dict_4['prezzoAsk']          = [ 'double'  ,       'False', 'False', 'False' ]
table_dict_4['volume']             = [ 'double'  ,       'False', 'False', 'False' ]
table_dict_4['yieldToMaturityBid'] = [ 'double'  ,       'False', 'False', 'False' ]
table_dict_4['yieldToMaturityMid'] = [ 'double'  ,       'False', 'False', 'False' ]
table_dict_4['yieldToMaturityAsk'] = [ 'double'  ,       'False', 'False', 'False' ]
table_dict_4['discountMarginBid']  = [ 'double'  ,       'False', 'False', 'False' ]
table_dict_4['discountMarginMid']  = [ 'double'  ,       'False', 'False', 'False' ]
table_dict_4['discountMarginAsk']  = [ 'double'  ,       'False', 'False', 'False' ]
table_dict_4['dataScarico']        = [ 'DATE'  ,         'False', 'False', 'False' ]
table_dict_4['LastUpdate']         = [ 'DATE'  ,         'False', 'False', 'False' ]
table_dict_4['calctype']           = [ 'VARCHAR(30)'  ,  'False', 'False', 'False' ]
table_dict_4['calctypedes']        = [ 'VARCHAR(30)'  ,  'False', 'False', 'False' ]
table_dict_4['datagodimento']      = [ 'DATE'  ,         'False', 'False', 'False' ]
table_dict_4['issuedate']          = [ 'DATE'  ,         'False', 'False', 'False' ]
table_dict_4['issueprice']         = [ 'double'  ,       'False', 'False', 'False' ]
table_dict_4['isinflation']        = [ 'VARCHAR(30)'  ,  'False', 'False', 'False' ]
table_dict_4['principalfactor']    = [ 'VARCHAR(30)'  ,  'False', 'False', 'False' ]


# ------------ DProOption --------------------------------------------------------------
table_dict_5 = {}
table_dict_5['TipoDato']            = ['VARCHAR(20)',   'False', 'False', 'False']
table_dict_5['MaturityInt']         = ['INT(4)',        'False', 'False', 'False']
table_dict_5['Mercato']             = ['INT(4)',        'False', 'False', 'False']
table_dict_5['FareScarico']         = ['INT(4)',        'False', 'False', 'False']
table_dict_5['Descrizione']         = ['VARCHAR(100)',  'False', 'False', 'False']
table_dict_5['Currency']            = ['VARCHAR(20)',   'False', 'False', 'False']
table_dict_5['Contributor']         = ['VARCHAR(20)',   'False', 'False', 'False']
table_dict_5['BloombergTicker']     = ['VARCHAR(50)',   'True', 'False',  'False']
table_dict_5['Strike']              = ['double',        'False', 'False', 'False']
table_dict_5['maturityDate']        = ['DATE',          'False', 'False', 'False']
table_dict_5['TipoTicker']          = ['VARCHAR(10)',   'False', 'False', 'False']
table_dict_5['TipoOpzione']         = ['VARCHAR(50)',   'False', 'False', 'False']
table_dict_5['TickerSottostante']   = ['VARCHAR(50)',   'False', 'False', 'False']

# ------------ DProCFS --------------------------------------------------------------
table_dict_6 = {}
table_dict_6['TipoDato']            = ['VARCHAR(20)',   'False', 'False', 'False']
table_dict_6['MaturityInt']         = ['INT(4)',        'False', 'False', 'False']
table_dict_6['FareScarico']         = ['INT(4)',        'False', 'False', 'False']
table_dict_6['Descrizione']         = ['VARCHAR(100)',  'False', 'False', 'False']
table_dict_6['Currency']            = ['VARCHAR(20)',   'False', 'False', 'False']
table_dict_6['Contributor']         = ['VARCHAR(20)',   'False', 'False', 'False']
table_dict_6['BloombergTicker']     = ['VARCHAR(50)',   'True', 'False',  'False']
table_dict_6['Tenor']               = ['INT(4)',        'False', 'False', 'False']
table_dict_6['Strike']              = ['double',        'False', 'False', 'False']
table_dict_6['TipoTicker']          = ['VARCHAR(10)',   'False', 'False', 'False']
table_dict_6['TipoOpzione']         = ['VARCHAR(50)',   'False', 'False', 'False']
table_dict_6['deltaBp']             = ['INT(11)',       'False', 'False', 'False']

#----------------------------------------------------------------------------------


#----------------------------------------------------------------------------------


table_dict = {}
table_dict['DProTS_master']    = table_dict_1
table_dict['DProCurve']        = table_dict_2
table_dict['DProCDS']          = table_dict_3
table_dict['bond_master']      = table_dict_4
table_dict['DProOptions']      = table_dict_5
table_dict['DProCFS']          = table_dict_6

