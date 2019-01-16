

from_bond_table_to_datastream = {}
from_bond_table_to_datastream['id']                  = 'ISIN'
from_bond_table_to_datastream['isin']                = 'ISIN'
from_bond_table_to_datastream['Descrizione']         ='NAME'
from_bond_table_to_datastream['PUCollatType']        ='DEBT SENIORITY'
from_bond_table_to_datastream['CouponType']          ='COUPON TYPE'
from_bond_table_to_datastream['MaturityType']        ='AMORTISATION'
from_bond_table_to_datastream['issueprice']          ='ISSUE PRICE'
from_bond_table_to_datastream['ValoreRimborso']      ='REDEMPTION VALUE'
from_bond_table_to_datastream['issuedate']           ='ISSUE DATE'
from_bond_table_to_datastream['Scadenza']            ='REDEMPTION DATES'
from_bond_table_to_datastream['DayCountConvention']  ='ACCRUAL BASIS'
from_bond_table_to_datastream['FrequenzaCoupon']     ='COUPONS PER YEAR'
from_bond_table_to_datastream['Spread']              ='FRN COUP'
from_bond_table_to_datastream['CedolaInCorso']       ='COUPON'
from_bond_table_to_datastream['prezzoMid']           ='COMPOSITE PR MID'
from_bond_table_to_datastream['prezzoBid']           ='COMPOSITE PR BID'
from_bond_table_to_datastream['prezzoAsk']           ='COMPOSITE PR ASK'
from_bond_table_to_datastream['yieldToMaturityBid']  ='RED. YIELD'


from_bond_table_to_datastream['Divisa']  ='Currency'
from_bond_table_to_datastream['Indicizzazione']  ='NULL'
from_bond_table_to_datastream['PUAdjustmentRule']  ='NULL'
from_bond_table_to_datastream['Rating']  ='NULL'

from_bond_table_to_datastream['Contributor']  ='NULL'
from_bond_table_to_datastream['Emittente']  ='NULL'
from_bond_table_to_datastream['tipoStrumento']  ='NULL'
from_bond_table_to_datastream['FixingMethod']  ='NULL'


from_bond_table_to_datastream['data']  ='Data'
from_bond_table_to_datastream['Settore']  ='NULL'
from_bond_table_to_datastream['CleanDirty']  ='TIPO QUOTAZIONE'
from_bond_table_to_datastream['volume']  ='NULL'

from_bond_table_to_datastream['yieldToMaturityMid']  ='NULL'
from_bond_table_to_datastream['yieldToMaturityAsk']  ='NULL'


from_bond_table_to_datastream['discountMarginBid']  ='NULL'
from_bond_table_to_datastream['discountMarginMid']  ='NULL'
from_bond_table_to_datastream['discountMarginAsk']  ='NULL'


from_bond_table_to_datastream['dataScarico']  ='Data'
from_bond_table_to_datastream['LastUpdate']  ='NULL'

from_bond_table_to_datastream['calctype']  ='NULL'
from_bond_table_to_datastream['calctypedes']  ='NULL'
from_bond_table_to_datastream['datagodimento']  ='NULL'

from_bond_table_to_datastream['isinflation']  ='NULL'
from_bond_table_to_datastream['principalfactor']  ='NULL'



min_field_set = ['BloombergTicker','Contributor','Data','Datatype','ValoreBid','ValoreAsk','ValoreMid','LastUpdate','DataScarico','TipoDato','id']


cds_default_values = {}

"""
data_dict_to_insert_tmp['Contributor'][0]       = anag_dict_res['Contributor']
data_dict_to_insert_tmp['Datatype'][0]          = 'Livello'
data_dict_to_insert_tmp['ValoreBid'][0]         = 0.0
data_dict_to_insert_tmp['ValoreAsk'][0]         = 0.0
data_dict_to_insert_tmp['TipoDato'][0]          = anag_dict_res['TipoDato']
data_dict_to_insert_tmp['id'][0]                = data_values['Ticker'][i] + anag_dict_res['TipoTicker']

data_dict_to_insert_tmp['ValoreMid'][0]         = data_values['Valore'][i]
data_dict_to_insert_tmp['LastUpdate'][0]        = data_values['Data'][i]
data_dict_to_insert_tmp['DataScarico'][0]       = data_values['Data'][i]
data_dict_to_insert_tmp['Data'][0]              = data_values['Data'][i]
data_dict_to_insert_tmp['BloombergTicker'][0]   = data_values['Ticker'][i]
"""

"""
from_cds_table_to_datastream = {}
from_cds_table_to_datastream['id']                  = 'ISIN'
from_cds_table_to_datastream['isin']                = 'ISIN'
from_cds_table_to_datastream['Descrizione']         ='NAME'
from_cds_table_to_datastream['PUCollatType']        ='DEBT SENIORITY'
from_cds_table_to_datastream['CouponType']          ='COUPON TYPE'
from_cds_table_to_datastream['MaturityType']        ='AMORTISATION'
from_cds_table_to_datastream['issueprice']          ='ISSUE PRICE'
from_cds_table_to_datastream['ValoreRimborso']      ='REDEMPTION VALUE'
from_cds_table_to_datastream['issuedate']           ='ISSUE DATE'
from_cds_table_to_datastream['Scadenza']            ='REDEMPTION DATES'
from_cds_table_to_datastream['DayCountConvention']  ='ACCRUAL BASIS'
from_cds_table_to_datastream['FrequenzaCoupon']     ='COUPONS PER YEAR'
from_cds_table_to_datastream['Spread']              ='FRN COUP'
from_cds_table_to_datastream['CedolaInCorso']       ='COUPON'
from_cds_table_to_datastream['prezzoMid']           ='COMPOSITE PR MID'
from_cds_table_to_datastream['prezzoBid']           ='COMPOSITE PR BID'
from_cds_table_to_datastream['prezzoAsk']           ='COMPOSITE PR ASK'
from_cds_table_to_datastream['yieldToMaturityBid']  ='RED. YIELD'
"""




