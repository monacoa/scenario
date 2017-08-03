
# campo excel: [campo db, valore di default]

field_map = {'Isin':			['isin', None],
			'Descrizione':		['Descrizione', 'Sconosciuto'],
			'Seniority':		['PUCollatType', 'Senior unsecured'],
			'Tipo tasso':	    ['CouponType', 'Fixed'],
			'Tipo quotazione':	['CleanDirty', 'Clean'],
			'Tipo rimborso':	[None, 'Bullet'],
			'Prezzo emissione':	['issueprice', 100.0],
			'Prezzo rimborso / Inflation Ratio':	['ValoreRimborso', 'principalfactor'],
			'Data emissione':	['issuedate'],
			'Data scadenza':	['Scadenza'],
			'Tempo scadenza': 	[None, None],
			'Giorni di fixing':	[None, 0],
			'Tipo fixing':		['FixingMethod', None],
			'Basis':			['DayCountConvention', 'ACT/ACT'],
			'Adjustment':		['PUAdjustmentRule', 'Modified following'],	
			'Periodicita cedola (mesi)':['FrequenzaCoupon', 6],
			'Tenor del tasso floater (anni)':['FrequenzaCoupon', 6],	
			'Tasso cedolare annuo (Fisso/spread)':	['Spread', None],
			'Cedola in corso':	['CedolaInCorso', None],
			'Prezzo-MID':	['prezzoMid', None],
			'Prezzo-BID':	['prezzoBid', None],
			'Prezzo-ASK':	['prezzoBid', None],
			'YTM/DM (MID)':	['yieldToMaturityMid', None],
			'YTM/DM (BID)':	['discountMarginBid', None],
			'YTM/DM (ASK)':	['yieldToMaturityAsk', None],
			'Tasso di riferiemnto':	[None, None],
			'Tasso repo': 			[None, None],
			'Data prezzo di mercato':	['LastUpdate', None],
			'Contributor':	['Contributor',None],
			'Peso':	[None, 1],
			'Indicizzazione':	['Indicizzazione', None]}



















