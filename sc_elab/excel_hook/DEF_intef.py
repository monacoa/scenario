import sys
import os

FORMATT = "gg/mm/aaaa"

fitt_translate ={'0' : 'linear', '1':'AVD', '2':'SVE', '3':'CIR', '4':'NS' }

nameSheetCurve         = 'SwapCurve'
nameSheetBootstrap     = 'ElabSwapCurve'
nameSheetCDSBootSurv   = 'Surv Prob CDS'
nameSheetCDSBootMarginal = 'Marginal default'
nameSheetCDSBootHazard = 'Hazard rate'
nameSheetCDSBootSpread = 'Spread CDS'
nameSheetCDSBootZcSpread = 'ZC Spread CDS'
nameSheetCDSBootPySpread = 'PY Spread CDS'
nameSheetBVolBoot       = 'Spot volatilities'

nameSheetBondSurv = 'Surv Prob BOND'
nameSheetBondSpread = 'Spread BOND'
nameSheetBondResults = 'Bond results'
nameSheetBondFittig    = 'Bond fitting'
nameSheetResBondFittig    = 'Results Bond fitting'
nameSheetBondOptPrms = 'Opt prms'


nameSheetCDS      = 'CDS'
nameSheetBond     = 'Bond data'
nameSheetCalib    = 'Calib data'
nameSheetCalibRes = 'Calib result'

#INSERITA PER IL CARICAMENTO DATI, ANCORA DA UTILIZZARE
nameSheetLoadAnagrafica = 'Anagrafica'
nameSheetLoadDati = 'Dati'

#INSERITA PER LO SCARICO DELLA MATRICE SWAPTION
nameSheetScaricoSwaption='Matrix_Swaption'

#INSERITA PER L'ELABORAZIONE DELLE MATRICI
nameSheetElabMatrix = 'Matrix_Elab'
nameSheetElabMatrixResult = 'Matrix_Elab_Result'

#INSERITA PER LO SCARICO DELLE VOLATILITA FLAT IMPLICITE NEI CAP E FLOOR
nameSheetCapFloorVolatilities = 'Cap Floor Vols'
