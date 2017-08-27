import sys
import os

FORMATT = "gg/mm/aaaa"

fitt_translate ={'0' : 'linear', '1':'AVD', '2':'SVE', '3':'CIR', '4':'NS' }

nameSheetCurve         = 'SwapCurve'
nameSheetBootstrap     = 'ElabSwapCurve'
nameSheetCDSBootSurv   = 'Surv Prob CDS'
nameSheetCDSBootSpread = 'Spread CDS'

nameSheetBondSurv = 'Surv Prob BOND'
nameSheetBondSpread = 'Spread BOND'
nameSheetBondResults = 'Bond results'
nameSheetBondFittig    = 'Bond fitting'
nameSheetResBondFittig    = 'Results Bond fitting'

nameSheetCDS     = 'CDS'
nameSheetBond     = 'Bond data'