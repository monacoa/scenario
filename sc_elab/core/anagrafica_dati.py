import numpy as np

MaturityFromIntToString = dict(zip(360*np.arange(1,70),
                                   [s + 'Y' for s in map(str, range(1, 70))]))
MaturityFromIntToString[1] ='O/N'
MaturityFromIntToString[2] ='T/N'
MaturityFromIntToString[7] ='1W'
MaturityFromIntToString[14] ='2W'
MaturityFromIntToString[15] ='2W'
MaturityFromIntToString[21] ='3W'
MaturityFromIntToString[30] ='1M'
MaturityFromIntToString[60] ='2M'
MaturityFromIntToString[90] ='3M'
MaturityFromIntToString[120] ='4M'
MaturityFromIntToString[150] ='5M'
MaturityFromIntToString[180] ='6M'
MaturityFromIntToString[210] ='7M'
MaturityFromIntToString[240] ='8M'
MaturityFromIntToString[270] ='9M'
MaturityFromIntToString[300] ='10M'
MaturityFromIntToString[330] ='11M'
MaturityFromIntToString[450] ='15M'
MaturityFromIntToString[540] ='18M'
MaturityFromIntToString[630] ='21M'


MaturityFromStringToYear = dict(zip([s + 'Y' for s in map(str, range(1, 70))],np.arange(1.,70.)))
MaturityFromStringToYear['1M']  = 1./12
MaturityFromStringToYear['2M']  = 1./6
MaturityFromStringToYear['3M']  = 0.25
MaturityFromStringToYear['4M']  = 1./3
MaturityFromStringToYear['5M']  = 5./12
MaturityFromStringToYear['6M']  = 0.5
MaturityFromStringToYear['7M']  = 7./12
MaturityFromStringToYear['8M']  = 2./3
MaturityFromStringToYear['9M']  = 0.75
MaturityFromStringToYear['10M'] = 5./6
MaturityFromStringToYear['11M'] = 11./12
MaturityFromStringToYear['18M'] = 1.5