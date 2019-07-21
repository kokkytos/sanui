import csv
import os
import pandas as pd
import numpy as np

files = [
        'Amman/corr.table.Amman.csv',
        'Ankara/corr.table.Ankara.csv',
        'Athens/corr.table.Athens.csv',
        'Brussels/corr.table.Brussels.csv',        
        'Cairo/corr.table.Cairo.csv',
        'London/corr.table.London.csv',
        'Milano/corr.table.Milano.csv',
        'Moscow/corr.table.Moscow.csv',
        'Paris/corr.table.Paris.csv',
        'Frankfurt/corr.table.Frankfurt.csv',
        ]

# correlation table with VIIRS
corviirs = pd.DataFrame(columns=[ os.path.basename(elem).split(".")[2] for elem in files])
corviirs['']=[
        'NDVI',
        'SeaWinds',
        'stable lights',
        'VANUI',
        "TVANUI",
        "New method"          
]


corviirs.set_index('', inplace =True)


corolscal = pd.DataFrame.copy(corviirs)

for mycsv in files:
	with open(mycsv, 'rb') as csvfile:
                df = pd.read_csv(mycsv)
                df.set_index('Unnamed: 0', inplace = True)
                
                mypath = os.path.basename(mycsv).split(".")[2]
                
                ena_ndvi =  df.at["1-NDVI[0-1]",'VIIRS']
                seawinds = df.at["SeaWinds",'VIIRS']
                stable_lights = df.at["OLS (stable lights)",'VIIRS']
                vanui = df.at["VANUI",'VIIRS']
                tvanui = df.at["TVANUI",'VIIRS']
                new_method = df.at["SANUI",'VIIRS']       
                
                corviirs[mypath] = ['%.2f' % elem for elem in [ena_ndvi, seawinds,stable_lights, vanui, tvanui, new_method ]]


                
                ena_ndvi =  df.at["1-NDVI[0-1]",'OLS (radiance calibrated)']
                seawinds = df.at["SeaWinds",'OLS (radiance calibrated)']
                stable_lights = df.at["OLS (stable lights)",'OLS (radiance calibrated)']
                vanui = df.at["VANUI",'OLS (radiance calibrated)']
                new_method = df.at["SANUI",'OLS (radiance calibrated)']       
                tvanui = df.at["TVANUI",'OLS (radiance calibrated)']

                corolscal[mypath] = ['%.2f' % elem for elem in [ena_ndvi, seawinds,stable_lights, vanui, tvanui, new_method ]]



corviirs.to_csv('viirs.csv')
corolscal.to_csv('ols_cal.csv')
                
print corviirs
