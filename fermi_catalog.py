import numpy as np
import pandas as pd

# The LCR is based on the Fermi 4FGL-DR3 catalog:
# https://arxiv.org/pdf/2201.11184
# Here we collect functions to extract information from teh Fermi LAT catalog
# For the 4FGL query you can pull any column available in the 4FGL catalog 
# (same keyword names as in Table 12 of https://arxiv.org/pdf/1902.10045.pdf)

"""
--------------------------------------------------------------------------
Read the csv file with the 4FGL catalog.
variables available for each source:
Source_Name, 
Flux1000, 
Energy_Flux100, 
PL_Index, 
Variability_Index,      Variability Index
Frac_Variability,       Fractional Variability
Uni_Frac_variability    1-sigma err on the Fractional Variability
CLASS1 					type of the object 
ASSOC1                  most propbable association to known source
--------------------------------------------------------------------------
"""

df_4fgldr3 = pd.read_csv('4fgl-dr3_LCR.csv', sep = ',', dtype={'Source_Name': str,
                                                               'Flux1000': float, 
																'Energy_Flux100': float, 
																'PL_Index': float, 
																'Variability_Index': float,
																'Frac_Variability': float,
																'Uni_Frac_variability': float,
																'CLASS1': str, 
																'ASSOC1': str} )

var_index_th = 21.67


def select_allblazars():
	print('-----------------------------')
	print('Selecting BL Lac type blazars')
	blz_df = df_4fgldr3.loc[(df_4fgldr3['Variability_Index'] >= var_index_th) & \
						    ((df_4fgldr3['CLASS1'] == 'bll') | (df_4fgldr3['CLASS1'] == 'BLL') | \
						    (df_4fgldr3['CLASS1'] == 'fsrq') | (df_4fgldr3['CLASS1'] == 'FSRQ') | \
						    (df_4fgldr3['CLASS1'] == 'bcu'))]
	print('Total number of variable objects:', len(blz_df))
	return blz_df

def select_bll():
	print('-----------------------------')
	print('Selecting BL Lac type blazars')
	bll_df = df_4fgldr3.loc[(df_4fgldr3['Variability_Index'] >= var_index_th) & \
						   ((df_4fgldr3['CLASS1'] == 'bll') | (df_4fgldr3['CLASS1'] == 'BLL'))]
	print('Total number of variable objects:', len(bll_df))
	return bll_df

def select_fsrq():
	print('-----------------------------')
	print('Selecting FSRQ type blazars')
	fsrq_df = df_4fgldr3.loc[(df_4fgldr3['Variability_Index'] >= var_index_th) & \
						    ((df_4fgldr3['CLASS1'] == 'fsrq') | (df_4fgldr3['CLASS1'] == 'FSRQ'))]
	print('Total number of variable objects:', len(fsrq_df))
	return fsrq_df

def select_bcu():
	print('-----------------------------')
	print('Selecting unknown type blazars')
	bcu_df = df_4fgldr3.loc[(df_4fgldr3['Variability_Index'] >= var_index_th) &
						    (df_4fgldr3['CLASS1'] == 'bcu')]
	print('Total number of variable objects:', len(bcu_df))
	return bcu_df


if __name__ == '__main__':
		
	print(select_bll())
	print(select_fsrq())
	print(select_bcu())

	print(select_allblazars())