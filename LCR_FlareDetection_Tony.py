import numpy as np
import pandas as pd


# This function takes the dataframe of sources of variable objects, and the name of a source.
# It will sort through every entry in the dataframe, and detect consecutive entries with a flux level higher than the average.
# It counts the duration of these flaring states, and notes them in an array 'final', which it returns upon completion.
# Each row represents a flaring state, with characteristics of the flaring state in the following order:
# Name of the source,, duration of the flaring state, average flux, average flux error, and the redshift of the source.
df_4lacdr3 = pd.read_csv('4lac_redshifts.csv', delimiter='\t', comment='#')
def flare_detector(name,dataframe):
    final = [0,0,0,0,0]
    temp =[]
    temp_err=[]
    duration = []
    average=[]
    # This section is for setting up our dataframes which will be used in the following loop.
    # We primarily are eliminating outliers and NaN values that were marked with -3333.
    sourcematrix = dataframe.loc[dataframe['source_name'] == name].reset_index(drop=True)
    fluxmatrix = sourcematrix[sourcematrix['photon_flux2']!=-3333]['photon_flux2'].reset_index(drop=True)
    averageflux = np.mean(sourcematrix[sourcematrix['photon_flux2']!=-3333]['photon_flux2'])
    errormatrix = sourcematrix[sourcematrix['photon_flux_error2']!=-3333]['photon_flux_error2'].reset_index(drop=True)
    boolmatrix = fluxmatrix/averageflux > 1  
    boolmatrix = boolmatrix.reset_index(drop=True)
    if len(df_4lacdr3.loc[(df_4lacdr3['redshift']>-3333)&(df_4lacdr3['4FGL_name']==name)]) > 0:
        sourceredshift = df_4lacdr3.loc[(df_4lacdr3['redshift']>-3333)&(df_4lacdr3['4FGL_name']==name)].values[0][1]
    else:
        sourceredshift = -3333
    print(sourceredshift)
    
    for i in range(len(fluxmatrix)):
        # If the flux value for index i is greater than the average, then this condition is met.
        if boolmatrix[i] == True:
            temp.append(fluxmatrix[i])
            temp_err.append(errormatrix[i])
        # If the flux value is not above average, and there are values stored in the temp array, we gather them and send them to 'final'.
        elif len(temp) != 0:
            if i != 0:
                duration = ((len(temp)))
                average = np.nanmean(temp)
                error = ((sum(k*k for k in temp_err))**0.5)*(1/(len(temp_err)))
                temp = []
                temp_err = []
                if len(final) == 0:
                    final = [name,duration,average,error,sourceredshift]
                else:
                    final = np.vstack([final,[name,duration,average,error,sourceredshift]])
        # Otherwise, we add a row of 0's to show that the source has no flaring states.
        #else:
        #    final = np.vstack([final,[0,0,0,0,0]])


    return final