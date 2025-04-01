
import numpy as np
import HopFinder
import HOP
import LC
import LC_Set
import matplotlib.pyplot as plt

def quiescent_background_finder(sourcelightcurve,method='forward'):
    # Determines the "quiescent background" of a given lightcurve
    # Takes the average of all flux values not within HOP objects.
    # !!!!! MUST USE get_bblocks AND find_hop METHODS BEFORE RUNNING THIS FUNCTION !!!!!!
    qui = sourcelightcurve
    mask = []
    j=0
    start = 0
    end = 0
    if sourcelightcurve.hops == None:
        quiescent_background = np.nanmean(qui.flux)
        qui_err = np.sum((np.array(qui.flux_error)**2)*(1/(len(np.array(qui.flux_error))))**2)**0.5
        print('no flares detected, quiescent background is all bins')
        return quiescent_background,qui_err

    if method == 'forward':
        #forward. Compares time current time bin and next time bin when finding tmid.
        for j in range(0,len(sourcelightcurve.hops)):

            for i in range(0,len(sourcelightcurve.time)-1):
                tmid = (sourcelightcurve.time[i]+sourcelightcurve.time[i+1])/2
                if tmid == sourcelightcurve.hops[j].start_time:
                    start = i
                if tmid == sourcelightcurve.hops[j].end_time:
                    end = i
            if end < start:
                break
            #print(start,end)
            mask.append((start,end))
        maskindices = np.concatenate([np.arange(start,end) for start,end in mask])
        qui.flux = np.delete(qui.flux,maskindices)
        qui.time = np.delete(qui.time,maskindices)
        
        baseaverage = []
        weights = []
        tempavg = []
        
        for i in range(0,len(qui.flux)-1):

            if (qui.time[i+1]-qui.time[i]!= 604800 or i==len(qui.flux)):
                if np.isnan(np.mean(tempavg)):
                    pass
                else:
                    baseaverage.append(np.mean(tempavg))
                    weights.append(len(tempavg))
                    tempavg =[]
            else:
                tempavg.append(qui.flux[i])

        #print(baseaverage,weights)
        qui.flux_error = np.delete(qui.flux_error,maskindices)


        quiescent_background = np.nanmean(np.average(baseaverage,weights = weights))
        qui_err = np.sum((np.array(qui.flux_error)**2)*(1/(len(np.array(qui.flux_error))))**2)**0.5
        print(quiescent_background,qui_err)

    if qui_err == None:
        qui_err = 0

    return quiescent_background, qui_err

def LCTimeRange(sourcearray, timerangestart,timerangeend):
    sourcearray = sourcearray[timerangestart:timerangeend]
    return sourcearray

def quiescent_flare_plot(cadence_df,sourcename=None,sourcenum=0,percent = 0.1, MJDREFI=51910, MJDREFF=7.428703703703703e-4,bkg_err = False,factor = 1):
    # Depending on what we are doing analysis on/how our time is binned, we can change the MJDREFI to other values.
    
    SecsInDay = 86400

    if sourcename == None:
        sourcearray = cadence_df[cadence_df['source_name'] == cadence_df['source_name'][sourcenum]]
        titlestring=cadence_df['source_name'][sourcenum]
    else:
        sourcearray = cadence_df[cadence_df['source_name'] == sourcename].reset_index(drop=True)
        titlestring=sourcename

    sourcearray = sourcearray[sourcearray['photon_flux2']!=-3333].reset_index(drop=True)
    time = sourcearray['tmin']+(MJDREFI/SecsInDay)
    photon_flux = sourcearray['photon_flux2'] * factor
    errors = sourcearray['photon_flux_error2']

    sourcelightcurve = LC.LightCurve(time,photon_flux,errors,time_format='mjd')
    maxflux = np.max(photon_flux)
    minflux =np.min(photon_flux)
    delta_flux = maxflux - minflux
    delta_flux_percent = delta_flux * percent
    thresholdflux = minflux + delta_flux_percent


    # Finding first set of flares using threshold flux.
    sourcelightcurve.get_bblocks(gamma_value=0.05)
    sourcelightcurve.find_hop(method = 'baseline', lc_edges ='add', baseline = thresholdflux)

    #if sourcelightcurve.hops == None:
    #    print("No flares detected for "+str(titlestring))
    #    return None


    # Finding quiescent background.
    quiescent_background, qui_err = quiescent_background_finder(sourcelightcurve,'forward')


    # Using quiescent background to find flares again.
    sourcelightcurve = LC.LightCurve(time,photon_flux,errors,time_format='mjd')


    #sourcelightcurve.flux = np.subtract(sourcelightcurve.flux,quiescent_background)
    sourcelightcurve.get_bblocks(gamma_value=0.05)
    #sourcelightcurve.get_bblocks_above(threshold = 0)
    sourcelightcurve.find_hop(method = 'baseline', lc_edges ='add',baseline = quiescent_background)

    
    # Plotting the Lightcurve itself.
    plt.figure(figsize=(16,9))
    plt.xlabel("Time (s)")
    plt.ylabel('Photon Flux (Photons/$cm^2\u22c5s^{-1}$)')
    plt.title("Photon Flux vs Time" ' (Source: ' +str(titlestring)+ ')' )
    sourcelightcurve.plot_bblocks(size=2)
    sourcelightcurve.plot_hline(value = quiescent_background, color='green',label='detection threshold',lw=3,linestyle = 'dashed')

    if bkg_err == True:
        y1 = quiescent_background + qui_err
        y2 = quiescent_background - qui_err
        plt.fill_between(range(len(sourcelightcurve.time)),y1,y2,alpha = 0.3)

    sourcelightcurve.plot_hop()
    plt.legend()