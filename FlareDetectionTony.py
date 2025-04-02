import numpy as np
import HopFinder
import HOP
import LC
import LC_Set
import matplotlib.pyplot as plt

def quiescent_background_finder(sourcelightcurve,method):
    # Determines the "quiescent background" of a given lightcurve
    # Takes the average of all flux values not within HOP objects.
    # !!!!! MUST USE get_bblocks AND find_hop METHODS BEFORE RUNNING THIS FUNCTION !!!!!!
    qui = sourcelightcurve
    mask = []
    j=0


    if method == 'forward':
    
        #forward. Compares time current time bin and next time bin when finding tmid.
        for j in range(0,len(sourcelightcurve.hops)):
            for i in range(0,len(sourcelightcurve.time)-1):
                tmid = (sourcelightcurve.time[i]+sourcelightcurve.time[i+1])/2
                if tmid == sourcelightcurve.hops[j].start_time:
                    start = i
                if tmid == sourcelightcurve.hops[j].end_time:
                    end = i
            mask.append((start,end))
        maskindices = np.concatenate([np.arange(start,end) for start,end in mask])
        qui.flux = np.delete(qui.flux,maskindices)
        quiescent_background = np.mean(qui.flux)

    if method == 'back':
    #back. Compares time current time bin and previous bin when finding tmid.
        for j in range(0,len(sourcelightcurve.hops)):
            for i in range(0,len(sourcelightcurve.time)-1):
                tmid = (sourcelightcurve.time[i-1]+sourcelightcurve.time[i])/2
                if tmid == sourcelightcurve.hops[j].start_time:
                    start = i
                if tmid == sourcelightcurve.hops[j].end_time:
                    end = i
            mask.append((start,end))
        maskindices = np.concatenate([np.arange(start,end) for start,end in mask])
        qui.flux = np.delete(qui.flux,maskindices)
        quiescent_background = np.mean(qui.flux)

    return quiescent_background


def LCTimeRange(sourcearray, timerangestart,timerangeend):
    sourcearray = sourcearray[timerangestart:timerangeend]
    return sourcearray
    

def quiescent_flare_plot(cadence_df,sourcename,sourcenum,percent):
    
    MJDREFI = 51910
    MJDREFF = 7.428703703703703e-4
    SecsInDay = 86400

    if sourcename == None:
        sourcearray = cadence_df[cadence_df['source_name'] == cadence_df['source_name'][sourcenum]]
        titlestring=cadence_df['source_name'][sourcenum]
    else:
        sourcearray = cadence_df[cadence_df['source_name'] == sourcename].reset_index(drop=True)
        titlestring=sourcename

    sourcearray = sourcearray[sourcearray['photon_flux2']!=-3333].reset_index(drop=True)
    


    time = sourcearray['tmin']+(MJDREFI/SecsInDay)
    photon_flux = sourcearray['photon_flux2']
    errors = sourcearray['photon_flux_error2']
    sourcelightcurve=LC.LightCurve(time,photon_flux,errors,time_format='mjd')
    titlestring=cadence_df['source_name'][sourcenum]

    maxflux=np.max(photon_flux)
    minflux =np.min(photon_flux)
    delta_flux = maxflux - minflux
    delta_flux_percent = delta_flux * percent
    thresholdflux = minflux+delta_flux_percent


    # Finding first set of flares using threshold flux.
    sourcelightcurve.get_bblocks(gamma_value=0.05)
    sourcelightcurve.find_hop(method = 'baseline', lc_edges ='neglect', baseline = thresholdflux)

    # Finding quiescent background.
    quiescent_background = quiescent_background_finder(sourcelightcurve,'forward')


    # Using quiescent background to find flares again.
    sourcelightcurve = LC.LightCurve(time,photon_flux,errors,time_format='mjd')
    #sourcelightcurve.flux = np.subtract(sourcelightcurve.flux,quiescent_background)
    sourcelightcurve.get_bblocks(gamma_value=0.05)
    #sourcelightcurve.get_bblocks_above(threshold = 0)
    sourcelightcurve.find_hop(method = 'baseline', lc_edges ='neglect',baseline = quiescent_background)

    # Plotting the Lightcurve itself.
    plt.figure(figsize=(16,9))
    plt.xlabel("Time (s)")
    plt.ylabel('Photon Flux (Photons/$cm^2\u22c5s^{-1}$)')
    plt.title("Photon Flux vs Time" ' (Source: ' +str(titlestring)+ ')' )
    sourcelightcurve.plot_bblocks(size=2)
    sourcelightcurve.plot_hline(value = quiescent_background, color='green',label='detection threshold',lw=3,linestyle = 'dashed')
    sourcelightcurve.plot_hop()
    plt.legend()