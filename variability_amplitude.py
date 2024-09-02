import numpy as np

#Based on Equations (25) and (26) from K. K. Singh & P. J. Meintjes, Astronomische Nachrichten, Volume 341, Issue 713, pp. 713-725
#Error calculated with Equation (27-28) of the same paper
#See also Romero, Cellone & Combi, A&ASS, 135, 477 (1999)
#Variability amplitude Amp is a parameter used for quantifying the variability of a dataset, analogous to the Fvar, but neglecting the errors of the data. Therefore, even when errors are large and Fvar cannot be computed, Amp values can still be derived under this caveat
#It estimates the amount of variability as the difference between the maximum and minimum values with respect to the mean
#The only limitation for Amp happens when the difference between the maximum and minimum is smaller than twice the squared error of the data, meaning that we cannot claim significant variability
#Since we have different errors for each point, instead of twice the squared error we can the sum of the squared errors of the maximum and minimum flux points

def amp(flux,flux_err):
    flux_max = np.max(flux)
    flux_min = np.min(flux)

    max_pos = np.argmax(flux)
    min_pos = np.argmin(flux)

    sigma_max = flux_err[max_pos]
    sigma_min = flux_err[min_pos]

    diff = flux_max - flux_min

    mean_flux = np.mean(flux)

    #Error of the mean, i.e. standard deviation divided by the squared root of the number of data points
    sigma_mean = np.std(flux_err)/np.sqrt(len(flux_err)) 

    if diff**2 < 2 * sigma_mean**2:
        amp, amp_percent, amp_percent_err = np.NaN, np.NaN, np.NaN

    else:    
        amp = np.sqrt(diff**2 - (sigma_max**2 + sigma_min**2))

        amp_percent = amp * 100/mean_flux

        amp_percent_err = 100 * diff/(mean_flux * amp) * np.sqrt((sigma_max/mean_flux)**2 + (sigma_min/mean_flux)**2 + (sigma_mean/diff)**2 * amp**4)

    return amp_percent, amp_percent_err

        

    

    