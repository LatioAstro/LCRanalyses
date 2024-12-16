import numpy as np

#Based on Equation (34) from K. K. Singh & P. J. Meintjes, Astronomische Nachrichten, Volume 341, Issue 713, pp. 713-725
#Error calculated with Equation (37) of the same paper
#Fvar quantifies the variability of a dataset in terms of its variance with respect to the mean
#It accounts for the errors in the measurements by subtracting the sum of the squared errors to the total variance of the data
#This leads sometimes to squared roots of negative values (when the sum of the squared errors is larger than the variance), 
#meaning that we cannot statistically claim that the source is significantly variable
#Alternatively, we can also quantify the variability through a similar parameter called the varibility amplitude (see the 
#corresponding py file with its calculation), that neglects the errors, assigning a value of the variability to all datasets even 
#when errors are large

def fvar(flux,flux_err):
	#Based on Equation (34) from K. K. Singh & P. J. Meintjes, Astronomische Nachrichten, Volume 341, Issue 713, pp. 713-725
	#Error calculated with Equation (37) of the same paper

	N=len(flux)
	if N < 1:
		fvar, fvar_err = np.NaN, np.NaN
	else:
		variance=np.var(flux) #variance
		mean_squared_error=(1/N)*np.sum(flux_err**2)
		mean_flux=np.mean(flux)
		# print("N, variance, mean_squared_error, mean_flux", N, variance, mean_squared_error, mean_flux)

		num = variance-mean_squared_error
		den = mean_flux**2

		if num<0:
			fvar, fvar_err = np.NaN, np.NaN
		else:
			fvar=np.sqrt(num/den)

			fvar_err=np.sqrt(fvar**2 + np.sqrt((2/N) * (mean_squared_error/den)**2 + (mean_squared_error/N) * (2*fvar/mean_flux)**2 )) - fvar

		# print('num, den, fvar, fvar_err', num, den, fvar, fvar_err)


	return fvar,fvar_err


def fvar_percent(flux,flux_err):
	#Based on Equation (34) from K. K. Singh & P. J. Meintjes, Astronomische Nachrichten, Volume 341, Issue 713, pp. 713-725
	#Error calculated with Equation (37) of the same paper
	#Fvar value (and error) multiplied by 100 to express the fractional variability as a "percentage-scaled" value

	fvar, fvar_err = fvar(flux,flux_err)
	
	return fvar*100,fvar_err*100

