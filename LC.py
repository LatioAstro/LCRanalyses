
import numpy as np 
from matplotlib import pyplot as plt
import pickle
import astropy
import astropy.stats.bayesian_blocks as bblocks
#https://docs.astropy.org/en/stable/api/astropy.stats.bayesian_blocks.html
from HopFinder import *

import logging
logging.basicConfig(level=logging.ERROR)
"""
set logging to the desired level
logging options:
DEBUG:      whatever happens will be thrown at you
INFO:       confirmation that things are working as expected
WARNING:    sth unexpected happened
ERROR:      sth didn't work, abort mission
""" 

def load_lc_npy(path):
    """
    to load LCs that are saved as numpy array through pickle with lc.save_npy()
    ATTENTION: LC.py does not ge updated! 
    """
    pickle_file = open(path, "rb") #open file and read bytes
    return pickle.load(pickle_file)

def load_lc_csv(path):
    a = np.genfromtxt(path)
    lc = LightCurve(a[0], a[1], a[2])
    return lc

def flux_puffer(flux, threshold, threshold_error):
    """
    ATTENTION! This returns artificial flux values! Use cautiously if at all..
    Set every flux bin under threshold to threshold before initializing light curve
    Apply Bayesian blocks -> detect significant variations wrt threshold = flares?
    """
    flux_new = np.where(flux > threshold, flux, threshold)
    flux_error_new = np.where(flux > threshold, flux_error, th_error)
    return(flux_new, flux_error_new)

def fix_data(time, flux, flux_error):
    """
    ATTENTION! this deletes bins, if there is np.nan in flux(_error) or duplicate in time
    """
    flux_ = flux[np.invert(np.isnan(flux)) * np.invert(np.isnan(flux_error))]
    flux_error_ = flux_error[np.invert(np.isnan(flux)) * np.invert(np.isnan(flux_error))]
    time_ = time[np.invert(np.isnan(flux)) * np.invert(np.isnan(flux_error))]
    logging.info('Deleted ' + str(len(flux) - len(flux_)) + ' np.nan values.')
    unique_time, unique_time_id = np.unique(time_, return_index=True)
    good_flux = flux_[unique_time_id]
    good_flux_error = flux_error_[unique_time_id]
    logging.info('Deleted ' + str(len(time_) - len(unique_time)) + ' time duplicates')
    return(unique_time, good_flux, good_flux_error)


def get_gti_iis(time, n_gaps, n_pick):
    # get index of good time intervals (divide LC into secitons in case there are n_gaps gaps in data; like FACT)
    # biggest time gaps (length in s) in chronological order
    diff = np.array([t - s for s, t in zip(time, time[1:])])
    diff1 = np.sort(diff)
    ii = [x for x in range(len(diff)) if diff[x] in diff1[-n_gaps:]] #index of the 10 longest gaps
    GTI_start_ii = np.array(ii)+1
    GTI_start_ii = np.insert(GTI_start_ii,0,0)
    GTI_end_ii = np.array(ii)
    GTI_end_ii = np.append(GTI_end_ii, len(time)-1)
    if n_pick:
        # only consider the n_pick longest gtis 
        gap_len = np.array([t - s for s,t in zip(GTI_start_ii, GTI_end_ii)])
        gap_len1 = np.sort(gap_len)
        ii = [x for x in range(len(gap_len)) if gap_len[x] in gap_len1[-n_pick:]] # n_gaps = considered gaps (longest not gaps)
        GTI_start_ii_ = GTI_start_ii[ii]
        GTI_end_ii_ = GTI_end_ii[ii]
        return GTI_start_ii_, GTI_end_ii_
    else:
        return GTI_start_ii, GTI_end_ii

def make_gti_lcs(lc, n_gaps, n_pick=None):
    """
    Divide one lc with n_gaps gaps into several lcs with good coverage.
    """
    gti_starts, gti_ends = get_gti_iis(lc.time, n_gaps, n_pick)
    if n_pick is None:
        n_pick = n_gaps + 1 #select all 
    chunks = []
    for g in range(n_pick):
        gti_lc = LightCurve(lc.time[gti_starts[g]:gti_ends[g]+1], 
                            lc.flux[gti_starts[g]:gti_ends[g]+1], 
                            lc.flux_error[gti_starts[g]:gti_ends[g]+1],
                            name=lc.name, z=lc.z)
        chunks.append(gti_lc)
    return(np.array(chunks))


#--------------------------------------------------------------------------------------------------
class LightCurve:
    """
    Light Curve Class
    ------------------
    Create a light curve based on input data: time, flux, flux_error
    """
    def __init__(self, time, flux, flux_error, time_format=None, name=None, z=None, 
                 telescope=None, cadence=None):
        self.time = np.array(time)
        self.flux = np.array(flux)
        self.flux_error = np.array(flux_error)
        self.time_format = time_format
        self.name = name
        self.z = z
        self.telescope = telescope
        self.cadence = cadence
        if len(time) != len(flux) or len(time) != len(flux_error):
            raise ValueError('Input arrays do not have same length')
        if len(flux[np.isnan(flux)]) > 0 or len(flux_error[np.isnan(flux_error)]) > 0:
            raise TypeError('flux or flux_error contain np.nan values')
        if len(time) != len(np.unique(time)):
            raise ValueError('time contains duplicate values')
        if time_format:
            """ format of the astropy.time.Time object """
            self.astropy_time = astropy.time.Time(time, format=time_format)

    def __repr__(self):
        #this is everything you need to know about this instance (eg used for == method)
        # __str__() draws from this (eg used for print)
        return f'LightCurve (bins = {len(self.flux)}, name = {self.name}, '+ \
               f'cadende = {self.cadence}, telescope = {self.telescope}, z = {self.z})'
               #could be extended to make 100% sure this is unambiguous (eg with bblocks)
    
    def __len__(self):
        return len(self.time)

    def __getitem__(self, inbr):
        """
        overwriting getitem method to access
        * one LC bin like: lc[i]
        * a slice of an LC (= return new LC) like: lc[i:ii]
        * selected bins of an LC like: lc[[i1, i2, i3, i4]]
        """
        if type(inbr) is int:
            return np.array([self.time[inbr], self.flux[inbr], self.flux_error[inbr]])
        elif type(inbr) is slice: 
            return LightCurve(self.time[inbr], self.flux[inbr], 
                              self.flux_error[inbr], self.time_format, 
                              self.name, self.z, self.telescope, self.cadence)
        elif type(inbr) is list:
            #can't be implemented with 'int or list' -> confusion with slice
            return np.array([self.time[inbr], self.flux[inbr], self.flux_error[inbr]])

    def select_by_time(self, t_min, t_max):
        # select certain part of the light curve based on start and end TIME
        i_s = np.where(self.time >= t_min)[0][0]
        i_e = np.where(self.time >= t_max)[0][0]
        return self.__getitem__(slice(i_s, i_e, None))

    def save_npy(self, path):
        """
        save LC as numpy arraw with pickle
        to open this see load_lc_npy() at top 
        ATTENTION: LC.py does not ge updated! 
        """
        pickle_file = open(path,"wb") #open file and write in bytes
        pickle.dump(self, pickle_file)
        pickle_file.close()

    def save_csv(self, path, bblocks=False):
        #to open this see load_lc_csv() at top 
        if bblocks is True:
            data = np.array([self.time, self.flux, self.flux_error, self.block_pbin])
            np.savetxt(path, data, comments='#time, flux, flux_error, block_pbin')
        else:
            data = np.array([self.time, self.flux, self.flux_error])
            np.savetxt(path, data, comments='#time, flux, flux_error')

    def plot_lc(self, data_color='k', ax=None, new_time_format=None, size=1, **kwargs):
        if ax is None:
            ax = plt.gca()
        ax.errorbar(x=self.time, y=self.flux, yerr=self.flux_error, ecolor=data_color, 
                     elinewidth=1, linewidth=0, marker='+', markersize=3*size, 
                     color=data_color, **kwargs)
        if self.time_format and new_time_format is not None:
            axtop = ax.twiny()
            axtop.set_xticks(ax.get_xticks())
            axtop.set_xbound(ax.get_xbound())
            axtop.set_xlim(ax.get_xlim())
            format_labels = astropy.time.Time([t for t in ax.get_xticks()], format=self.time_format)
            if new_time_format == 'isot': 
                new_labels = [format_labels.to_value(format='isot')[i].split('T')[0] for i in range(len(format_labels))]
                axtop.set_xticklabels(new_labels) #= yyyy-mm-dd
            elif new_time_format == 'decimalyear':
                new_labels = format_labels.to_value(format='decimalyear')
                axtop.set_xticklabels(np.round(new_labels, 1)) #= yyyy.y
            else:
                new_labels = format_labels.to_value(format=new_time_format)
                axtop.set_xticklabels(new_labels)
            plt.sca(ax) #go back to initial bottom axis

    def plot_hline(self, value, ax=None, **kwargs):
        if ax is None:
            ax = plt.gca()
        ax.hlines(value, xmin=min(self.time), xmax=max(self.time), **kwargs)

    def plot_vline(self, value, ax=None, **kwargs):
        if ax is None:
            ax = plt.gca()
        ax.vlines(value, ymin=min(self.flux), ymax=max(self.flux), **kwargs)

    def plot_shade(self, start_time, end_time, ax=None, **kwargs):
        if ax is None:
            ax = plt.gca()
        x = np.linspace(start_time, end_time)
        y = np.ones(len(x)) * np.max(self.flux)
        y1 = np.ones(len(x)) * np.min(self.flux)
        ax.fill_between(x, y, y1, step="mid", alpha=0.2, zorder=0, **kwargs)

    def plot_grid(self, spacing=10, ax=None, **kwargs):
        if ax is None:
            ax = plt.gca()
        rounded_start = np.round(np.min(self.time)/spacing) * spacing
        rounded_end = np.round(np.max(self.time)/spacing) * spacing
        ax.set_xticks(np.arange(rounded_start, rounded_end, spacing), minor=True)
        ax.grid(which='minor', **kwargs)

    #----------------------------------------------------------------------------------------------
    def get_bblocks(self, gamma_value=None, p0_value=0.05): 
        """
        Bayesian block algorithm (https://ui.adsabs.harvard.edu/abs/2013arXiv1304.2818S/abstract)
        fitness is set to 'measures' since we assume Gaussian error for flux measurements
        from astropy (https://docs.astropy.org/en/stable/api/astropy.stats.bayesian_blocks.html)
        Returns edges of blocks (significant changes of flux) in units of time (e.g. MJD)
        Edges are converted to edge_index (based on position in time array)
        -> See GitHub description and Jupyter Notebook for more information
        block_val are the flux values of all blocks based on the mean of all flux bins within
        block_val_error is the corresponding error computed with Gaussian error propagation
        block_pbin has the same shape as flux and is filled with corrsponding block values
        """
        # get Bayesian block edges for light curve
        self.edges = bblocks(t=self.time, x=self.flux, sigma=self.flux_error, fitness='measures',
                             gamma=gamma_value, p0=p0_value)
        logging.debug('got edges for light curve')

        if len(self.edges) <= 2:
            logging.warning('light curve is constant; only one bayesian block found.')
            self.block_pbin = np.ones(len(self.flux)) * np.mean(self.flux)
            self.block_val = np.array([np.mean(self.flux)])
            self.block_val_error = np.array([np.std(self.flux)])
            self.edge_index = np.array([0, -1])
            self.edges = np.array([self.time[0], self.time[-1]])
            return(self.block_pbin, self.block_val, self.block_val_error, self.edge_index,
                   self.edges)
        # get edge_index
        self.edge_index = np.array([np.where(self.time >= self.edges[i])[0][0] 
                                    for i,_ in enumerate(self.edges)])
        #change last entry such that loop over [j:j+1] gives all blocks
        self.edge_index[-1] += 1

        # determine flux value (mean) and error (Gaussian propagation) for each block
        self.block_val = np.zeros(len(self.edge_index)-1)  
        self.block_val_error = np.zeros(len(self.edge_index)-1) 
        for j in range(len(self.edge_index)-1):
            self.block_val[j] = np.mean(self.flux[self.edge_index[j]: self.edge_index[j+1]])
            self.block_val_error[j] = (np.sqrt(np.sum(self.flux_error[self.edge_index[j]:
                                                     self.edge_index[j+1]]**2))
                                      / (self.edge_index[j+1]-self.edge_index[j]))

        # create block-per-bin array corresponding to flux
        self.block_pbin = np.zeros(len(self.flux))
        for k,_ in enumerate(self.block_val):
            self.block_pbin[self.edge_index[k] : self.edge_index[k+1]] = self.block_val[k]
        logging.debug('got block parameters for light curve')

        return(self.block_pbin, self.block_val, self.block_val_error, self.edge_index, self.edges)

    #----------------------------------------------------------------------------------------------
    def get_bblocks_above(self, threshold, pass_gamma_value=None, pass_p0_value=None):
        """
        ATTENTION! This returns artificial flux values! Use cautiously if at all..
        Note: get_bblocks has to be applied first
        Determine Bayesian blocks for light curve but set all blocks that are lower than threshold
        to that threshold (i.e. set small block_val to threshold and neglect edges under threshold)
        -> leaves only significant variations wrt threshold = flares?
        """ 
        self.block_pbin = np.where(self.block_pbin > threshold, self.block_pbin, threshold)
        self.block_val = np.where(self.block_val > threshold, self.block_val, threshold)
        block_mask = np.ones(len(self.block_val), dtype = bool)
        edge_mask = np.ones(len(self.edges), dtype=bool)
        for i in range(len(self.block_val)-1):
            if self.block_val[i] == threshold and self.block_val[i+1] == threshold:
                block_mask[i+1] = False
                edge_mask[i+1] = False
        self.block_val = self.block_val[block_mask]
        self.block_val_error = self.block_val_error[block_mask]
        self.edge_index = self.edge_index[edge_mask]
        self.edges = self.edges[edge_mask]
        return(self.block_pbin, self.block_val, self.block_val_error, self.edge_index, self.edges)

    #----------------------------------------------------------------------------------------------
    def plot_bblocks(self, bb_color='steelblue', data_color='k', data_label='obs flux',
                     size=1, ax=None, new_time_format=None):
        if ax is None:
                ax = plt.gca()
        try:   
            ax.step(self.time, self.block_pbin, where='mid', linewidth=1*size, label='bblocks', 
            	     color=bb_color, zorder=1000)
            #ax.errorbar(x=self.time, y=self.flux, yerr=self.flux_error, label=data_label, 
            # 	         ecolor=data_color, elinewidth=1*size, linewidth=0, marker='+', 
            #             markersize=3*size, color=data_color)
            self.plot_lc(data_color=data_color, ax=ax, label=data_label, 
                         new_time_format=new_time_format, size=size)
        except AttributeError:
            raise AttributeError('Initialize Bayesian blocks with .get_bblocks() first!')

    #----------------------------------------------------------------------------------------------
    def bb_i(self, t):
        """
        Convert time to index of corresponding Bayesian block (e.g. block_value of peak_time)
        use bb_i_start/bb_i_end to make sure you get the block left/right outside of hop
        this works fine for flip, halfclap, and sharp but *NOT for BASELINE* (-> block inside hop)
        """
        if t == self.edges[0]:
            return(int(0))
        else:
            block_index = [
                e for e in range(len(self.edges)-1) if t > self.edges[e] and t <= self.edges[e+1]]
            return(int(block_index[0]))

    def bb_i_start(self,t):
        """
        if time = edge -> take block on the left
        ATTENTION: for baseline method this is first block of hop -> use bb_i() instead (works)
        """
        block_index = [
            e for e in range(len(self.edges)-1) if t >= self.edges[e] and t < self.edges[e+1]]
        return(int(block_index[0]))

    def bb_i_end(self,t):
        """
        if time = edge -> take block on the right
        ATTENTION: for baseline method this is last block of hop - use bb_i() instead (TBD)
        """
        block_index = [
            e for e in range(len(self.edges)-1) if t > self.edges[e] and t <= self.edges[e+1]]
        return(int(block_index[0]))
    
    #----------------------------------------------------------------------------------------------
    def find_hop(self, method='half', lc_edges='neglect', baseline=None):
        if method == 'baseline':
            if baseline is None:
                self.baseline = np.mean(self.flux)
            else:
                self.baseline = baseline
            hopfinder = HopFinderBaseline(lc_edges)
        if method == 'half':
            hopfinder = HopFinderHalf(lc_edges)
        if method == 'sharp':
            hopfinder = HopFinderSharp(lc_edges)
        if method == 'flip':
            hopfinder = HopFinderFlip(lc_edges)
        self.hops = hopfinder.find(self)
        return self.hops

    #----------------------------------------------------------------------------------------------
    def plot_hop(self, ax=None, **kwargs):
        """
        Plot shaded area for all hops in light curve
        """
        if self.hops is None:
            return # no hop in this lc
        if ax is None:
            ax = plt.gca()
        for i,hop in enumerate(self.hops):
            x = np.linspace(hop.start_time, hop.end_time)
            y = np.ones(len(x)) * np.max(self.flux)
            y1 = np.min(self.flux)
            if i == 0:
                ax.fill_between(x, y, y1, step="mid", color='lightsalmon', alpha=0.2,
                                 label='hop', zorder=0)
            if i == 1:
                ax.fill_between(x, y, y1, step="mid", color='orchid', alpha=0.2, label='hop',
                                 zorder=0)
            elif i % 2:
                ax.fill_between(x, y, y1, step="mid", color='orchid', alpha=0.2, zorder=0)
            elif i != 0:
                ax.fill_between(x, y, y1, step="mid", color='lightsalmon', alpha=0.2, zorder=0)
        #ax.set_title(lc.name, hop.method)

    #----------------------------------------------------------------------------------------------
    def plot_all_hop(self, gamma_value=None, p0_value=0.05, lc_edges='neglect'):
        """
        Plot all HOP methods in one figure for comparison
        """
        fig = plt.figure(0,(15,9))
        plt.suptitle('All HOP methods', fontsize=16)

        ax0 = fig.add_subplot(511)
        self.find_hop('baseline')
        self.plot_bblocks()
        self.plot_hop()
        plt.ylabel('baseline')

        ax1 = fig.add_subplot(512)
        self.find_hop('half')
        self.plot_bblocks()
        self.plot_hop()
        plt.ylabel('half')

        ax2 = fig.add_subplot(513)
        self.find_hop('flip')
        self.plot_bblocks()
        self.plot_hop()
        plt.ylabel('flip')

        ax3 = fig.add_subplot(514)
        self.find_hop('sharp')
        self.plot_bblocks()
        self.plot_hop()
        plt.ylabel('sharp')
        fig.subplots_adjust(hspace=0)


