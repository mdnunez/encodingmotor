# pdm7_resultsmodel1.py - Plots results of preregistered JAGS models
#
# Copyright (C) 2018 Michael D. Nunez, <mdnunez1@uci.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Record of Revisions
#
# Date            Programmers                         Descriptions of Change
# ====         ================                       ======================
# 08/21/17      Michael Nunez                             Original code
# 11/01/18      Michael Nunez                Results of model 'n200suball_EEG_lapse'

# Imports
import os
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy import stats
# Initialize ipython matplotlib plotting graphics
get_ipython().magic('pylab')

# Definitions
def jellyfish(possamps):  # jellyfish plots
    """Plots posterior distributions of given posterior samples in a jellyfish
    plot. Jellyfish plots are posterior distributions (mirrored over their
    horizontal axes) with 99% and 95% credible intervals (currently plotted
    from the .5% and 99.5% & 2.5% and 97.5% percentiles respectively.
    Also plotted are the median and mean of the posterior distributions"

    Parameters
    ----------
    possamps : ndarray of posterior chains where the last dimension is
    the number of chains, the second to last dimension is the number of samples
    in each chain, all other dimensions describe the shape of the parameter
    """

    # Number of chains
    nchains = possamps.shape[-1]

    # Number of samples per chain
    nsamps = possamps.shape[-2]

    # Number of dimensions
    ndims = possamps.ndim - 2

    # Number of variables to plot
    nvars = np.prod(possamps.shape[0:-2])

    # Index of variables
    varindx = np.arange(nvars).reshape(possamps.shape[0:-2])

    # Reshape data
    alldata = np.reshape(possamps, (nvars, nchains, nsamps))
    alldata = np.reshape(alldata, (nvars, nchains * nsamps))

    # Plot properties
    LineWidths = np.array([2, 5])
    teal = np.array([0, .7, .7])
    blue = np.array([0, 0, 1])
    orange = np.array([1, .3, 0])
    Colors = [teal, blue]

    # Initialize ylabels list
    ylabels = ['']

    for v in range(0, nvars):
        # Create ylabel
        whereis = np.where(varindx == v)
        newlabel = ''
        for l in range(0, ndims):
            newlabel = newlabel + ('_%i' % whereis[l][0])

        ylabels.append(newlabel)

        # Compute posterior density curves
        kde = stats.gaussian_kde(alldata[v, :])
        bounds = stats.scoreatpercentile(alldata[v, :], (.5, 2.5, 97.5, 99.5))
        for b in range(0, 2):
            # Bound by .5th percentile and 99.5th percentile
            x = np.linspace(bounds[b], bounds[-1 - b], 100)
            p = kde(x)

            # Scale distributions down
            maxp = np.max(p)

            # Plot jellyfish
            upper = .25 * p / maxp + v + 1
            lower = -.25 * p / maxp + v + 1
            lines = plt.plot(x, upper, x, lower)
            plt.setp(lines, color=Colors[b], linewidth=LineWidths[b])
            if b == 1:
                # Mark mode
                wheremaxp = np.argmax(p)
                mmode = plt.plot(np.array([1., 1.]) * x[wheremaxp],
                                 np.array([lower[wheremaxp], upper[wheremaxp]]))
                plt.setp(mmode, linewidth=3, color=orange)
                # Mark median
                mmedian = plt.plot(np.median(alldata[v, :]), v + 1, 'ko')
                plt.setp(mmedian, markersize=10, color=[0., 0., 0.])
                # Mark mean
                mmean = plt.plot(np.mean(alldata[v, :]), v + 1, '*')
                plt.setp(mmean, markersize=10, color=teal)

    # Display plot
    plt.setp(plt.gca(), yticklabels=ylabels, yticks=np.arange(0, nvars + 1))


def diagnostic(insamples):
    """
    Returns Rhat (measure of convergence, less is better with an approximate
    1.10 cutoff) and Neff, number of effective samples).

    Reference: Gelman, A., Carlin, J., Stern, H., & Rubin D., (2004).
              Bayesian Data Analysis (Second Edition). Chapman & Hall/CRC:
              Boca Raton, FL.


    Parameters
    ----------
    insamples: dic
        Sampled values of monitored variables as a dictionary where keys
        are variable names and values are numpy arrays with shape:
        (dim_1, dim_n, iterations, chains). dim_1, ..., dim_n describe the
        shape of variable in JAGS model.

    Returns
    -------
    dict:
        Rhat for each variable. Prints Maximum Rhat
    """

    result = {}  # Initialize dictionary
    maxrhats = np.zeros((len(insamples.keys())), dtype=float)
    allkeys ={} # Initialize dictionary
    keyindx = 0
    for key in insamples.keys():
        if key[0] != '_':
            result[key] = {}

            possamps = insamples[key]

            nchains = possamps.shape[-1]
            nsamps = possamps.shape[-2]
            # Mean of each chain
            chainmeans = np.mean(possamps, axis=-2)
            # Global mean of each parameter
            globalmean = np.mean(chainmeans, axis=-1)
            result[key]['mean'] = globalmean
            globalmeanext = np.expand_dims(
                globalmean, axis=-1)  # Expand the last dimension
            globalmeanext = np.repeat(
                globalmeanext, nchains, axis=-1)  # For differencing
            # Between-chain variance
            between = np.sum(np.square(chainmeans - globalmeanext),
                             axis=-1) * nsamps / (nchains - 1.)
            # Mean of the variances of each chain
            within = np.mean(np.var(possamps, axis=-2), axis=-1)
            # Total estimated variance
            totalestvar = (1. - (1. / nsamps)) * \
                within + (1. / nsamps) * between
            # Rhat (Gelman-Rubin statistic)
            temprhat = np.sqrt(totalestvar / within)
            maxrhats[keyindx] = np.nanmax(temprhat) # Ignore NANs
            allkeys[keyindx] = key
            result[key]['rhat'] = temprhat
            keyindx += 1
            # Possible number of effective samples?
            # Geweke statistic?
    print "Maximum Rhat: %3.2f of variable %s" % (np.max(maxrhats),allkeys[np.argmax(maxrhats)])
    return result



# Save location
thisdir = os.getcwd()
os.chdir('..')
os.chdir('..')
os.chdir('Data')
os.chdir('jagsout')
loadloc = os.getcwd()

# Initialize ipython matplotlib plotting graphics
get_ipython().magic('pylab')

# modelfile = loadloc + '/' + 'model_n200suball_EEG_lapseOct_31_18_18_30.mat'
# modelfile = loadloc + '/' + 'model_n200subter_EEG_lapseNov_01_18_12_45.mat'
modelfile = loadloc + '/' + 'model_p300suball_EEG_lapseNov_03_18_00_38.mat'

os.chdir(thisdir)

# Plot values
fontsize = 20

# Plot figures

samples = sio.loadmat(modelfile)
diags = diagnostic(samples)

plt.figure()
jellyfish(samples['tercond'])
plt.title('Condition-level residual non-decision time', fontsize=fontsize)
plt.figure()
jellyfish(samples['deltacond'])
plt.title('Condition-level residual drift rate parameter', fontsize=fontsize)
plt.figure()
jellyfish(samples['alphacond'])
plt.title('Condition-level residual boundary separation', fontsize=fontsize)
plt.figure()
jellyfish(samples['probcond'])
plt.title('Condition-level resdiual probability of mind-wandering', fontsize=fontsize)


plt.figure()
jellyfish(samples['EEGcond'])
plt.title('Condition-level ERP latencies', fontsize=fontsize)

plt.figure()
jellyfish(samples['effectcond'])
plt.title('Condition-level ERP effects', fontsize=fontsize)

# # Calculate Bayes Factors
# possamps = samples['n200gammacond']
# # Number of chains
# nchains = possamps.shape[-1]
# # Number of samples per chain
# nsamps = possamps.shape[-2]
# # Number of dimensions
# ndims = possamps.ndim - 2
# # Number of variables to plot
# nvars = np.prod(possamps.shape[0:-2])
# # Index of variables
# varindx = np.arange(nvars).reshape(possamps.shape[0:-2])
# # Reshape data
# alldata = np.reshape(possamps, (nvars, nchains, nsamps))
# alldata = np.reshape(alldata, (nvars, nchains * nsamps))

# axes1 = plt.figure().add_subplot(111)
# jellyfish(samples['n200gammacond'][0])

# plt.plot([1, 1], [0, 7], 'r--', linewidth=2)
# plt.plot([0, 0], [0, 7], 'b--', linewidth=2)
# plt.xlim((-4,5))
# # plt.title('Condition-level n200 effects', fontsize=fontsize)
# yticks = axes1.get_yticks().tolist()
# yticks[0] = ''
# yticks[1] = 'R Hand, .6 sec'
# yticks[2] = 'R Hand, .9 sec'
# yticks[3] = 'R Hand, 1.5 sec'
# yticks[4] = 'L Hand, .6 sec'
# yticks[5] = 'L Hand, .9 sec'
# yticks[6] = 'L Hand, 1.5 sec'
# axes1.set_yticklabels(yticks)
# plt.xlabel(
#     'Effect of N200 latency on non-decision time (ms increase per ms increase in N200 latency)', fontsize=fontsize)
# axes1.tick_params(axis='both', which='major', labelsize=fontsize)
# # Calculate Bayes Factors using Savage-Dickey density ratio
# bf = dict()
# for k in range(0, 6):
#     kde = stats.gaussian_kde(alldata[k, :])
#     # Prior density of effect parameters
#     denom = stats.norm.pdf(1, loc=1, scale=3)
#     bf[k] = kde(1) / denom
#     plt.annotate('BF: %3.2f' % bf[k], (4.25, k + .95), fontsize=fontsize)
#     # Check kde0 integrates to one
#     # sum = 0
#     # for n in np.arange(-200, 200, .1):
#     #     sum = sum + kde(n)*.1
#     #
#     # (sum == 1)?
# fig = plt.gcf()
# fig.set_size_inches(20, 12)
# figManager = plt.get_current_fig_manager() #Maximize screen
# figManager.window.showMaximized()
# plt.savefig(saveloc + '/nondecisiontime_effects_random.png', dpi=300)

# plt.figure()
# jellyfish(samples['n200gammacond'][1])
# plt.title('Effect of N200 latency on drift-rate', fontsize=fontsize)
# plt.figure()
# jellyfish(samples['n200gammacond'][2])
# plt.title('Effect of N200 latency on boundary separation',
#           fontsize=fontsize)


# axes2 = plt.figure().add_subplot(111)
# ultdata = np.reshape(samples['n200gammault'][0,0,:],(nchains * nsamps))
# kde = stats.gaussian_kde(ultdata)
# bounds = stats.scoreatpercentile(ultdata, (.5, 2.5, 97.5, 99.5))
# # Plot properties
# LineWidths = np.array([2, 5])
# teal = np.array([0, .7, .7])
# blue = np.array([0, 0, 1])
# orange = np.array([1, .3, 0])
# Colors = [teal, blue]
# for b in range(0, 2):
#     # Bound by .5th percentile and 99.5th percentile
#     x = np.linspace(bounds[b], bounds[-1 - b], 100)
#     pdf = kde(x)
#     lines = plt.plot(x, pdf)
#     plt.setp(lines, color=Colors[b], linewidth=LineWidths[b])
#     if b == 1:
#         # Mark mode
#         wheremaxp = np.argmax(pdf)
#         mmode = plt.plot(np.array([0.5, 0.5]) * x[wheremaxp],
#                          np.array([pdf[wheremaxp], pdf[wheremaxp]]))
#         plt.setp(mmode, linewidth=3, color=orange)
#         # Mark median
#         mmedian = plt.plot(np.median(ultdata), 0.5, 'ko')
#         plt.setp(mmedian, markersize=10, color=[0., 0., 0.])
#         # Mark mean
#         mmean = plt.plot(np.mean(ultdata), 0.5, '*')
#         plt.setp(mmean, markersize=10, color=teal)
# plt.plot([1, 1], [0, 1], 'r--', linewidth=2)
# plt.plot([0, 0], [0, 1], 'b--', linewidth=2)
# plt.xlim((-4,5))
# # Prior density of effect parameters
# denom = stats.norm.pdf(1, loc=1, scale=3)
# bf[6] = kde(1) / denom
# plt.annotate('BF: %3.2f' % bf[6], (2.5, 0.5), fontsize=fontsize)

# fig = plt.gcf()
# plt.xlabel(
# 'Effect of N200 latency on non-decision time (ms increase per ms increase in N200 latency)', fontsize=fontsize)
# axes2.tick_params(axis='both', which='major', labelsize=fontsize)
# fig.set_size_inches(20, 12)
# figManager = plt.get_current_fig_manager() #Maximize screen
# figManager.window.showMaximized()
# plt.savefig(saveloc + '/nondecisiontime_expeffects_random.png', dpi=300)
