# pdm7_plotERPs.py - Plots N200, P300 and RP waveforms split by condition
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
# 11/02/18     Michael Nunez               Coverted from pdm5b_plotN200_split.py


## References:
# http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=6
# https://matplotlib.org/gallery/text_labels_and_annotations/custom_legends.html

# Imports
from __future__ import division
import numpy as np
from scipy import stats
import scipy.io as sio
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from IPython import get_ipython  # Run magic functions from script
# Initialize ipython matplotlib plotting graphics
get_ipython().magic('pylab')

# Parameters

fontsize = 16
fontsize2 = 20
color1 = "#D95F02"
color2 = "#1B9E77"
color3 = "#7570B3"
color4 = "#E6AB02"
color5 = "#66A61E"
color6 = "#E7298A"
color7 = "#A6761D"
colororder = (color1, color2, color3,color4, color5, color6)

# Initial
loadloc = '/home/michael/data10/michael/pdm/exp7data/subjects/allerps.mat'


# Load all data
print 'Loading all ERP data...'
erps = sio.loadmat(loadloc)

# Load condition-level data
sesconddata = np.genfromtxt('../sesconddata.csv', delimiter=',')    # file name
n200latvec = sesconddata[1:, 8]
meanRT = sesconddata[1:, 0]
accuracy = sesconddata[1:, 2]
condition = sesconddata[1:, 3]
condRT = np.empty(6)

# Plot all N200 waveforms
print 'Plotting the N200 data'
plt.figure(figsize=(12, 6))
for n in range(0,erps['N200timecourse'].shape[1]):
    plt.plot(erps['N200window'][0,:], np.squeeze(erps['N200timecourse'][:,n,:]),color=colororder[n])
zeroline = plt.plot(np.array([0., 0.]), np.array([-.15, .15]))
plt.setp(zeroline, linewidth=3, linestyle=':',color='k')
zeroline2 = plt.plot(np.array([-50., 275.]), np.array([0, 0]))
plt.setp(zeroline2, linewidth=3, linestyle=':',color='k')
plt.xlim(-50, 275)
plt.ylim(-.15, .15)

# Plot N200 deflection distributions
kde_n200lat = stats.gaussian_kde(n200latvec*1000)
x_n200lat = np.linspace(125, 275, 100)
p_n200lat = kde_n200lat(x_n200lat)
p_n200lat = .15 * p_n200lat / np.max(p_n200lat)
plt.fill_between(x_n200lat, p_n200lat, np.zeros((100)), color=color7, alpha=0.25)
plt.plot(x_n200lat, p_n200lat, color='k', alpha=0.25)

# Format plot
frame = plt.gca()
frame.axes.get_yaxis().set_visible(False)
plt.tick_params(axis='both', which='major', labelsize=fontsize)
plt.xlabel('Time (ms) after Gabor onset', fontsize=fontsize)

# Custom legend
custom_lines = [Line2D([0], [0], color=color1, lw=4),
                Line2D([0], [0], color=color2, lw=4),
                Line2D([0], [0], color=color3, lw=4),
                Line2D([0], [0], color=color4, lw=4),
                Line2D([0], [0], color=color5, lw=4),
                Line2D([0], [0], color=color6, lw=4)]
plt.legend(custom_lines, ['Right .6 s', 'Right .9 s', 'Right 1.5 s',
    'Left .6 s', 'Left .9 s', 'Left 1.5 s'],loc=2)

# Save the figure
print 'Saving the figure...'
plt.savefig('N200_waveforms.png', dpi=300, format='png')

# Plot all N200 & P300 average waveforms
print 'Plotting the average N200 and P300 waveforms'
plt.figure(figsize=(18, 6))
for n in range(0,erps['RPtimecourse'].shape[1]):
    plt.plot(erps['N200window'][0,:], np.mean(np.squeeze(erps['N200timecourse'][:,n,:]),axis=1),
        color=colororder[n],linewidth=3,linestyle='--')
    plt.plot(erps['P300window'][0,:], 3*np.mean(np.squeeze(erps['P300timecourse'][:,n,:]),axis=1),
        color=colororder[n],linewidth=3)
    condRT[n] = np.mean(meanRT[(condition == n+1)]) 
    plt.plot(np.array([condRT[n], condRT[n]])*1000, np.array([-.15, .15]),
        color=colororder[n],linewidth=5,linestyle=':')
zeroline = plt.plot(np.array([0., 0.]), np.array([-.15, .15]))
plt.setp(zeroline, linewidth=5, linestyle=':',color='k')
zeroline2 = plt.plot(np.array([-50., 800.]), np.array([0, 0]))
plt.setp(zeroline2, linewidth=5, linestyle=':',color='k')
plt.xlim(-50, 800)
plt.ylim(-.15, .15)

# Format plot
frame = plt.gca()
frame.axes.get_yaxis().set_visible(False)
plt.tick_params(axis='both', which='major', labelsize=fontsize2)
plt.xlabel('Time (ms) after Gabor onset', fontsize=fontsize2)

# Custom legend
custom_lines = [Line2D([0], [0], color=color1, lw=4),
                Line2D([0], [0], color=color2, lw=4),
                Line2D([0], [0], color=color3, lw=4),
                Line2D([0], [0], color=color4, lw=4),
                Line2D([0], [0], color=color5, lw=4),
                Line2D([0], [0], color=color6, lw=4),
                Line2D([0], [0], color='k', lw=4,linestyle='--'),
                Line2D([0], [0], color='k', lw=4),
                Line2D([0], [0], color='k', lw=4,linestyle=':')]
plt.legend(custom_lines, ['Right .6 s', 'Right .9 s', 'Right 1.5 s',
    'Left .6 s', 'Left .9 s', 'Left 1.5 s','N200 waveform','P300 waveform','Mean RT'],
    fontsize=fontsize, loc=2, ncol=2, frameon=False)

# Save the figure
print 'Saving the figure...'
plt.savefig('Stimlocked_waveforms.png', dpi=300, format='png')

# Plot all RP and beta desynchronization waveforms
print 'Plotting the average RP waveform and beta desynchronization waveforms'
plt.figure(figsize=(15, 6))
for n in range(0,erps['RPtimecourse'].shape[1]):
    plt.plot(erps['RPwindow'][0,:], np.mean(np.squeeze(erps['RPtimecourse'][:,n,:]),axis=1),
        color=colororder[n],linewidth=3)
    plt.plot(erps['BETAwindow'][0,:], .1*np.mean(np.squeeze(erps['BETAtimecourse'][:,n,:])-.65,axis=1),
        color=colororder[n],linewidth=3,linestyle='--')
    condRT[n] = np.mean(meanRT[(condition == n+1)]) 
    plt.plot(np.array([condRT[n], condRT[n]])*-1000, np.array([-.1, .025]),
        color=colororder[n],linewidth=5,linestyle=':')
zeroline = plt.plot(np.array([0., 0.]), np.array([-.1, .025]))
plt.setp(zeroline, linewidth=5, linestyle=':',color='k')
zeroline2 = plt.plot(np.array([-700., 50.]), np.array([0, 0]))
plt.setp(zeroline2, linewidth=5, linestyle=':',color='k')
plt.xlim(-700, 50)
plt.ylim(-.1, .025)

# Format plot
frame = plt.gca()
frame.axes.get_yaxis().set_visible(False)
plt.tick_params(axis='both', which='major', labelsize=fontsize2)
plt.xlabel('Time (ms) before response', fontsize=fontsize2)

# Custom legend
custom_lines = [Line2D([0], [0], color='k', lw=4),
                Line2D([0], [0], color='k', lw=4,linestyle='--'),
                Line2D([0], [0], color='k', lw=4,linestyle=':')]
plt.legend(custom_lines, ['Readiness Potential (RP)', 'Beta (14-22 Hz) at C3 (right hand) or C4 (left hand)', 'Mean Gabor Onset (-RT)'],
    fontsize=fontsize2, loc=4, frameon=False)

# Save the figure
print 'Saving the figure...'
plt.savefig('Resplocked_waveforms.png', dpi=300, format='png')