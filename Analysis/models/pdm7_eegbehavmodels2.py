# pdm7_eegbehavmodels2.py - Runs JAGS models of type 'ter_EEG_lapse'
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
# 08/17/17      Michael Nunez                             Original code
# 10/31/18      Michael Nunez               Edits for flexible model building
# 11/02/18      Michael Nunez              Fix data indexing

# Imports
from __future__ import division
import numpy as np
import numpy.ma as ma
import scipy.io as sio
import pyjags
from scipy import stats
from time import strftime
import os
from pdm7_jagsmodels import *

# Flag
EEGmeasure = 'n200sub'
# EEGmeasure = 'betasub'
# EEGmeasure = 'p300sub'
# EEGmeasure = 'RPsub'


# Save location
thisdir = os.getcwd()
os.chdir('..')
os.chdir('..')
os.chdir('Data')
saveloc = os.getcwd()

# Load trial data
trialdata = np.genfromtxt('trialdata.csv', delimiter=',')    # file name

# Find data to remove
testdata = np.logical_not(np.array(trialdata[1:, 17], dtype=bool))
artifact = np.array(trialdata[1:, 18], dtype=bool)
allremove = np.logical_or(testdata,artifact)
allremove = np.logical_or(allremove,np.isnan(trialdata[1:, 0]))
# removeRT = np.array(trialdata[1:, 16], dtype=bool)
# allremove = np.logical_or(testdata, removeRT)
allkeep = np.logical_not(allremove)

y = trialdata[1:, 0] * (trialdata[1:, 1] * 2 - 1)
y = y[allkeep]
N = y.shape[0]
condition = trialdata[1:, 2]
condition = condition[allkeep]
session = trialdata[1:, 5]
session = session[allkeep]
nses = np.unique(session).shape[0]
nconds = 6

# Load condition-level data
sesconddata = np.genfromtxt('sesconddata.csv', delimiter=',')    # file name
smallcond = sesconddata[1:, 3]
smallses = sesconddata[1:, 6]
n200latvec = sesconddata[1:, 8]
rmn200vec = np.array(sesconddata[1:, 9], dtype=bool)
minBETAvec = sesconddata[1:, 11]
rmBETAvec = np.array(sesconddata[1:, 12], dtype=bool)
p300latvec = sesconddata[1:, 14]
rmp300vec = np.array(sesconddata[1:, 15], dtype=bool)
RPlatvec = sesconddata[1:, 17]
rmRPvec = np.array(sesconddata[1:, 18], dtype=bool)

minrt = np.empty((nconds, nses))
maxrt = np.empty((nconds, nses))
n200sub = np.empty((nconds, nses))
betasub = np.empty((nconds, nses))
p300sub = np.empty((nconds, nses))
RPsub = np.empty((nconds, nses))

for k in range(0, nconds):
    for j in range(0, nses):
        trialwhere = (session == (j + 1)) & (condition == (k + 1))
        sescondwhere = (smallses == (j + 1)) & (smallcond == (k + 1))
        n200sub[k, j] = n200latvec[sescondwhere]
        if rmn200vec[sescondwhere]:
            n200sub[k, j] = np.nan
        betasub[k, j] = minBETAvec[sescondwhere]
        if rmBETAvec[sescondwhere]:
            betasub[k, j] = np.nan
        p300sub[k, j] = p300latvec[sescondwhere]
        if rmp300vec[sescondwhere]:
            p300sub[k, j] = np.nan
        RPsub[k, j] = RPlatvec[sescondwhere]
        if rmRPvec[sescondwhere]:
            RPsub[k, j] = np.nan
        # Initialize non-decision time with mininum RT for each subject and
        # condition
        minrt[k, j] = np.min(np.abs(y[trialwhere]))
        maxrt[k, j] = np.max(np.abs(y[trialwhere]))
n200sub = ma.masked_invalid(n200sub)  # Remove NANs for pyjags
betasub = ma.masked_invalid(betasub)  # Remove NANs for pyjags
p300sub = ma.masked_invalid(p300sub)  # Remove NANs for pyjags
RPsub = ma.masked_invalid(RPsub)  # Remove NANs for pyjags

# Input for mixture modeling
Ones = np.ones(N)
Constant = 10

# pyjags code

# Make sure $LD_LIBRARY_PATH sees /usr/local/lib
pyjags.modules.load_module('wiener')
pyjags.modules.load_module('dic')
pyjags.modules.list_modules()

nchains = 6
burnin = 2000  # Note that scientific notation breaks pyjags
nsamps = 50000

# Track these variables
trackvars1 = ['alphasub', 'deltasub', 'tersub', 'probsub',
             'alphacond', 'deltacond', 'tercond', 'probcond', 'EEGcond', 'effectcond',
             'alphasubsd', 'deltasubsd', 'tersubsd','probsubsd', 'EEGsubsd', 'effectsd',
             'effectult']

initials = []
for c in range(0, nchains):
    chaininit = {
        'alphasub': np.random.uniform(.5, 2., size=(nconds, nses)),
        'deltasub': np.random.uniform(-4., 4., size=(nconds, nses)),
        'tersub': np.random.uniform(0., .3, size=(nconds, nses)),
        'probsub': np.random.uniform(0., 1., size=(nconds, nses)),
        'effectcond': np.random.uniform(-2., 2., size=(1, nconds)),
        'alphacond': np.random.uniform(.5, 2., size=(1, nconds)),
        'deltacond': np.random.uniform(-4., 4., size=(1, nconds)),
        'tercond': np.random.uniform(0., .3, size=(1, nconds)),
        'probcond': np.random.uniform(0., 1., size=(1, nconds)),
        'EEGsubsd': np.random.uniform(.01, .1),
        'alphasubsd': np.random.uniform(.01, 1.),
        'deltasubsd': np.random.uniform(.01, 3.),
        'tersubsd': np.random.uniform(.01, .1),
        'probsubsd': np.random.uniform(.01, .1),
        'effectsd': np.random.uniform(.01, 3.),
        'effectult': np.random.uniform(-2., 2.)
    }
    for k in range(0, nconds):
        for j in range(0, nses):
            chaininit['tersub'][k, j] = np.random.uniform(0., minrt[k, j]/2)
    initials.append(chaininit)

# Run JAGS model

# Choose JAGS model type
modelname = 'ter_EEG_lapse'

thismodel = jagsmodels[modelname]

# Save model
timestart = strftime('%b') + '_' + strftime('%d') + '_' + \
    strftime('%y') + '_' + strftime('%H') + '_' + strftime('%M')
modelfile = saveloc + '/jagsmodels/model_' + EEGmeasure + modelname + timestart + '.jags'
f = open(modelfile, 'w')
f.write(thismodel)
f.close()
print 'Fitting model %s ...' % (EEGmeasure + modelname + timestart)

exec('thismeasure = %s' % EEGmeasure)
threaded = pyjags.Model(file=modelfile, init=initials,
                        data=dict(N=N, y=y, nses=nses, nconds=nconds, EEGsession=session,
                                EEGsub=thismeasure, condition=condition, Constant=Constant, Ones=Ones, maxrt=maxrt),
                        chains=nchains, adapt=burnin, threads=6,
                        progress_bar=True)


samples = threaded.sample(nsamps, vars=trackvars1, thin=10)

savestring = saveloc + '/jagsout/model_' + \
    EEGmeasure + modelname + timestart + ".mat"

print 'Saving results to: \n %s' % savestring

sio.savemat(savestring, samples)