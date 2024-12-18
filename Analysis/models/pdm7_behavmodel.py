# pdm7_behavmodel.py - Runs JAGS model 'behavior_only'
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
# 10/29/18      Michael Nunez                             Original code
# 10/30/18      Michael Nunez               Condition-level mind wandering parameters

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

# Initialize non-decision time with mininum RT for each subject and condition
# Use maximum RT for the bounds on the lapse process, modeled by a uniform distribution
minrt = np.empty((nconds, nses))
maxrt = np.empty((nconds, nses))
for k in range(0, nconds):
    for j in range(0, nses):
        where = (session == (j + 1)) & (condition == (k + 1))
        minrt[k, j] = np.min(np.abs(y[where]))
        maxrt[k, j] = np.max(np.abs(y[where]))


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

# # Track these variables
# trackvars1 = ['alphasub', 'deltasub', 'tersub',
#              'alphacond', 'deltacond', 'tercond',
#              'alphasubsd', 'deltasubsd', 'tersubsd','probsub']


# initials = []
# for c in range(0, nchains):
#     chaininit = {
#         'alphasub': np.random.uniform(.5, 2., size=(nconds, nses)),
#         'deltasub': np.random.uniform(-4., 4., size=(nconds, nses)),
#         'tersub': np.random.uniform(0., .3, size=(nconds, nses)),
#         'alphacond': np.random.uniform(.5, 2., size=(1, nconds)),
#         'deltacond': np.random.uniform(-4., 4., size=(1, nconds)),
#         'tercond': np.random.uniform(0., .3, size=(1, nconds)),
#         'alphasubsd': np.random.uniform(.01, 1.),
#         'deltasubsd': np.random.uniform(.01, 3.),
#         'tersubsd': np.random.uniform(.01, .1),
#         'probsub': np.random.dirichlet((1,1),size=(nconds, nses))
#     }
#     for k in range(0, nconds):
#         for j in range(0, nses):
#             chaininit['tersub'][k, j] = np.random.uniform(0., minrt[k, j]/2)
#     initials.append(chaininit)

# Track these variables
trackvars1 = ['alphasub', 'deltasub', 'tersub', 'probsub',
             'alphacond', 'deltacond', 'tercond', 'probcond',
             'alphasubsd', 'deltasubsd', 'tersubsd','probsubsd', 'component_chosen']

initials = []
for c in range(0, nchains):
    chaininit = {
        'alphasub': np.random.uniform(.5, 2., size=(nconds, nses)),
        'deltasub': np.random.uniform(-4., 4., size=(nconds, nses)),
        'tersub': np.random.uniform(0., .3, size=(nconds, nses)),
        'probsub': np.random.uniform(0., 1., size=(nconds, nses)),
        'alphacond': np.random.uniform(.5, 2., size=(1, nconds)),
        'deltacond': np.random.uniform(-4., 4., size=(1, nconds)),
        'tercond': np.random.uniform(0., .3, size=(1, nconds)),
        'probcond': np.random.uniform(0., 1., size=(1, nconds)),
        'alphasubsd': np.random.uniform(.01, 1.),
        'deltasubsd': np.random.uniform(.01, 3.),
        'tersubsd': np.random.uniform(.01, .1),
        'probsubsd': np.random.uniform(.01, .1)
    }
    for k in range(0, nconds):
        for j in range(0, nses):
            chaininit['tersub'][k, j] = np.random.uniform(0., minrt[k, j]/2)
    initials.append(chaininit)


# Run JAGS model

# Choose JAGS model type
modelname = 'behavior_only'

thismodel = jagsmodels[modelname]

# Save model
timestart = strftime('%b') + '_' + strftime('%d') + '_' + \
    strftime('%y') + '_' + strftime('%H') + '_' + strftime('%M')
modelfile = saveloc + '/jagsmodels/model_' + timestart + '.jags'
f = open(modelfile, 'w')
f.write(thismodel)
f.close()
print 'Fitting model %s ...' % (modelname + timestart)

threaded = pyjags.Model(file=modelfile, init=initials,
                        data=dict(N=N, y=y, nses=nses, nconds=nconds, EEGsession=session,
                                  condition=condition, Constant=Constant, Ones=Ones, maxrt=maxrt),
                        chains=nchains, adapt=burnin, threads=6,
                        progress_bar=True)


samples = threaded.sample(nsamps, vars=trackvars1, thin=10)

savestring = saveloc + '/jagsout/model_' + \
    modelname + timestart + ".mat"

print 'Saving results to: \n %s' % savestring

sio.savemat(savestring, samples)

