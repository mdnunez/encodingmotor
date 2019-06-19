# pdm7_jagsmodels.py - Contains preregistered JAGS models used in the analysis
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
# 10/25/18      Michael Nunez            Use of all_n1lat_random_lapse model
# 10/30/18      Michael Nunez                    behavior_only
# 11/01/18      Michael Nunez            all_EEG_lapse and ter_EEG_lapse models


# JAGS Models
jagsmodels = dict()


# Basic hierarchical drift-diffusion model (DDM) with lapse trials
jagsmodels['behavior_only'] = '''
model {
	##########
	#Fixed Parameters
	beta <- .5
	##########
	#Between-EEGsession variability in non-decision time
	tersubsd ~ dgamma(.2,1)

	#Between-EEGsession variability in drift
	deltasubsd ~ dgamma(1,1)

	#Between-EEGsession variability in boundary separation
	alphasubsd ~ dgamma(1,1)

	#Between-EEGsession variability in mind-wandering probability
	probsubsd ~ dgamma(.2,1)

	##########
	#Block-level parameters
	##########
	for (k in 1:nconds) {

	        #Condition-level non-decision time
			tercond[1,k] ~ dnorm(.3, pow(.25,-2))

			#Condition-level drift rate
			deltacond[1,k] ~ dnorm(1, pow(2, -2))

			#Condition-level boundary separation
			alphacond[1,k] ~ dnorm(1, pow(.5,-2))

			#Condition-level mind wandering probability
			probcond[1,k] ~ dnorm(.2, pow(.1,-2))

		#EEGsession-level parameters
		for (ses in 1:nses) {
			#EEGsession-level non-decision time
			tersub[k,ses] ~ dnorm(tercond[1,k],
				pow(tersubsd, -2))T(0,.7)

			#EEGsession-level drift rate
			deltasub[k,ses] ~ dnorm(deltacond[1,k],
				pow(deltasubsd, -2))T(-9, 9)

			#EEGsession-level boundary separation
			alphasub[k,ses] ~ dnorm(alphacond[1,k],
				pow(alphasubsd, -2))T(.1, 3)

			#EEGsession-level mind wandering probability
        	probsub[k, ses] ~ dnorm(probcond[1,k],
        		pow(probsubsd, -2))T(0, 1)
        	DDMprobsub[k, ses] <- 1 - probsub[k, ses]

		}
	}
	##########
	# Wiener likelihoods
	for (i in 1:N) {
        # Log density for DDM process
        ld_comp[i, 1] <- dlogwiener(y[i],alphasub[condition[i],EEGsession[i]],
        tersub[condition[i],EEGsession[i]],
        beta,
        deltasub[condition[i],EEGsession[i]])

        # Log density for lapse trials (negative max RT to positive max RT)
        ld_comp[i, 2] <- logdensity.unif(y[i], -maxrt[condition[i],EEGsession[i]], maxrt[condition[i],EEGsession[i]])

        # Select one of these two densities (Mixture of nonlapse and lapse trials)
        density[i] <- exp(ld_comp[i, component_chosen[i]] - Constant)
		
        # Generate a likelihood for the MCMC sampler using a trick to maximize density value
        Ones[i] ~ dbern(density[i])

        # Probability of mind wandering trials (lapse trials)
        component_chosen[i] ~ dcat( c(DDMprobsub[condition[i],EEGsession[i]], probsub[condition[i],EEGsession[i]]) )
	}
}
'''


# 4 Parameter Model with random effects of EEG latency on all parameters across subjects
#  Properly accounting for lapse trials
jagsmodels['all_EEG_lapse'] = '''
model {
	##########
	#Fixed Parameters
	beta <- .5
	##########
	#Between-EEGsession variability in non-decision time
	tersubsd ~ dgamma(.2,1)

	#Between-EEGsession variability in drift
	deltasubsd ~ dgamma(1,1)

	#Between-EEGsession variability in boundary separation
	alphasubsd ~ dgamma(1,1)

	#Between-EEGsession variability in mind-wandering probability
	probsubsd ~ dgamma(.2,1)

	#Between-EEGsession variability in EEG latency
	EEGsubsd ~ dgamma(.2,1)

	for (f in 1:4) {
		effectult[1,f] ~ dnorm(1,pow(3,-2)) #'Informative' prior
		effectsd[1,f] ~ dgamma(1,1)
	}

	##########
	#Block-level parameters
	##########
	for (k in 1:nconds) {

	        #Condition-level residual non-decision time
			tercond[1,k] ~ dnorm(.3, pow(.25,-2))

			#Condition-level residual drift rate
			deltacond[1,k] ~ dnorm(1, pow(2, -2))

			#Condition-level residual boundary separation
			alphacond[1,k] ~ dnorm(1, pow(.5,-2))

			#Condition-level residual mind wandering probability
			probcond[1,k] ~ dnorm(.2, pow(.1,-2))

			#Condition-level EEG latency
			EEGcond[1,k] ~ dnorm(.3, pow(.1,-2))

			#Condition-level effects of EEG input on non-decision time, drift rate, and boundary separation
			for (f in 1:4) {
				effectcond[f,k] ~ dnorm(effectult[1,f],pow(effectsd[1,f],-2))
			}
		#EEGsession-level parameters
		for (ses in 1:nses) {
			#EEGsession-level non-decision time
			tersub[k,ses] ~ dnorm(tercond[1,k]
				+ effectcond[1,k]*EEGsub[k,ses],
				pow(tersubsd, -2))T(0,.7)

			#EEGsession-level drift rate
			deltasub[k,ses] ~ dnorm(deltacond[1,k]
				+ effectcond[2,k]*EEGsub[k,ses],
				pow(deltasubsd, -2))T(-9, 9)

			#EEGsession-level boundary separation
			alphasub[k,ses] ~ dnorm(alphacond[1,k]
				+ effectcond[3,k]*EEGsub[k,ses],
				pow(alphasubsd, -2))T(.1,3)

			#EEGsession-level mind wandering probability
        	probsub[k, ses] ~ dnorm(probcond[1,k]
        		+ effectcond[4,k]*EEGsub[k,ses],
        		pow(probsubsd, -2))T(0, 1)
        	DDMprobsub[k, ses] <- 1 - probsub[k, ses]

        	#EEGsession-level EEG latency (impute missing values)
			EEGsub[k, ses] ~ dnorm(EEGcond[1,k],
				pow(EEGsubsd, -2))

		}
	}
	##########
	# Wiener likelihoods
	for (i in 1:N) {
        # Log density for DDM process
        ld_comp[i, 1] <- dlogwiener(y[i],alphasub[condition[i],EEGsession[i]],
        tersub[condition[i],EEGsession[i]],
        beta,
        deltasub[condition[i],EEGsession[i]])

        # Log density for lapse trials (negative max RT to positive max RT)
        ld_comp[i, 2] <- logdensity.unif(y[i], -maxrt[condition[i],EEGsession[i]], maxrt[condition[i],EEGsession[i]])

        # Select one of these two densities (Mixture of nonlapse and lapse trials)
        density[i] <- exp(ld_comp[i, component_chosen[i]] - Constant)
		
        # Generate a likelihood for the MCMC sampler using a trick to maximize density value
        Ones[i] ~ dbern(density[i])

        # Probability of mind wandering trials (lapse trials)
        component_chosen[i] ~ dcat( c(DDMprobsub[condition[i],EEGsession[i]], probsub[condition[i],EEGsession[i]]) )
	}
}
'''

# 4 Parameter Model with random effects of EEG latency on all parameters across subjects
#  Properly accounting for lapse trials
jagsmodels['ter_EEG_lapse'] = '''
model {
	##########
	#Fixed Parameters
	beta <- .5
	##########
	#Between-EEGsession variability in non-decision time
	tersubsd ~ dgamma(.2,1)

	#Between-EEGsession variability in drift
	deltasubsd ~ dgamma(1,1)

	#Between-EEGsession variability in boundary separation
	alphasubsd ~ dgamma(1,1)

	#Between-EEGsession variability in mind-wandering probability
	probsubsd ~ dgamma(.2,1)

	#Between-EEGsession variability in EEG latency
	EEGsubsd ~ dgamma(.2,1)

	effectult ~ dnorm(1,pow(3,-2)) #'Informative' prior
	effectsd ~ dgamma(1,1)

	##########
	#Block-level parameters
	##########
	for (k in 1:nconds) {

	        #Condition-level residual non-decision time
			tercond[1,k] ~ dnorm(.3, pow(.25,-2))

			#Condition-level drift rate
			deltacond[1,k] ~ dnorm(1, pow(2, -2))

			#Condition-level boundary separation
			alphacond[1,k] ~ dnorm(1, pow(.5,-2))

			#Condition-level mind wandering probability
			probcond[1,k] ~ dnorm(.2, pow(.1,-2))

			#Condition-level EEG latency
			EEGcond[1,k] ~ dnorm(.3, pow(.1,-2))

			#Condition-level effects of EEG input on non-decision time, drift rate, and boundary separation
			effectcond[1,k] ~ dnorm(effectult,pow(effectsd,-2))

		#EEGsession-level parameters
		for (ses in 1:nses) {
			#EEGsession-level non-decision time
			tersub[k,ses] ~ dnorm(tercond[1,k]
				+ effectcond[1,k]*EEGsub[k,ses],
				pow(tersubsd, -2))T(0,.7)

			#EEGsession-level drift rate
			deltasub[k,ses] ~ dnorm(deltacond[1,k],
				pow(deltasubsd, -2))T(-9, 9)

			#EEGsession-level boundary separation
			alphasub[k,ses] ~ dnorm(alphacond[1,k],
				pow(alphasubsd, -2))T(.1,3)

			#EEGsession-level mind wandering probability
        	probsub[k, ses] ~ dnorm(probcond[1,k],
        		pow(probsubsd, -2))T(0, 1)
        	DDMprobsub[k, ses] <- 1 - probsub[k, ses]

        	#EEGsession-level EEG latency (impute missing values)
			EEGsub[k, ses] ~ dnorm(EEGcond[1,k],
				pow(EEGsubsd, -2))

		}
	}
	##########
	# Wiener likelihoods
	for (i in 1:N) {
        # Log density for DDM process
        ld_comp[i, 1] <- dlogwiener(y[i],alphasub[condition[i],EEGsession[i]],
        tersub[condition[i],EEGsession[i]],
        beta,
        deltasub[condition[i],EEGsession[i]])

        # Log density for lapse trials (negative max RT to positive max RT)
        ld_comp[i, 2] <- logdensity.unif(y[i], -maxrt[condition[i],EEGsession[i]], maxrt[condition[i],EEGsession[i]])

        # Select one of these two densities (Mixture of nonlapse and lapse trials)
        density[i] <- exp(ld_comp[i, component_chosen[i]] - Constant)
		
        # Generate a likelihood for the MCMC sampler using a trick to maximize density value
        Ones[i] ~ dbern(density[i])

        # Probability of mind wandering trials (lapse trials)
        component_chosen[i] ~ dcat( c(DDMprobsub[condition[i],EEGsession[i]], probsub[condition[i],EEGsession[i]]) )
	}
}
'''

# 3 Parameter Model with random effects of minimum beta desynchronization on all parameters across sessions
jagsmodels['all_minBETA_random'] = '''
model {
	##########
	#Fixed Parameters
	beta <- .5
	##########
	#Between-session variability in non-decision time
	tersessd ~ dgamma(.2,1)

	#Between-session variability in drift
	deltasessd ~ dgamma(1,1)

	#Between-session variability in boundary separation
	alphasessd ~ dgamma(1,1)

	#Between-session variability in N200 latency
	betasessd ~ dgamma(.2,1)

	for (f in 1:3) {
		betagammault[1,f] ~ dnorm(1,pow(3,-2)) #'Informative' prior
		betagammasd[1,f] ~ dgamma(1,1)
	}

	##########
	#Block-level parameters
	##########
	for (k in 1:6) {
		#Condition-level N200 latency
		betacond[1,k] ~ dnorm(.2, pow(.1,-2))

        #Condition-level residual non-decision time
		tercond[1,k] ~ dnorm(.3, pow(.25,-2))

		#Condition-level residual drift rate
		deltacond[1,k] ~ dnorm(1, pow(2, -2))

		#Condition-level residual boundary separation
		alphacond[1,k] ~ dnorm(1, pow(.5,-2))

		#Condition-level effects of N200 latency on non-decision time, drift rate, and boundary separation
		for (f in 1:3) {
			betagammacond[f,k] ~ dnorm(betagammault[1,f],pow(betagammasd[1,f],-2))
		}

		#session-level parameters
		for (ses in 1:nses) {
			#session-level non-decision time
			terses[k,ses] ~ dnorm(tercond[1,k]
				+ betagammacond[1,k]*betases[k,ses],
				pow(tersessd, -2))T(0,1)

			#session-level drift rate
			deltases[k,ses] ~ dnorm(deltacond[1,k]
				+ betagammacond[2,k]*betases[k,ses],
				pow(deltasessd, -2))T(-9, 9)

			#session-level boundary separation
			alphases[k,ses] ~ dnorm(alphacond[1,k]
				+ betagammacond[3,k]*betases[k,ses],
				pow(alphasessd, -2))T(.1,3)

		    #session-level N200 latency
			betases[k,ses] ~ dnorm(betacond[1,k],
				pow(betasessd, -2))
		}
	}
	##########
	#Wiener likelihoods
	for (i in 1:N) {
		y[i] ~ dwiener(alphases[condition[i],session[i]],
		terses[condition[i],session[i]],
		beta,
		deltases[condition[i],session[i]])
	}
}
'''

# 3 Parameter Model with random effects of trial-averaged N200 on all parameters across sessions, additive effects of hand and time-to-respond
jagsmodels['all_n200lat_additive'] = '''
model {
	##########
	#Fixed Parameters
	beta <- .5
	##########
	#Between-session variability in non-decision time
	tersessd ~ dgamma(.2,1)

	#Between-session variability in drift
	deltasessd ~ dgamma(1,1)

	#Between-session variability in boundary separation
	alphasessd ~ dgamma(1,1)

	#Between-session variability in N200 latency
	n200sessd ~ dgamma(.2,1)

	##########
	#Block-level parameters
	##########

	#Condition-level effects of N200 latency on non-decision time
	for (f in 1:3) {
		n200gammacond[f,1] ~ dnorm(1,pow(3,-2))
		for (k in 2:6) {
			n200gammacond[f,k] ~ dnorm(0,pow(1,-2))
		}
	}

	for (k in 1:6) {
		#Condition-level N200 latency
		n200cond[1,k] ~ dnorm(.2, pow(.1,-2))

        #Condition-level residual non-decision time
		tercond[1,k] ~ dnorm(.3, pow(.25,-2))

		#Condition-level residual drift rate
		deltacond[1,k] ~ dnorm(1, pow(2, -2))

		#Condition-level residual boundary separation
		alphacond[1,k] ~ dnorm(1, pow(.5,-2))

		#session-level parameters
		for (ses in 1:nses) {
			#session-level non-decision time
			terses[k,ses] ~ dnorm(tercond[1,k]
				+ n200gammacond[1,1]*n200ses[k,ses]
				+ n200gammacond[1,2]*((k==2) || (k==3) || (k>=5))*n200ses[k,ses]
				+ n200gammacond[1,3]*((k==3) || (k==6))*n200ses[k,ses]
				+ n200gammacond[1,4]*(k<4)*n200ses[k,ses]
				+ n200gammacond[1,5]*(k<4)*((k==2) || (k==3) || (k>=5))*n200ses[k,ses]
				+ n200gammacond[1,6]*(k<4)*((k==3) || (k==6))*n200ses[k,ses],
				pow(tersessd, -2))T(0,1)

			#session-level drift rate
			deltases[k,ses] ~ dnorm(deltacond[experiment[ses]+1,k]
				+ n200gammacond[2,1]*n200ses[k,ses]
				+ n200gammacond[2,2]*((k==2) || (k==3) || (k>=5))*n200ses[k,ses]
				+ n200gammacond[2,3]*((k==3) || (k==6))*n200ses[k,ses]
				+ n200gammacond[2,4]*(k<4)*n200ses[k,ses]
				+ n200gammacond[2,5]*(k<4)*((k==2) || (k==3) || (k>=5))*n200ses[k,ses]
				+ n200gammacond[2,6]*(k<4)*((k==3) || (k==6))*n200ses[k,ses],
				pow(deltasessd, -2))T(-9, 9)

			#session-level boundary separation
			alphases[k,ses] ~ dnorm(alphacond[experiment[ses]+1,k]
				+ n200gammacond[3,1]*n200ses[k,ses]
				+ n200gammacond[3,2]*((k==2) || (k==3) || (k>=5))*n200ses[k,ses]
				+ n200gammacond[3,3]*((k==3) || (k==6))*n200ses[k,ses]
				+ n200gammacond[3,4]*(k<4)*n200ses[k,ses]
				+ n200gammacond[3,5]*(k<4)*((k==2) || (k==3) || (k>=5))*n200ses[k,ses]
				+ n200gammacond[3,6]*(k<4)*((k==3) || (k==6))*n200ses[k,ses],
				pow(alphasessd, -2))T(.1,3)

		    #session-level N200 latency
			n200ses[k,ses] ~ dnorm(n200cond[1,k],
				pow(n200sessd, -2))
		}
	}
	##########
	#Wiener likelihoods
	for (i in 1:N) {
		y[i] ~ dwiener(alphases[condition[i],session[i]],
		terses[condition[i],session[i]],
		beta,
		deltases[condition[i],session[i]])
	}
}
'''

# 3 Parameter Model with random effects of trial-averaged N200 and minimum beta desynchronization on all parameters across sessions
jagsmodels['all_n200lat_minBETA_random'] = '''
model {
	##########
	#Fixed Parameters
	beta <- .5
	##########
	#Between-session variability in non-decision time
	tersessd ~ dgamma(.2,1)

	#Between-session variability in drift
	deltasessd ~ dgamma(1,1)

	#Between-session variability in boundary separation
	alphasessd ~ dgamma(1,1)

	#Between-session variability in N200 latency
	n200sessd ~ dgamma(.2,1)

	#Between-session variability in minimum beta desynchronization time
	betasessd ~ dgamma(.2,1)

	for (f in 1:3) {
		n200gammault[1,f] ~ dnorm(1,pow(3,-2)) #'Informative' prior
		n200gammasd[1,f] ~ dgamma(1,1)
		betagammault[1,f] ~ dnorm(1,pow(3,-2)) #'Informative' prior
		betagammasd[1,f] ~ dgamma(1,1)
	}

	##########
	#Block-level parameters
	##########
	for (k in 1:6) {
		#Condition-level N200 latency
		n200cond[1,k] ~ dnorm(.2, pow(.1,-2))

		#Condition-level minimum beta desynchronization time
		betacond[1,k] ~ dnorm(.4, pow(.1,-2))

        #Condition-level residual non-decision time
		tercond[1,k] ~ dnorm(.3, pow(.25,-2))

		#Condition-level residual drift rate
		deltacond[1,k] ~ dnorm(1, pow(2, -2))

		#Condition-level residual boundary separation
		alphacond[1,k] ~ dnorm(1, pow(.5,-2))

		#Condition-level effects of N200 latency on non-decision time, drift rate, and boundary separation
		for (f in 1:3) {
			n200gammacond[f,k] ~ dnorm(n200gammault[1,f],pow(n200gammasd[1,f],-2))
			betagammacond[f,k] ~ dnorm(betagammault[1,f],pow(betagammasd[1,f],-2))
		}

		#session-level parameters
		for (ses in 1:nses) {
			#session-level non-decision time
			terses[k,ses] ~ dnorm(tercond[1,k]
				+ n200gammacond[1,k]*n200ses[k,ses]
				+ betagammacond[1,k]*betases[k,ses],
				pow(tersessd, -2))T(0,1)

			#session-level drift rate
			deltases[k,ses] ~ dnorm(deltacond[1,k]
				+ n200gammacond[2,k]*n200ses[k,ses]
				+ betagammacond[2,k]*betases[k,ses],
				pow(deltasessd, -2))T(-9, 9)

			#session-level boundary separation
			alphases[k,ses] ~ dnorm(alphacond[1,k]
				+ n200gammacond[3,k]*n200ses[k,ses]
				+ betagammacond[3,k]*betases[k,ses],
				pow(alphasessd, -2))T(.1,3)

		    #session-level N200 latency
			n200ses[k,ses] ~ dnorm(n200cond[1,k],
				pow(n200sessd, -2))

			#session-level minimum beta desynchronization time
			betases[k,ses] ~ dnorm(betacond[1,k],
				pow(betasessd, -2))
		}
	}
	##########
	#Wiener likelihoods
	for (i in 1:N) {
		y[i] ~ dwiener(alphases[condition[i],session[i]],
		terses[condition[i],session[i]],
		beta,
		deltases[condition[i],session[i]])
	}
}
'''
