# pdm7_posteriordistributions1.R - Plots overlapping posterior distributions from Model 1
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

## Record of Revisions
#   Date           Programmers               Description of change
#   ====        =================            =====================
#  11/02/18        Michael Nunez     Converted from pdm5b_posteriordistributions1.R

## Necessary packages
library(ggplot2)
library(ggjoy)
library(viridis)
library(R.matlab)

loadloc = '../jagsout'

## Code

# Read in the reaction times
# samples = readMat(paste(loadloc,
#   '/model_n200suball_EEG_lapseOct_31_18_18_30.mat',sep=""))
samples = readMat(paste(loadloc,
  '/model_p300suball_EEG_lapseNov_03_18_00_38.mat',sep=""))

maineffect = as.vector(samples[9]$effectult[1,2,,])
print(sprintf('Length of P300 effect posterior samples is %d',length(maineffect)))
maineffectCI = quantile(maineffect, prob=c(.025,.5,.975))
print(sprintf('Effect (median posterior and 95%% credible interval) of trial-averaged P300 latency on drift rate: %.2f, CI: [%.2f, %.2f]',maineffectCI[2],maineffectCI[1],maineffectCI[3]));
maineffect_density = density(maineffect)
x11()
plot(maineffect_density)


nsamps = length(maineffect)
alleffects <- data.frame(matrix(ncol = 2, nrow = nsamps*7))
colnames(alleffects) <- c("Effect", "Condition")

alleffects$Effect[1:nsamps] = as.vector(samples[9]$effectult[1,3,,])
alleffects$Condition[1:nsamps] = sprintf('G')

alleffects$Effect[(nsamps+1):(nsamps*2)] = as.vector(samples[11]$effectcond[3,6,,])
alleffects$Condition[(nsamps+1):(nsamps*2)] = sprintf('A')

alleffects$Effect[(nsamps*2+1):(nsamps*3)] = as.vector(samples[11]$effectcond[3,5,,])
alleffects$Condition[(nsamps*2+1):(nsamps*3)] = sprintf('B')

alleffects$Effect[(nsamps*3+1):(nsamps*4)] = as.vector(samples[11]$effectcond[3,4,,])
alleffects$Condition[(nsamps*3+1):(nsamps*4)] = sprintf('C')

alleffects$Effect[(nsamps*4+1):(nsamps*5)] = as.vector(samples[11]$effectcond[3,3,,])
alleffects$Condition[(nsamps*4+1):(nsamps*5)] = sprintf('D')

alleffects$Effect[(nsamps*5+1):(nsamps*6)] = as.vector(samples[11]$effectcond[3,2,,])
alleffects$Condition[(nsamps*5+1):(nsamps*6)] = sprintf('E')

alleffects$Effect[(nsamps*6+1):(nsamps*7)] = as.vector(samples[11]$effectcond[3,1,,])
alleffects$Condition[(nsamps*6+1):(nsamps*7)] = sprintf('F')

# N200intercepts <- data.frame(matrix(ncol = 3, nrow = nsamps*6))
# colnames(N200intercepts) <- c("Intercept", "N200Lat", "Condition")

# N200intercepts$Intercept[1:nsamps] = as.vector(samples[2]$tercond[2,1,,])*1000
# N200intercepts$Condition[1:nsamps] = sprintf('A')
# N200intercepts$Intercept[(nsamps+1):(nsamps*2)] = as.vector(samples[2]$tercond[2,2,,])*1000
# N200intercepts$Condition[(nsamps+1):(nsamps*2)] = sprintf('B')
# N200intercepts$Intercept[(nsamps*2+1):(nsamps*3)] = as.vector(samples[2]$tercond[2,3,,])*1000
# N200intercepts$Condition[(nsamps*2+1):(nsamps*3)] = sprintf('C')
# N200intercepts$Intercept[(nsamps*3+1):(nsamps*4)] = as.vector(samples[2]$tercond[1,1,,])*1000
# N200intercepts$Condition[(nsamps*3+1):(nsamps*4)] = sprintf('D')
# N200intercepts$Intercept[(nsamps*4+1):(nsamps*5)] = as.vector(samples[2]$tercond[1,2,,])*1000
# N200intercepts$Condition[(nsamps*4+1):(nsamps*5)] = sprintf('E')
# N200intercepts$Intercept[(nsamps*5+1):(nsamps*6)] = as.vector(samples[2]$tercond[1,3,,])*1000
# N200intercepts$Condition[(nsamps*5+1):(nsamps*6)] = sprintf('F')

# N200intercepts$N200Lat[1:nsamps] = as.vector(samples[8]$n1cond[2,1,,])*1000
# N200intercepts$N200Lat[(nsamps+1):(nsamps*2)] = as.vector(samples[8]$n1cond[2,2,,])*1000
# N200intercepts$N200Lat[(nsamps*2+1):(nsamps*3)] = as.vector(samples[8]$n1cond[2,3,,])*1000
# N200intercepts$N200Lat[(nsamps*3+1):(nsamps*4)] = as.vector(samples[8]$n1cond[1,1,,])*1000
# N200intercepts$N200Lat[(nsamps*4+1):(nsamps*5)] = as.vector(samples[8]$n1cond[1,2,,])*1000
# N200intercepts$N200Lat[(nsamps*5+1):(nsamps*6)] = as.vector(samples[8]$n1cond[1,3,,])*1000

condlabs = c("Left 1.5 s", "Left .9 s", "Left .6 s", "Right 1.5 s", "Right .9 s", "Right .6 s", "Overall (hierarchical) effect") 

# http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=7
cbbPalette <- c("#000000", "#E7298A", "#66A61E", "#E6AB02", "#7570B3", "#1B9E77", "#D95F02", "#A6761D")

# + xlim(-2,4) + 
## Plot
png(paste('Model1_alleffects.png',sep=""),units="in",width=10,height=10,res=300)
plot1 = ggplot(alleffects,aes(x=Effect,y=Condition,fill=Condition)) +
  geom_joy(scale=3, alpha=.67) +
  xlab('Effect of trial-averaged P300 latency on boundary \n (evidence unit increase per sec increase in P300)') +
  ylab('') +
  # ggtitle('ERP latency effect posterior distributions') +
  annotate("segment", x=0, xend=0, y=1, yend=10, colour = 'red',size=3,alpha=.5) +
  theme(text = element_text(size=20), legend.position="none") +
  scale_fill_manual(values = tail(cbbPalette,n=7)) +
  scale_y_discrete(labels=condlabs)
# #Print Bayes Factors
# BF = vector()
# for (n in seq(1,7)) {
# 	temp_density = density(alleffects$Effect[(nsamps*(n-1)+1):(nsamps*n)])
# 	numerator = approx(temp_density$x,temp_density$y,xout=1)
# 	denominator = dnorm(1,mean=1,sd=3)
# 	BF[n] = numerator$y / denominator
# 	if (n==1) {
# 		yplace = 7.5
# 	} else {
# 		yplace = n - 0.5
# 	}
# 	plot1 = plot1 + annotate("text",x=3.5,y=yplace,label=sprintf('BF1: %3.2f',BF[n]),size=5,colour = 'blue')
# }
plot(plot1)
dev.off()


# png(paste('Model1_ResidualNDT.png',sep=""),units="in",width=10,height=10,res=300)
# plot2 = ggplot(N200intercepts,aes(x=Intercept,y=Condition,fill=Condition)) +
#   geom_joy(scale=3, alpha=.67) + xlim(-200,600) + 
#   xlab('Residual non-decision time after N200 effect (ms)') +
#   ylab('') +
#   # ggtitle('ERP latency effect posterior distributions') +
#   annotate("segment", x=0, xend=0, y=1, yend=10, colour = 'red',size=3,alpha=.5) +
#   theme(text = element_text(size=20), legend.position="none") +
#   scale_fill_manual(values = tail(cbbPalette,n=7)) +
#   scale_y_discrete(labels=condlabs) +
#   scale_x_continuous(breaks=c(-200,0, 200, 400, 600))
# plot(plot2)
# dev.off()

# png(paste('Model1_N200Lat.png',sep=""),units="in",width=10,height=10,res=300)
# plot3 = ggplot(N200intercepts,aes(x=N200Lat,y=Condition,fill=Condition)) +
#   geom_joy(scale=3, alpha=.67) + xlim(100,300) + 
#   xlab('Hierarchical Mean N200 peak latencies (ms)') +
#   ylab('') +
#   # ggtitle('ERP latency effect posterior distributions') +
#   theme(text = element_text(size=20), legend.position="none") +
#   scale_fill_manual(values = tail(cbbPalette,n=7)) +
#   scale_y_discrete(labels=condlabs) +
#   scale_x_continuous(breaks=c(100, 150, 200, 250, 300))
# plot(plot3)
# dev.off()