# pdm5b_scatterplots.R - Creates scatter plots of N200 latencies versus RT percentiles
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
#  10/29/18        Michael Nunez         Converted from pdm5b_scatterplots.R
#  11/02/18        Michael Nunez        Rearrange colors, increase font size

## To do: Generate plotly

## References:
# http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3

## Necessary packages
library(ggplot2)
library(R.matlab)

## Code
fontsize = 22;

# Get linear model equation and R^2_adj as a string
# SOURCE: http://goo.gl/K4yh
lm_eqn <- function(lmresults){
    eq <- substitute(italic(y) == a + b %.% italic(x)*"  "~~italic(R)[adj]^2~"="~r2adj,
         list(a = format(coef(lmresults)[1], digits=0), 
              b = format(coef(lmresults)[2], digits = 2), 
             r2adj = format(summary(lmresults)$adj.r.squared, digits = 2)))
    as.character(as.expression(eq));                 
}

# http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=7
cbbPalette <- c("#D95F02", "#1B9E77", "#7570B3", "#E6AB02", "#66A61E", "#E7298A")
shapePalette <- c(21, 24, 22, 25, 21, 23)

# Read the tables and samples
sesconddata = read.csv(
  '../sesconddata.csv')

colnames(sesconddata)
for (i in seq(1,dim(sesconddata)[1])) {
	if (sesconddata$Condition[i] < 4) {
		expstr = 'Right'
	} else {
		expstr = 'Left'
	}
	if ( (sesconddata$Condition[i] == 1) | (sesconddata$Condition[i] == 4) ) {
		condstr = '.6 s'
	} else if ( (sesconddata$Condition[i] == 2) | (sesconddata$Condition[i] == 5) ) {
		condstr = '.9 s'
	} else {
		condstr = '1.5 s'
	}
	sesconddata$Block[i] = sprintf('%s %s',expstr,condstr)
}
sesconddata$Block = factor(sesconddata$Block, levels = c("Right .6 s", "Right .9 s", "Right 1.5 s","Left .6 s", "Left .9 s", "Left 1.5 s"))

# Remove low accuracy observations (indications that participants did not perform the cognitive task) and boundary N200 latency observations
sesdata = sesconddata[(sesconddata$Accuracy > .6) & (!sesconddata$RemoveN200),]

#Convert to milliseconds
sesdata$N200latencies = sesdata$N200latencies*1000
sesdata$meanRT = sesdata$meanRT*1000


## Plots
lm1 = lm(meanRT ~ N200latencies,sesdata)
png('n200_meanrt.png',units="in",width=10,height=10,res=300)
plot1 = ggplot(sesdata,aes(x=N200latencies, y=meanRT)) +
geom_point(aes(shape=Block, fill=Block),size=3) +
geom_smooth(method="lm", colour='black',size=1.5) +
xlim(150,275) + 
coord_cartesian(ylim=c(400,900)) +
theme(axis.text=element_text(size=fontsize),axis.title=element_text(size=fontsize,face='bold'), legend.text = element_text(size=fontsize),
    legend.title=element_text(size=fontsize,face='bold')) +
labs(x='Trial-averaged N200 peak-latency (ms)',y='Mean RT (ms)',color='Block') +
scale_fill_manual(values = cbbPalette) + 
scale_shape_manual(values = shapePalette) + 
geom_abline(slope=1,intercept=lm1$coefficients[1]-20,colour="#A6761D",size=2, linetype=2) +
annotate("text",x=237,y=900,label="Accuracy < 60% removed",size=8)+
annotate("text",x=166,y=900,label="N = 201",size=8)+
annotate("text",x = 200, y = 875, label = lm_eqn(lm1), parse = TRUE,size=8)
plot(plot1)
dev.off()

# Remove low accuracy observations (indications that participants did not perform the cognitive task) and boundary P300 latency observations
P300data = sesconddata

#Convert to milliseconds
P300data$P300latencies = P300data$P300latencies*1000
P300data$meanRT = P300data$meanRT*1000


## Plots
lm2 = lm(meanRT ~ P300latencies,P300data)
png('p300_meanrt.png',units="in",width=10,height=10,res=300)
plot2 = ggplot(P300data,aes(x=P300latencies, y=meanRT)) +
geom_point(aes(shape=Block, fill=Block),size=3) +
geom_smooth(method="lm", colour='black',size=1.5) +
xlim(275,800) + 
coord_cartesian(ylim=c(400,900)) +
theme(axis.text=element_text(size=fontsize),axis.title=element_text(size=fontsize,face='bold'), legend.text = element_text(size=fontsize),
    legend.title=element_text(size=fontsize,face='bold')) +
labs(x='Trial-averaged P300 peak-latency (ms)',y='Mean RT (ms)',color='Block') +
scale_fill_manual(values = cbbPalette) + 
scale_shape_manual(values = shapePalette) + 
geom_abline(slope=1,intercept=0,colour="#A6761D",size=2, linetype=2) +
# annotate("text",x=650,y=900,label="Accuracy < 60% removed",size=8)+
annotate("text",x=350,y=900,label="N = 324",size=8)+
annotate("text",x = 500, y = 875, label = lm_eqn(lm2), parse = TRUE,size=8)
plot(plot2)
dev.off()
