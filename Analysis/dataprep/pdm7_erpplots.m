function eeg = pdm7_erpplots(varargin)
%PDM7_ERPPLOTS - Script that plots ERP data
%
% Copyright (C) 2018 Michael D. Nunez, <mdnunez1@uci.edu>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Record of Revisions
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  8/14/17        Michael Nunez                 Original code
%  8/15/17        Michael Nunez            Add beta desynch plots
%  03/21/17       Michael Nunez        Plot top 4 SVD components
%  07/16/18       Michael Nunez       Load _sfinal.mat data
%  07/18/18       Michael Nunez          Improve N200 window, improve beta estimation
%  10/10/18		  Michael Nunez           Add eegstruct input
%  10/11/18       Michael Nunez             Add alpha plot
%  10/18/18       Michael Nunez           Plot response locked ERPs
%  10/19/18       Michael Nunez          Remove response-locked theta
%  10/23/18       Michael Nunez        Collapse across similar conditions to improve ERP estimates
%  10/29/18       Michael Nunez        Collapse across similar conditions to improve ERP estimates per hand
%  10/30/18       Michael Nunez           Export eeg structure

%To do:
%1) Plot response time curves
%2) Calculate stimulus-locked power

%% Initial
[~,sub,ses,eeg]=parsevar(varargin,'sub','getsub','ses','getses','eegstruct',[]);
if strcmp(sub,'getsub')
    subpath = fileparts(pwd);
    [~,sub] = fileparts(subpath);
end
fprintf('Subject: %s \n',sub);

if strcmp(ses,'getses')
    [~,ses] = fileparts(pwd);
end
fprintf('Session: %s \n',ses);


%Load cleaned data
if isempty(eeg),
	fprintf('%%Loading the cleaned data...\n');
	eeg = load(sprintf('%s_%s_sfinal.mat',sub,ses));
end

ntrials = size(eeg.data,3);
nchans = size(eeg.data,2);

disp('Assuming all non-zero and non-NaN data are good...');
thevars=squeeze(var(eeg.data));
artifact=thevars==0 | isnan(thevars);

if isfield(eeg,'artifact')
     disp('Rejecting data from .artifact field matrix...');
     artifact = eeg.artifact | artifact;
end

goodtrials=setdiff(1:ntrials,find(sum(artifact)==nchans));
goodchans=setdiff(1:nchans,find(sum(artifact,2)==ntrials));

%Plot ERPs
fprintf('Plotting stimulus-locked ERPs...\n');
f1 = figure('units','normalized','outerposition',[0 0 1 1]);
cortplotx(eeg.stimlocked.window, squeeze(eeg.stimlocked.timecourse(:,:)), 'Linewidth',2);
line1 = line([0 0],ylim);
set(line1, 'Linewidth', 2, 'Color', 'k', 'LineStyle', '--');
line2 = line([125 125],ylim);
set(line2, 'Linewidth', 2, 'Color', 'r', 'LineStyle', '--');
line3 = line([275 275],ylim);
set(line3, 'Linewidth', 2, 'Color', 'r', 'LineStyle', '--');
legend({'Right Hand: .6 seconds', 'Right Hand: .9 seconds', 'Right Hand: 1.5 seconds', ...
    'Left Hand: .6 seconds', 'Left Hand: .9 seconds', 'Left Hand: 1.5 seconds'}, 'Fontsize', 20,'Location','SouthEast');
xlabel('Time after Gabor stimulus (ms)','Fontsize',20);

f3 = figure('units','normalized','outerposition',[0 0 1 1]);
figure(f3);
hortontopo(eeg.stimlocked.weights(:),eeg.hm,'weights',1);

f8 = figure('units','normalized','outerposition',[0 0 1 1]);
cortplotx(eeg.stimlocked.P300window, squeeze(eeg.stimlocked.P300timecourse(:,:)), 'Linewidth',2);
% Pz = 62; %electrode Pz
% cortplotx(eeg.stimlocked.P300window, squeeze(eeg.stimlocked.P300erp(:,Pz,:)), 'Linewidth',2);
line1 = line([0 0],ylim);
set(line1, 'Linewidth', 2, 'Color', 'k', 'LineStyle', '--');
line2 = line([800 800],ylim);
set(line2, 'Linewidth', 2, 'Color', 'r', 'LineStyle', '--');
line3 = line([275 275],ylim);
set(line3, 'Linewidth', 2, 'Color', 'r', 'LineStyle', '--');
legend({'Right Hand: .6 seconds', 'Right Hand: .9 seconds', 'Right Hand: 1.5 seconds', ...
    'Left Hand: .6 seconds', 'Left Hand: .9 seconds', 'Left Hand: 1.5 seconds'}, 'Fontsize', 20,'Location','SouthEast');
set(gca,'YTick',[],'Fontsize',20);
xlabel('Time after Gabor stimulus (ms)','Fontsize',20);

f9 = figure('units','normalized','outerposition',[0 0 1 1]);
hortontopo(eeg.stimlocked.P300weights(:),eeg.hm,'weights',1);



fprintf('Plotting response-locked ERPs...\n');
f7 = figure('units','normalized','outerposition',[0 0 1 1]);
cortplotx(eeg.resplocked.erpwindow, squeeze(eeg.resplocked.timecourse(:,:)), 'Linewidth',2);
line4 = line([0 0],ylim);
set(line4, 'Linewidth', 2, 'Color', 'k', 'LineStyle', '--');
legend({'Right Hand: .6 seconds', 'Right Hand: .9 seconds', 'Right Hand: 1.5 seconds', ...
    'Left Hand: .6 seconds', 'Left Hand: .9 seconds', 'Left Hand: 1.5 seconds'}, 'Fontsize', 20,'Location','SouthWest');
set(gca,'YTick',[],'Fontsize',20);
xlabel('Time after Gabor stimulus (ms)','Fontsize',20);

f10 = figure('units','normalized','outerposition',[0 0 1 1]);
for h=1:2, %per condition (i.e. block)

	figure(f10);
	subplot(1,2,h);
	hortontopo(eeg.resplocked.weights(:,h),eeg.hm,'weights',1);

end

fprintf('Plotting response-locked beta-desynchronization...\n');
f2 = figure('units','normalized','outerposition',[0 0 1 1]);
murhythm = find((eeg.resplocked.freqs==16 | eeg.resplocked.freqs==20)); %Window 14-22 Hz
C3 = 36; %left motor cortex electrode (expected right hand desynchronization location)
C4 = 104; %right motor cortex electrode (expected left hand desynchronization location)
smoothingwindowsize = 25;
filterorder = 1;
betainput = [squeeze(mean(eeg.resplocked.powdata(:,murhythm,C3,1:3),2)) squeeze(mean(eeg.resplocked.powdata(:,murhythm,C4,4:6),2))];
smoothbeta = sgolayfilt(betainput,filterorder,smoothingwindowsize); %Frequency range 14-22 Hz
% smoothbeta = betainput;
baselinesamples = [205:256];
bpower = mean(smoothbeta(baselinesamples,:));
bpowermat = ones(size(smoothbeta,1),1)*bpower;
plotdata = (smoothbeta - bpowermat)./bpowermat;
cortplotx(eeg.resplocked.window,plotdata, 'Linewidth',2);
line1 = line([0 0],ylim);
set(line1, 'Linewidth', 2, 'Color', 'k', 'LineStyle', '--');
line2 = line([-600 -600],ylim);
set(line2, 'Linewidth', 2, 'Color', 'r', 'LineStyle', '--');
line3 = line([-200 -200],ylim);
set(line3, 'Linewidth', 2, 'Color', 'r', 'LineStyle', '--');
legend({'Right Hand: .6 seconds', 'Right Hand: .9 seconds', 'Right Hand: 1.5 seconds', ...
	'Left Hand: .6 seconds', 'Left Hand: .9 seconds', 'Left Hand: 1.5 seconds'}, 'Fontsize', 20,'Location','SouthWest');
set(gca,'YTick',[],'Fontsize',20);
xlabel('Time before response (ms)','Fontsize',20);

fprintf('Plotting response-locked alpha-desynchronization...\n');
f5 = figure('units','normalized','outerposition',[0 0 1 1]);
alpharhythm = find((eeg.resplocked.freqs==8 | eeg.resplocked.freqs==12)); %Window 6-14 Hz
Pz = 62; %electrode Pz
Oz = 75; %electrode Oz
smoothingwindowsize = 25;
filterorder = 1;
alphainput = squeeze(mean(eeg.resplocked.powdata(:,alpharhythm,Pz,:),2));
smoothalpha = sgolayfilt(alphainput,filterorder,smoothingwindowsize); %Frequency range 6-14 Hz
% smoothbeta = betainput;
baselinesamples = [205:256];
apower = mean(smoothalpha(baselinesamples,:));
apowermat = ones(size(smoothalpha,1),1)*apower;
plotalpha = (smoothalpha - apowermat)./apowermat;
cortplotx(eeg.resplocked.window,plotalpha, 'Linewidth',2);
line4 = line([0 0],ylim);
set(line1, 'Linewidth', 2, 'Color', 'k', 'LineStyle', '--');
% line5 = line([-600 -600],ylim);
% set(line5, 'Linewidth', 2, 'Color', 'r', 'LineStyle', '--');
% line6 = line([-200 -200],ylim);
% set(line6, 'Linewidth', 2, 'Color', 'r', 'LineStyle', '--');
legend({'Right Hand: .6 seconds', 'Right Hand: .9 seconds', 'Right Hand: 1.5 seconds', ...
	'Left Hand: .6 seconds', 'Left Hand: .9 seconds', 'Left Hand: 1.5 seconds'}, 'Fontsize', 20,'Location','SouthWest');
set(gca,'YTick',[],'Fontsize',20);
xlabel('Time before response (ms)','Fontsize',20);

% fprintf('Plotting response-locked theta-synchronization...\n');
% f6 = figure('units','normalized','outerposition',[0 0 1 1]);
% thetarhythm = find((eeg.resplocked.freqs==4)); %Window 2-6 Hz
% Fz = 11; %electrode Fz
% smoothingwindowsize = 25;
% filterorder = 1;
% thetainput = squeeze(mean(eeg.resplocked.powdata(:,thetarhythm,Pz,:),2));
% smooththeta = sgolayfilt(thetainput,filterorder,smoothingwindowsize); %Frequency range 6-14 Hz
% % smoothbeta = betainput;
% baselinesamples = [205:256];
% tpower = mean(smooththeta(baselinesamples,:));
% tpowermat = ones(size(smooththeta,1),1)*tpower;
% plottheta = (smooththeta - tpowermat)./apowermat;
% cortplotx(eeg.resplocked.window,plottheta, 'Linewidth',2);
% line4 = line([0 0],ylim);
% set(line1, 'Linewidth', 2, 'Color', 'k', 'LineStyle', '--');
% % line5 = line([-600 -600],ylim);
% % set(line5, 'Linewidth', 2, 'Color', 'r', 'LineStyle', '--');
% % line6 = line([-200 -200],ylim);
% % set(line6, 'Linewidth', 2, 'Color', 'r', 'LineStyle', '--');
% legend({'Right Hand: .6 seconds', 'Right Hand: .9 seconds', 'Right Hand: 1.5 seconds', ...
% 	'Left Hand: .6 seconds', 'Left Hand: .9 seconds', 'Left Hand: 1.5 seconds'}, 'Fontsize', 20,'Location','SouthWest');
% set(gca,'YTick',[],'Fontsize',20);
% xlabel('Time before response (ms)','Fontsize',20);


fprintf('Plotting stimulus-locked single-trial ERPs...\n');
fprintf('Plotting response time curves...\n');
f4 = figure('units','normalized','outerposition',[0 0 1 1]);
f11 = figure('units','normalized','outerposition',[0 0 1 1]);
for c=1:6, %per condition (i.e. block)

	condtrials = find(eeg.expinfo.condition == c);

	figure(f4);
	subplot(2,3,c);
	alltrials = intersect(goodtrials,condtrials);
	lines = cortplotx(eeg.stimlocked.window, squeeze(eeg.stimlocked.sterp(:,alltrials)), 'Linewidth',2);
	line1 = line([0 0],ylim);
	set(line1, 'Linewidth', 2, 'Color', 'k', 'LineStyle', '--');
	line2 = line([125 125],ylim);
	set(line2, 'Linewidth', 2, 'Color', 'r', 'LineStyle', '--');
	line3 = line([275 275],ylim);
	set(line3, 'Linewidth', 2, 'Color', 'r', 'LineStyle', '--');
	set(gca,'YTick',[],'Fontsize',20);

	figure(f11);
	subplot(2,3,c);
	rtacc = (eeg.expinfo.rt(alltrials)/eeg.sr).*(2*eeg.expinfo.correct(alltrials) - 1);
	[ks,xi] = ksdensity(rtacc(~isnan(rtacc)));
	lines = cortplotx(xi,ks,'Linewidth',4);
	set(gca,'YTick',[],'Fontsize',20);
	xlabel('Reaction time (secs)','Fontsize',20);
	ylabel(sprintf('%.2f %% Accuracy',nanmean(rtacc>0)*100),'Fontsize',20);

end

