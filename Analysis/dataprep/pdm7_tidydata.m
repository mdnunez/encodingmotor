function pdm7_tidydata(uniqueses,varargin)
%PDM7_TIDYDATA - Script that exports table data to csv file
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
%  8/15/17        Michael Nunez                 Original code
%  8/17/17        Michael Nunez      Write tables with removed data
%  11/6/17        Michael Nunez           Update file name
%  07/18/18       Michael Nunez          Improve N200 window, improve beta estimation
%                                        Do not remove any reaction times
%  10/25/18       Michael Nunez        Improve .csv creation
%  10/29/18       Michael Nunez        Find P300 peaks and RP minimums

% To do:
% 1) Calculate beta desynchronization at the found electrode locations for the Readiness Potential

%% Initial
[~,sub,ses]=parsevar(varargin,'sub','getsub','ses','getses');
if strcmp(sub,'getsub')
    subpath = fileparts(pwd);
    [~,sub] = fileparts(subpath);
end
fprintf('Subject: %s \n',sub);

if strcmp(ses,'getses')
    [originaldir,ses] = fileparts(pwd);
end
fprintf('Session: %s \n',ses);

%Load cleaned data
fprintf('Loading the cleaned data...\n');
eeg = load(sprintf('%s_%s_sfinal.mat',sub,ses));

%% Save data

encodingmotorloc = fileparts(which(mfilename));
fprintf('''encodingmotor'' repository location found at %s \n', encodingmotorloc);
cd(encodingmotorloc);
cd('..');
cd('..');
cd('Data');
datafolder = pwd;

trialcsvloc = sprintf('%s/trialdata.csv',datafolder);
trialcsvNANloc = sprintf('%s/trialdataNAN.csv',datafolder);
if exist(trialcsvloc) ~= 2,
	fprintf('Creating file %s...\n',trialcsvloc);
	fprintf('Creating file %s...\n',trialcsvNANloc);
	fileID1 = fopen(trialcsvloc,'w');
	fprintf(fileID1, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...)
		'RT',...
		'Accuracy',...
		'Condition',...
		'SubjectID',...
		'Subjects_session',...
		'Unique_Session',...
		'Noisetimes',...
		'Intertrial',...
		'Stimtime',...
		'Spatial_frequencies',...
		'Random_rotations',...
		'Speed_cutoff',...
		'Given_feedback',...
		'N200amplitudes',...
		'N200latencies',...
		'RemoveN200',...
		'RemoveRT',...
		'Training',...
		'Artifact');
	fclose(fileID1);
	fileID2 = fopen(trialcsvNANloc,'w');
	fprintf(fileID2, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...)
		'RT',...
		'Accuracy',...
		'Condition',...
		'SubjectID',...
		'Subjects_session',...
		'Unique_Session',...
		'Noisetimes',...
		'Intertrial',...
		'Stimtime',...
		'Spatial_frequencies',...
		'Random_rotations',...
		'Speed_cutoff',...
		'Given_feedback',...
		'N200amplitudes',...
		'N200latencies',...
		'RemoveN200',...
		'RemoveRT',...
		'Training',...
		'Artifact');
	fclose(fileID2);
end

sescondcsvloc = sprintf('%s/sesconddata.csv',datafolder);
sescondcsvNANloc = sprintf('%s/sesconddataNAN.csv',datafolder);
if exist(sescondcsvloc) ~= 2
	fprintf('Creating file %s...\n',sescondcsvloc);
	fprintf('Creating file %s...\n',sescondcsvNANloc);
	fileID3 = fopen(sescondcsvloc,'w');
	fprintf(fileID3,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
		'meanRT', ...
		'RT10prctile', ...
		'Accuracy', ...
		'Condition', ...
		'SubjectID', ...
		'Subjects_session', ...
		'Unique_Session', ...
		'N200amplitudes', ...
		'N200latencies', ...
		'RemoveN200', ...
		'minBETAamplitudes', ...
		'minBETAlatencies', ...
		'RemoveminBETA', ...
		'P300amplitudes', ...
		'P300latencies', ...
		'RemoveP300', ...
		'minRPamplitudes', ...
		'minRPlatencies', ...
		'RemoveminRP');
	fclose(fileID3);
	fileID4 = fopen(sescondcsvNANloc,'w');
	fprintf(fileID4,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
		'meanRT', ...
		'RT10prctile', ...
		'Accuracy', ...
		'Condition', ...
		'SubjectID', ...
		'Subjects_session', ...
		'Unique_Session', ...
		'N200amplitudes', ...
		'N200latencies', ...
		'RemoveN200', ...
		'minBETAamplitudes', ...
		'minBETAlatencies', ...
		'RemoveminBETA', ...
		'P300amplitudes', ...
		'P300latencies', ...
		'RemoveP300', ...
		'minRPamplitudes', ...
		'minRPlatencies', ...
		'RemoveminRP');
	fclose(fileID4);
end


%Find single-trial N200 estimates
thesesamps = find((eeg.stimlocked.window >= 125) & (eeg.stimlocked.window <= 275));
[stmins, stargmins] = min(eeg.stimlocked.sterp(thesesamps,:),[],1);
stn200 = stmins;
stn200lat = (eeg.stimlocked.window(thesesamps(stargmins)))/1000; % in seconds
removestn200 =  ( (stargmins == 1) | (stargmins == length(thesesamps)) ); %Remove boundary minimums

%Find trial-averaged N200 estimates
[mins, argmins] = min(eeg.stimlocked.timecourse(thesesamps,:),[],1);
n200 = mins;
n200lat = (eeg.stimlocked.window(thesesamps(argmins)))/1000; % in seconds
removen200 =  ( (argmins == 1) | (argmins == length(thesesamps)) ); %Remove boundary minimums

%Find trial-averaged P300 estimates
P300samps = find((eeg.stimlocked.P300window >= 275) & (eeg.stimlocked.P300window <= 800));
[maxs, argmaxs] = max(eeg.stimlocked.P300timecourse(P300samps,:),[],1);
p300 = maxs;
p300lat = (eeg.stimlocked.P300window(P300samps(argmaxs)))/1000; % in seconds
removep300 =  ( (argmaxs == 1) | (argmaxs == length(P300samps)) ); %Remove boundary maximums

%Find trial-averaged RP minimums
RPsamps = find((eeg.resplocked.erpwindow <= 0) & (eeg.resplocked.erpwindow >= -500));
[mins, argmins] = min(eeg.resplocked.timecourse(RPsamps,:),[],1);
minRP = mins;
minRPlat = abs(eeg.resplocked.erpwindow(RPsamps(argmins)))/1000; % in seconds
removeRP =  ( (argmins == 1) | (argmins == length(RPsamps)) ); %Remove boundary minimums

%Find trial-averaged beta desynchronization minimum estimates
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
baselinedbeta = (smoothbeta - bpowermat)./bpowermat;


respsamps = find((eeg.resplocked.window >= -600) & (eeg.resplocked.window <= -200));
[BETAmins, BETAargmins] = min(baselinedbeta(respsamps,:),[],1);
minBETA = BETAmins;
minBETAlat = -(eeg.resplocked.window(respsamps(BETAargmins)))/1000; % in seconds before response time
removeBETA =  ( (BETAargmins == 1) | (BETAargmins == length(respsamps)) ); %Remove boundary minimums

%Calculate and remove RT
rt = (eeg.expinfo.rt/1024)';
removert = (rt < .35) | isnan(rt);

meanrt = nan(6,1);
accuracy = nan(6,1);
rt10prctile = nan(6,1);
for c=1:6
	condtrials = find(eeg.expinfo.condition == c);
	goodtrials = intersect(find(~removert),condtrials);

	meanrt(c) = mean(rt(goodtrials));
	rt10prctile(c) = prctile(rt(goodtrials),10);
	accuracy(c) = mean(eeg.expinfo.correct(goodtrials));
end

%Find EEG artifact trials
ntrials = size(eeg.data,3);
nchans = size(eeg.data,2);

disp('Assuming all non-zero and non-NaN data are good...');
thevars=squeeze(var(eeg.data));
artifact=thevars==0 | isnan(thevars);

if isfield(eeg,'artifact')
     disp('Rejecting data from .artifact field matrix...');
     artifact = eeg.artifact | artifact;
end

artifacttrials = (sum(artifact)==nchans);

%Find trials to remove
removeall = artifacttrials' | removert;

trialM = [rt eeg.expinfo.correct', eeg.expinfo.condition', ...
	ones(ntrials,1)*str2num(sub(2:end)), ones(ntrials,1)*eeg.expinfo.sesnum, ones(ntrials,1)*uniqueses, ...
	eeg.expinfo.noisetimes', eeg.expinfo.intertrial', eeg.expinfo.stimtime', eeg.expinfo.spfs', ...
	eeg.expinfo.randrots', eeg.expinfo.speed_cutoff', eeg.expinfo.given_feedback', stn200', stn200lat', removestn200', removert, ...
	eeg.traintest.training', artifacttrials'];

dlmwrite(trialcsvloc, trialM,'-append','delimiter',',','precision',5);

sescondM = [meanrt, rt10prctile, accuracy (1:6)' ones(6,1)*str2num(sub(2:end)) ones(6,1)*eeg.expinfo.sesnum ones(6,1)*uniqueses ...
	n200' n200lat' removen200' minBETA' minBETAlat' removeBETA' p300' p300lat' removep300' minRP' minRPlat' removeRP'];

dlmwrite(sescondcsvloc, sescondM,'-append','delimiter',',','precision',5);

%Remove trials
removeall = artifacttrials' | removert | ~eeg.traintest.training' | removestn200';
rt(removeall) = nan;
eeg.expinfo.correct(removeall) = nan;
stn200(removeall) = nan;
stn200lat(removeall) = nan;

trialM = [rt eeg.expinfo.correct', eeg.expinfo.condition', ...
	ones(ntrials,1)*str2num(sub(2:end)), ones(ntrials,1)*eeg.expinfo.sesnum, ones(ntrials,1)*uniqueses, ...
	eeg.expinfo.noisetimes', eeg.expinfo.intertrial', eeg.expinfo.stimtime', eeg.expinfo.spfs', ...
	eeg.expinfo.randrots', eeg.expinfo.speed_cutoff', eeg.expinfo.given_feedback', stn200', stn200lat', removestn200', removert, ...
	eeg.traintest.training', artifacttrials'];

dlmwrite(trialcsvNANloc, trialM,'-append','delimiter',',','precision',5);

%Remove trial-averaged N200
n200(removen200) = nan;
n200lat(removen200) = nan;

%Remove trial-averaged P300
p300(removep300) = nan;
p300lat(removep300) = nan;

%Remove beta desynchronization minimums
minBETA(removeBETA) = nan;
minBETAlat(removeBETA) = nan;

sescondM = [meanrt, rt10prctile, accuracy (1:6)' ones(6,1)*str2num(sub(2:end)) ones(6,1)*eeg.expinfo.sesnum ones(6,1)*uniqueses ...
	n200' n200lat' removen200' minBETA' minBETAlat' removeBETA' p300' p300lat' removep300' minRP' minRPlat' removeRP'];

dlmwrite(sescondcsvNANloc, sescondM,'-append','delimiter',',','precision',5);

cd(originaldir);