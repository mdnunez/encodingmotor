%PDM7_FINDALLERPS - Script that saves all ERP data in one location
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
%  10/30/18        Michael Nunez                 Original code
%  11/01/18        Michael Nunez      Export beta desynchronization

%% Initial
subjectdir = '~/data10/michael/pdm/exp7data/subjects';
findsubs = dir(subjectdir);
tracksub = 1;
for p=1:length(findsubs),
	if ~isempty(regexp(findsubs(p).name,'s')) & (findsubs(p).isdir)
		subjects{tracksub} = findsubs(p).name;
		tracksub = tracksub + 1;
	end
end

sessions = {'ses1', 'ses2'};

%% Code
uniqueses = 1;

uniqueses=1;
for m = 1:numel(subjects),
	for n=1:numel(sessions),
		sub = subjects{m};
		ses = sessions{n};
		cd(sprintf('%s/%s/%s',subjectdir,sub,ses));
        
        fileloc = sprintf('%s_%s_sfinal.mat',sub,ses); 
        
        if (exist(fileloc) == 2),
        	fprintf('Loading data from %s ...\n',fileloc);
        	eeg = load(fileloc);
        	if uniqueses==1,
        		N200window = eeg.stimlocked.window;
        		P300window = eeg.stimlocked.P300window;
        		RPwindow = eeg.resplocked.erpwindow;
                BETAwindow = eeg.resplocked.window;
        	end
            N200timecourse(:,:,uniqueses) = eeg.stimlocked.timecourse;
            N200weights(:,uniqueses) = eeg.stimlocked.weights;
            P300timecourse(:,:,uniqueses) = eeg.stimlocked.P300timecourse;
            P300weights(:,uniqueses) = eeg.stimlocked.P300weights;
            RPtimecourse(:,:,uniqueses) = eeg.resplocked.timecourse;
            RPweights(:,:,uniqueses) = eeg.resplocked.weights;

            fprintf('Find the beta desynchronization time course...\n');
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

            BETAtimecourse(:,:,uniqueses) = plotdata;
            uniqueses = uniqueses + 1;
        end
	end
end

save(sprintf('%s/allerps.mat',subjectdir),'N200window','P300window','RPwindow','BETAwindow','BETAtimecourse',...
	'N200timecourse','N200weights','P300timecourse','P300weights','RPtimecourse','RPweights');