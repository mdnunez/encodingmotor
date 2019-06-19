%PDM7_BATCH3 - Script performs all other steps in batch
%
% Copyright (C) 2017 Michael D. Nunez, <mdnunez1@uci.edu>
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
%  8/18/17        Michael Nunez                 Original code
%  8/19/17        Michael Nunez                  Split subjects
%  8/21/17        Michael Nunez     Existence of _ica.mat, split subs
%  11/02/17       Michael Nunez        Different file naming scheme
%  11/06/17       Michael Nunez                New subjects

%% Initial
findsubs = dir('~/data10/michael/pdm/exp7data/subjects/');
subjects = {};
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

for m = 1:numel(subjects),
	for n=1:numel(sessions),
		sub = subjects{m};
		ses = sessions{n};
		cd(sprintf('~/data10/michael/pdm/exp7data/subjects/%s/%s',sub,ses));

		if (exist(sprintf('%s_%s_mcleaned.mat',sub,ses)) == 2),

			if exist(sprintf('%s_%s_mica.mat',sub,ses)) ~= 2
				%Load raw data
				fprintf('Loading %s %s filtered data...\n', sub, ses);
				eeg = load(sprintf('%s_%s_mcleaned.mat',sub,ses));

				%Independent component analysis on data
				fprintf('Running ICA on artscreened data...\n');
				ica = icasegdata(eeg,'ncomps',100,'nkeep',100,'fftfreq',50,'algorithm','fastica');

				%Save ICA data
				fprintf('Saving ICA data...\n');
				save(sprintf('%s_%s_mica.mat',sub,ses),'-struct','ica');

				%Automatic review of ICA components
				ica = icaadjust(ica);

				%Save ICA data with component evaluations
				fprintf('Overwriting ICA data with evaluated component information...\n');
				save(sprintf('%s_%s_mica.mat',sub,ses),'-struct','ica');
			else
				fprintf('Loading %s %s ICA data with evaluated component information...\n',sub,ses);
				ica = load(sprintf('%s_%s_mica.mat',sub,ses));
			end

			%ICs to Chan data: Keep components which are marked as "unsure"
			fprintf('Keeping all components marked as "good" or "unsure"...\n');
			cleaned = icatochan(ica,1); 

			%Remove negative reaction times
			cleaned.expinfo.rt(cleaned.expinfo.rt < 0) = nan;

			%Save cleaned data
			fprintf('Saving the %s %s ICA cleaned data...\n',sub,ses);
			save(sprintf('%s_%s_mfinal.mat',sub,ses),'-struct','cleaned');

			%Run analysis steps
			pdm7_analysissteps;

			%Place data into table
			pdm7_tidydata(uniqueses);
			uniqueses = uniqueses + 1;
		end
	end
end