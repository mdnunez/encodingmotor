%PDM7_BATCH2 - Script performs initial artscreen() step in batch
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
%  8/21/17        Michael Nunez                "
%  11/02/17       Michael Nunez           Append with _mcleaned.mat
%  11/03/17       Michael Nunez               New subjects

%To do:
%1) Use automatic EEG trial rejection scheme:
%Windows with amps exceeding 70 muV, with amplitude changes exceeding 100 muV within 100 ms, or with activity lower than 0.5 muV.

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

for m = 1:numel(subjects),
	for n=1:numel(sessions),
		sub = subjects{m};
		ses = sessions{n};
		cd(sprintf('~/data10/michael/pdm/exp7data/subjects/%s/%s',sub,ses));

        %If _mcleaned DOESN'T already exist and extracted data DOES exist
		if (exist(sprintf('%s_%s_mcleaned.mat',sub,ses)) ~= 2) & (exist(sprintf('%s_%s.mat',sub,ses)) == 2),
			%Load raw data
			fprintf('Loading %s %s raw data...\n',sub,ses);
			eeg = load(sprintf('%s_%s.mat',sub,ses));

			fprintf('Loading data into artscreen()...\n');
			eeg = artscreen(eeg);

			%Save artscreened data
			fprintf('Saving raw (filtered and artifact removed) data...\n');
			save(sprintf('%s_%s_mcleaned.mat',sub,ses),'-struct','eeg');
		end
	end
end