%PDM7_BATCH - Script performs data extraction in batch
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
%  11/03/17       Michael Nunez                New subjects

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

fail = zeros(length(subjects),length(sessions));

% %% Code

for m = 1:numel(subjects),
	for n=1:numel(sessions),
		sub = subjects{m};
		ses = sessions{n};
		cd(sprintf('~/data10/michael/pdm/exp7data/subjects/%s/%s',sub,ses));

		if exist(sprintf('%s_%s.mat',sub,ses)) ~= 2, %If file doesn't already exist
			fprintf('Extracting data for %s %s...\n',sub,ses);
			%Set up for artscreenEEG scripts
			try
				eeg = pdm7_dataext;
			catch me
				fprintf('Error in extracting data from %s %s!\n',sub,ses);
				fprintf('Reason: \n');
				fprintf('%s\n',me.message);
				fail(m,n) = 1;
				break;
			end

			%Save raw data
			fprintf('%%Saving raw data...\n');
			save(sprintf('%s_%s.mat',sub,ses),'-struct','eeg');
		end
	end
end