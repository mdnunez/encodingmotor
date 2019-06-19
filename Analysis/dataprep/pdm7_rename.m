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
%  11/02/17        Michael Nunez                 Original code

%% Initial
subjects = {'s94',  's112', 's116', 's143', 's144', 's145', 's146', 's147', 's148', 's149', ...
			's150', 's151', 's152', 's153', 's154', 's155', 's156', 's157', 's158', 's159'};

sessions = {'ses1', 'ses2'};

% %% Code

for m = 1:numel(subjects),
	for n=1:numel(sessions),
		sub = subjects{m};
		ses = sessions{n};
		cd(sprintf('~/data10/michael/pdm/exp7data/subjects/%s/%s',sub,ses));

		%Set up for artscreenEEG version 0.5.0 scripts
		if exist(sprintf('%s_%s_cleaned.mat',sub,ses)) == 2
			fprintf('Renaming files for %s %s\n',sub,ses);
			system(sprintf('mv %s_%s.mat %s_%s_mcleaned.mat',sub,ses,sub,ses));
			system(sprintf('mv %s_%s_ica.mat %s_%s_mica.mat',sub,ses,sub,ses));
			system(sprintf('mv %s_%s_cleaned.mat %s_%s_mfinal',sub,ses,sub,ses));
		end
	end
end