%PDM7_BATCH5 - Script performs pdm7_analysissteps() and pdm7_tidydata() in batch
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
%  8/18/17        Michael Nunez                 Converted from pdm7_batch4.m


%% Initial
findsubs = dir('/data10/michael/pdm/exp7data/subjects/');
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
		cd(sprintf('/data10/michael/pdm/exp7data/subjects/%s/%s',sub,ses));
        
        
        if (exist(sprintf('%s_%s_sfinal.mat',sub,ses)) == 2),
            pdm7_analysissteps;
            pdm7_tidydata(uniqueses);
            uniqueses = uniqueses + 1;
        end
	end
end
