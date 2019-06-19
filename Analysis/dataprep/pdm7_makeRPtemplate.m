%PDM7_MAKERPTEMPLATE - Script finds readiness potentials templates from prior experiments
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
%  10/12/18        Michael Nunez     Converted from avgerptopo_resplocked.m


%% Initial
load('~/data10/Stroke2AFC/young/young_erp_resplocked.mat');
thechans = [59 58 51 52 60 66 65 155 164 182 183 184 195 196];

%% Code

unaveragedRP = (allerp.as.left + allerp.as.right)/2;
RPtemplate = mean(mean(unaveragedRP(:,:,thechans),3),1); %average across subjects and channels
RPtemplatewindow = allerp.window;

save('../../Data/RPtemplate.mat','RPtemplate','RPtemplatewindow');
fprintf('Plotting RP template...\n');
figure;
plot(RPtemplatewindow,RPtemplate);
title('RP Waveform Template');