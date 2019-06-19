%PDM7_MAKEP300TEMPLATE - Script finds readiness potentials templates from prior experiments
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
%  10/18/18        Michael Nunez     Converted from pdm7_makeRPtemplate.m

%To do:
%1) Extract weights for each electrode like in pdm7_makeN200template()
%   Convert to topographic images then interpolate to 128 electrode positions


%% Initial
load('~/data10/Stroke2AFC/young/young_erp_stimlocked_RTthirds.mat');
parchans = [100 110 109 108 128 127 99 119]; %parietal chans

%% Code

unaveragedP300 = (allerp.as.fastest + allerp.as.medfast + allerp.as.slowest)/3;
P300template = mean(mean(unaveragedP300(:,:,parchans),3),1); %average across subjects and channels
P300templatewindow = allerp.window;

save('../../Data/P300template.mat','P300template','P300templatewindow');
fprintf('Plotting P300 template...\n');
figure;
plot(P300templatewindow,P300template);
title('P300 Waveform Template');