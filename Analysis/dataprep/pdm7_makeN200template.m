%PDM7_makeN200template - Script finds N200 topography template from prior experiments
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
%  7/16/18        Michael Nunez                 Original code
%  10/10/18       Michael Nunez         Save out waveform template as well


%% Initial
svdloc = '~/data10/michael/pdm/%s/%s/%s/erp_svd_%s_%s_v5.mat';
experiments = {'exp4data/subjects', 'exp5data/subjects/training'};
load('eginn128hm.mat');
sessions = {{'s1', 's2'}, {'ses1', 'ses2', 'ses3', 'ses4', 'ses5', 'ses6', 'ses7'}};
subjects = {'s59', 's64', 's68', 's80', 's82', 's93', ...
            's94', 's95', 's96', 's97', 's100', 's101', 's109', 's110'};

%% Code

%% Code
N200template = zeros(1,129);
N200wavetemplate = zeros(1100,1);

for e = 1:numel(experiments),
    for m = 1:numel(subjects),
    	for n=1:numel(sessions{e}),
    		sub = subjects{m};
    		ses = sessions{e}{n};
            filetoload = sprintf(svdloc,experiments{e},sub,ses,sub,ses);
            if exist(filetoload) == 2,
                fprintf('Loading %s ...\n',filetoload);
                load(filetoload);
    	        tempweights = zeros(1,129);
                tempweights(goodchans) = svd_rint_low.v(1,:);
                N200wavetemplate = N200wavetemplate + svd_rint_low.u(:,1);
                % figure(1);
                % hortontopo(tempweights,EGINN128);
                N200template = N200template + tempweights;
                tempweights(goodchans) = svd_rint_med.v(1,:);
                N200wavetemplate = N200wavetemplate + svd_rint_med.u(:,1);
                % figure(2);
                % hortontopo(tempweights,EGINN128);
                N200template = N200template + tempweights;
                tempweights(goodchans) = svd_rint_high.v(1,:);
                N200wavetemplate = N200wavetemplate + svd_rint_high.u(:,1);
                % figure(3);
                % hortontopo(tempweights,EGINN128);
                N200template = N200template + tempweights;
                
            end
    	end
    end
end

save('../../Data/N200template.mat','N200template','N200wavetemplate');
fprintf('Plotting N200 templates...\n');
figure;
hortontopo(N200template,EGINN128);
title('N200 Channel Template');
figure;
plot(-99:1000,N200wavetemplate);
title('N200 Waveform Template');