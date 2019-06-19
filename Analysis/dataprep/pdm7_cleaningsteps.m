function eeg = pdm7_cleaningsteps(varargin)
%PDM7_CLEANINGSTEPS - Script that follows the steps to clean data
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
%  7/18/17        Michael Nunez                 Original code
%  7/20/17        Michael Nunez  Does not change the current directory
%  8/14/17        Michael Nunez    Use automatic component evaluation

%% Initial
[~,sub,ses]=parsevar(varargin,'sub','getsub','ses','getses');
if strcmp(sub,'getsub')
    subpath = fileparts(pwd);
    [~,sub] = fileparts(subpath);
end
fprintf('Subject: %s \n',sub);

if strcmp(ses,'getses')
    [~,ses] = fileparts(pwd);
end
fprintf('Session: %s \n',ses);

% %% Code
% cd(sprintf('~/data10/michael/pdm/exp7data/pilot/%s/%s',sub,ses));

%Set up for artscreenEEG version 0.5.0 scripts
eeg = pdm7_dataext;

fprintf('%%Loading data into artscreen()...\n');
eeg = artscreen(eeg);

%Save ICA data
fprintf('%%Saving raw (filtered and artifact removed) data...\n');
save(sprintf('%s_%s.mat',sub,ses),'-struct','eeg');

%Independent component analysis on data
fprintf('%%Running ICA on artscreened data...\n');
ica = icasegdata(eeg,'ncomps',100,'nkeep',100,'fftfreq',50,'algorithm','fastica');

%Save ICA data
fprintf('%%Saving ICA data...\n');
save(sprintf('%s_%s_ica.mat',sub,ses),'-struct','ica');

%% Review ICA components
% fprintf('%%Review ICA components...\n');
% fprintf('%%Try ADJUST button...\n');
% ica = icareview(ica);

%Automatic review of ICA components
ica = icaadjust(ica);

%Save ICA data with component evaluations
fprintf('%%Overwriting ICA data with evaluated component information...\n');
save(sprintf('%s_%s_ica.mat',sub,ses),'-struct','ica');

%ICs to Chan data: Keep components which are marked as "unsure"
fprintf('%%Keeping all components marked as "good" or "unsure"...\n');
cleaned = icatochan(ica,1); 

%Remove any probable bad trials observed during "icareview"
fprintf('%%A second pass of artscreen() is recommended...\n');
cleaned = artscreen(cleaned);

%Save cleaned data
fprintf('%%Saving the cleaned data...\n');
save(sprintf('%s_%s_cleaned.mat',sub,ses),'-struct','cleaned');