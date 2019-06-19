%PDM7_ERPLOC - Script that generates average ERP response topographies
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
%  11/02/18        Michael Nunez             Adapted from pdm5b_n200loc

% Load data

fprintf('Loading previously generated ERP scalp topographies...\n');
erps = load('~/data10/michael/pdm/exp7data/subjects/allerps.mat');

% %% Plot SVD-localized N200 peak-latency potentials

% fprintf('Plotting the mean ERP...\n');
% figure;
% plotx(-99:1000,squeeze(mean(scalperp,1)));

fprintf('Loading head models...\n');
load ~/data10/finalheadmodel/HNL128_ELECTRODE.mat %This contains bad channels for Laplacian
load ~/data10/michael/headmodels/s05.mat
load eginn128hm.mat

Electrode.badchan.eeg = union([19],Electrode.badchan.eeg); %Add dead channel 19 to bad channel list
goodchan = setdiff(1:size(Electrode.Coordinate,1),Electrode.badchan.eeg);

Electrode = makegeolap(Electrode,Scalp); %Get Laplacian matrix
Electrode.directspline3D = splineinterp(.1,Electrode.Coordinate(goodchan,:),Scalp.Vertex); %Constant is based on Bill's topo.m
Electrode.spline3D = splineinterp(.1,Electrode.Coordinate,Scalp.Vertex);

viewpos = {{-60,8},{60,8},{0,0}};

%% Average plots

%Plot the mean response
fprintf('Plotting the N200 weights on a 3D scalp...\n');
scalpPotential = mean(erps.N200weights(goodchan,:),2)'*Electrode.directspline3D';

%Zero out periphery values
u1=find(Scalp.Vertex(:,3)<-40);
u3=find(Scalp.Vertex(:,2)>39 & Scalp.Vertex(:,3)<14);
scalpPotential(u1)=0;
scalpPotential(u3)=0;


f1=figure('units','normalized','outerposition',[0 0 1 1]);
h1 = drawmesh(Scalp);
setmesh(h1,'interp',-scalpPotential);
set(gca,'CLim', [-max(abs(scalpPotential)) max(abs(scalpPotential))]);
newmap = (bone + jet)/2;
colormap(newmap);
view(gca,viewpos{1}{:});

f2=figure('units','normalized','outerposition',[0 0 1 1]);
h2 = drawmesh(Scalp);
setmesh(h2,'interp',-scalpPotential);
set(gca,'CLim', [-max(abs(scalpPotential)) max(abs(scalpPotential))]);
newmap = (bone + jet)/2;
colormap(newmap);
view(gca,viewpos{2}{:});


fprintf('Saving the scalp potential...\n');
save('ScalpN200Potential.mat','scalpPotential');


fprintf('Plotting the P300 weights on a 3D scalp...\n');
scalpP300Potential = mean(erps.P300weights(goodchan,:),2)'*Electrode.directspline3D';

%Zero out periphery values
scalpP300Potential(u1)=0;
scalpP300Potential(u3)=0;


f1=figure('units','normalized','outerposition',[0 0 1 1]);
h1 = drawmesh(Scalp);
setmesh(h1,'interp',scalpP300Potential);
set(gca,'CLim', [-max(abs(scalpP300Potential)) max(abs(scalpP300Potential))]);
newmap = (bone + jet)/2;
colormap(newmap);
view(gca,viewpos{1}{:});

f2=figure('units','normalized','outerposition',[0 0 1 1]);
h2 = drawmesh(Scalp);
setmesh(h2,'interp',scalpP300Potential);
set(gca,'CLim', [-max(abs(scalpP300Potential)) max(abs(scalpP300Potential))]);
newmap = (bone + jet)/2;
colormap(newmap);
view(gca,viewpos{2}{:});

fprintf('Saving the scalp potential...\n');
save('ScalpP300Potential.mat','scalpP300Potential');

fprintf('Plotting the RP right hand weights on a 3D scalp...\n');
scalpRPrightPotential = mean(squeeze(erps.RPweights(goodchan,1,:)),2)'*Electrode.directspline3D';

%Zero out periphery values
scalpRPrightPotential(u1)=0;
scalpRPrightPotential(u3)=0;


f1=figure('units','normalized','outerposition',[0 0 1 1]);
h1 = drawmesh(Scalp);
setmesh(h1,'interp',-scalpRPrightPotential);
set(gca,'CLim', [-max(abs(scalpRPrightPotential)) max(abs(scalpRPrightPotential))]);
newmap = (bone + jet)/2;
colormap(newmap);
view(gca,viewpos{1}{:});

f2=figure('units','normalized','outerposition',[0 0 1 1]);
h2 = drawmesh(Scalp);
setmesh(h2,'interp',-scalpRPrightPotential);
set(gca,'CLim', [-max(abs(scalpRPrightPotential)) max(abs(scalpRPrightPotential))]);
newmap = (bone + jet)/2;
colormap(newmap);
view(gca,viewpos{2}{:});

fprintf('Saving the scalp potential...\n');
save('ScalpRPrightPotential.mat','scalpRPrightPotential');

fprintf('Plotting the RP left hand weights on a 3D scalp...\n');
scalpRPleftPotential = mean(squeeze(erps.RPweights(goodchan,2,:)),2)'*Electrode.directspline3D';

%Zero out periphery values
scalpRPleftPotential(u1)=0;
scalpRPleftPotential(u3)=0;


f1=figure('units','normalized','outerposition',[0 0 1 1]);
h1 = drawmesh(Scalp);
setmesh(h1,'interp',-scalpRPleftPotential);
set(gca,'CLim', [-max(abs(scalpRPleftPotential)) max(abs(scalpRPleftPotential))]);
newmap = (bone + jet)/2;
colormap(newmap);
view(gca,viewpos{1}{:});

f2=figure('units','normalized','outerposition',[0 0 1 1]);
h2 = drawmesh(Scalp);
setmesh(h2,'interp',-scalpRPleftPotential);
set(gca,'CLim', [-max(abs(scalpRPleftPotential)) max(abs(scalpRPleftPotential))]);
newmap = (bone + jet)/2;
colormap(newmap);
view(gca,viewpos{2}{:});

fprintf('Saving the scalp potential...\n');
save('ScalpRPleftPotential.mat','scalpRPleftPotential');


