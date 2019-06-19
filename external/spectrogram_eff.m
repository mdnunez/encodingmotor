function [fftdata, freqs, ftime] = spectrogram_eff(data,sr,maxfreq,fftwindowsize,fftstep,windowfn);
%SPECTROGRAM_EFF - AN EFFICIENT VERSION FOR MULTICHANNEL SPECTROGRAM
%INPUT
%data = time x channel x epoch
%sr = sampling rate
%maxfreq = maximum frequency to return
%fftwindowsize = number of samples into fft
%fftstep = how many samples to slide it before next fft
%windowfn = window function to apply before fft. 
%OUTPUTDATA
%fftdata = fourier coefficients limited to 50 Hz
%freqs = frequencies
%times = time points

% Copyright (C) 2017 Michael D. Nunez, <mdnunez1@uci.edu> & Ramesh Srinivasan
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

ntime = size(data,1);
nchan = size(data,2);
nepoch = size(data,3);

nwin = floor((ntime-fftwindowsize)/fftstep);
nfreq = ceil(maxfreq*fftwindowsize/sr)+1;
freqs = [0:nfreq-1]*sr/fftwindowsize;
ftime = [0:nwin-1]*fftstep+fftwindowsize;
fftdata = zeros(nwin,nfreq,nchan,nepoch);

fprintf('Calculating Fourier coefficients for %d windows ...\n',nwin);
fprintf('Calculating window 1');
for j = 1:nwin
	timepts = [j-1]*fftstep + [1:fftwindowsize];
	ddata = ndetrend(data(timepts,:,:),1);
	for k = 1:nepoch
		ddata(:,:,k) = (windowfn*ones(1,nchan)).*squeeze(ddata(:,:,k));
	end;
	fdata = fft(ddata,[],1);
	fftdata(j,:,:,:) = fdata(1:nfreq,:,:);
	if j ~= nwin,
		printiter(j);
	end
end;
fprintf('\n');