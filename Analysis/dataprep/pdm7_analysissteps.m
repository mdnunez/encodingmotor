function eeg = pdm7_analysissteps(varargin)
%PDM7_ANALYSISSTEPS - Script that follows the steps to analyze the EEG data
%
% Copyright (C) 2018 Michael D. Nunez, <mdnunez1@uci.edu> & Ramesh Srinivasan
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
%  8/14/17        Michael Nunez                 Original code
%  8/15/17        Michael Nunez   Adding power time course calculation
%  8/21/17        Michael Nunez            Remove missing RT trials
%  11/6/17        Michael Nunez     Update file names
%                  Major error fixed: Fixed the filtering sample rate!!!
%  03/21/17       Michael Nunez        Save out more SVD components
%  07/16/18       Michael Nunez       Load _sfinal.mat data
%  07/18/18       Michael Nunez          Save out original ERP
%                                       Use template to find matching component
%                                       Change filtering parameters
%  10/10/18       Michael Nunez     Find the N200 in only the first 10 components                   
%  10/12/18       Michael Nunez        Remove very slow trials (due to lack of samples)
%                        Match the waveform template instead of the scalp topography
%              Fix V(:, comp) mistake on line to compute eeg.stimlocked.timecourse(:,c)
%  10/18/18       Michael Nunez            Fixes, filter the RP canidate data
%  10/23/18       Michael Nunez        Collapse across similar conditions to improve ERP estimates
%  10/24/18       Michael Nunez         Calculate N200 per condition
%  10/26/18       Michael Nunez        Remove negative RTs
%  10/29/18       "                Collapse across hand to improve RP estimates


%Note: 1) N200 often is not captured by the first SVD component in this experiment
% 2) The point of deviation between hands may be informative for the time at which decision-making starts
%    E.g. the point at which the evidence accumulation tracks start to differ
%    Use subtraction to find these points

%To do:
%1) Extract weights for each electrode like in pdm7_makeN200template()
%   Convert to topographic images then interpolate to 128 electrode positions
%2) Calculate Lateralized Readiness Potential
%3) Calculate stimulus-locked power
%4) Calculate alpha power in N200 topography location, relate to the diffusion coefficient


%% Initial
[~,sub,ses,eeg]=parsevar(varargin,'sub','getsub','ses','getses','eegstruct',[]);
if strcmp(sub,'getsub')
    subpath = fileparts(pwd);
    [~,sub] = fileparts(subpath);
end
fprintf('Subject: %s \n',sub);

if strcmp(ses,'getses')
    [~,ses] = fileparts(pwd);
end
fprintf('Session: %s \n',ses);


%Load cleaned data
if isempty(eeg),
	fprintf('%%Loading the cleaned data...\n');
	eeg = load(sprintf('%s_%s_sfinal.mat',sub,ses));
end

%Filter data
% Lowpass the data
Fpass = 7;          % Passband Frequency
Fstop = 10;          % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 10;          % Stopband Attenuation (dB)
match = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, eeg.sr);
Hd = design(h, 'butter', 'MatchExactly', match);

% 2nd lowpass filter
% Lowpass the data
Fpass2 = 4;          % Passband Frequency
Fstop2 = 8;          % Stopband Frequency
Apass2 = 1;           % Passband Ripple (dB)
Astop2 = 10;          % Stopband Attenuation (dB)

% Construct an FDESIGN object and call its BUTTER method.
h2  = fdesign.lowpass(Fpass2, Fstop2, Apass2, Astop2, 1000);
Hd2 = design(h2, 'butter', 'MatchExactly', match);

%% Stimulus event-related potential calculations (N200 amplitude and latency)

erpwind = 1178:2304; %100 samples before stimulus onset to 1000 ms after stimulus onset, stimulus onset is sample 1281, sample rate is 1024
baseline = 1:103; %103 samples before the stimulus onset is the baseline (~100 ms)
N200window = 231:385; %between 125 and 275 milliseconds
P300window = 412:1436; %between 0 and 1000 milliseconds

ntrials = size(eeg.data,3);
nchans = size(eeg.data,2);

eeg.stimlocked.weights = nan(129,1);
eeg.stimlocked.timecourse = nan(length(erpwind),3);
eeg.stimlocked.sterp = nan(length(erpwind),ntrials);
eeg.stimlocked.erp = nan(length(erpwind),nchans);
eeg.stimlocked.P300weights = nan(129,1);
eeg.stimlocked.P300timecourse = nan(length(erpwind),3);
eeg.stimlocked.P300sterp = nan(length(erpwind),ntrials);
eeg.stimlocked.P300erp = nan(length(erpwind),nchans);

eeg.traintest.training = zeros(1,ntrials);
eeg.resplocked.powdata = nan(256,10,nchans,6);
eeg.stimlocked.N200window = N200window;
eeg.stimlocked.baseline = baseline;
eeg.stimlocked.erpwind = erpwind;

[goodchans, goodtrials] = goodbad(eeg);
goodtrials = intersect(goodtrials,find((eeg.expinfo.rt < 2560) & (eeg.expinfo.rt > 0))); %Only submit trials with reaction times less than 2560 samps (2500 ms)

%Load the N200 topography template
fprintf('Loading N200 waveform template...\n');
tempout = load('N200template');
N200template = tempout.N200template;
N200wavetemplate = resample(tempout.N200wavetemplate,1024,1000); %Convert from 1000 samps/sec to 1024 samps/sec

%Load the P300 topography template
fprintf('Loading P300 waveform template...\n');
tempout = load('P300template');
P300wavetemplate = resample(tempout.P300template,1024,1000); %Convert from 1000 samps/sec to 1024 samps/sec

%Load the RP topography template
fprintf('Loading RP waveform template...\n');
rpout = load('RPtemplate');
RPwindow = resample(rpout.RPtemplatewindow,1024,1000);
subwindow = 207:1230; %1000 ms before to response
RPwavetemplate = resample(rpout.RPtemplate,1024,1000); %Convert from 1000 samps/sec to 1024 samps/sec


fprintf('Calculating stimulus-locked ERPs for N200 waveforms...\n');
rng(str2num(sub(2:end))+1000*str2num(ses(4))); %Set seed with subject number and session number
trainindx = randperm(length(eeg.expinfo.condition));
traintrials = trainindx(1:round(.75*length(eeg.expinfo.condition)));
eeg.traintest.training(traintrials) = 1;

erptrials = intersect(goodtrials,traintrials);

% Carry out the filtering
fprintf('Low pass filtering stimulus-locked ERP data at 10 Hz...\n');
filtered10=filtfilthd(Hd,eeg.data(erpwind,:,:));

%Find the N200 waveform over all conditions
dataforerp = filtered10(:,:,erptrials);
dataforerp(:, setdiff(1:nchans,goodchans), :) = 0; %zero out bad channels
dataforerp = eegBaseline(dataforerp,baseline); %baseline the data
erp = mean(dataforerp,3); %Calculate the ERP
eeg.stimlocked.erp(:,:) = erp;

%Use singular value decomposition
[U,S,V] = svd(erp(N200window,:),'econ');

%Find the component that best matches the N200 template (and flip signs if necessary), equally weighting the waveform and scalp topography
templatematch = corr(N200wavetemplate(N200window),[U(:,1:10) -U(:,1:10)]) + corr(N200template',[V(:,1:10) -V(:,1:10)]);
% templatematch = N200template*[V(:,1:10) -V(:,1:10)];
comp = mod(argmax(templatematch),10); %mod() ensures a sign flip to match the template
fprintf('SVD Component %d of the ERP best matched the N200 template! \n',comp);
if argmax(templatematch) > 10,
	% Flip weights in the components
	V(:, comp) = -V(:, comp);
	U(:, comp) = -U(:, comp);
end

eeg.stimlocked.weights(:) = V(:, comp);
eeg.stimlocked.window = (-102:1024)/1.024;

% Calculate single-trial ERP estimates for all trials
alltrials = goodtrials;
dataforsterp = filtered10(:,:,alltrials);
allstdata=segtocat(dataforsterp);

% Place into single trials
cmpsig=allstdata*V(:, comp);
eeg.stimlocked.sterp(:,alltrials) =squeeze(cattoseg(cmpsig,length(erpwind)));

fprintf('Calculating stimulus-locked ERPs per cutoff condition for P300 waveforms...\n');


% Carry out the filtering
fprintf('Low pass filtering stimulus-locked ERP data at 4 Hz...\n');
filtered4=filtfilthd(Hd2,eeg.data(erpwind,:,:));

%Find the P300 waveform over all conditions
dataforerp = filtered4(:,:,erptrials);
dataforerp(:, setdiff(1:nchans,goodchans), :) = 0; %zero out bad channels
dataforerp = eegBaseline(dataforerp,baseline); %baseline the data
P300erp = mean(dataforerp,3); %Calculate the ERP
eeg.stimlocked.P300erp(:,:) = P300erp;

%Use singular value decomposition
[Up300,Sp300,Vp300] = svd(P300erp(103:1127,:),'econ');

%Find the component that best matches the P300 template (and flip signs if necessary)
P300templatematch = corr(P300wavetemplate(P300window)',[Up300(:,1:10) -Up300(:,1:10)]);
% templatematch = N200template*[V(:,1:10) -V(:,1:10)];
compP300 = mod(argmax(P300templatematch),10); %mod() ensures a sign flip to match the template
fprintf('SVD Component %d best matched the P300 template! \n',compP300);
if argmax(P300templatematch) > 10,
	% Flip weights in the components
	Vp300(:, compP300) = -Vp300(:, compP300);
	Up300(:, compP300) = -Up300(:, compP300);
end
eeg.stimlocked.P300weights(:) = Vp300(:, compP300);
eeg.stimlocked.P300window = (-102:1024)/1.024;

% Calculate single-trial ERP estimates for all trials
alltrials = goodtrials;
dataforstP300 = filtered4(:,:,alltrials);
allstP300data=segtocat(dataforstP300);

% Place into single trials
cmpsig2=allstP300data*Vp300(:, compP300);
eeg.stimlocked.P300sterp(:,alltrials) =squeeze(cattoseg(cmpsig2,length(erpwind)));

hands = {'right','left'};
for h=1:2,
	handtrials = intersect(erptrials,find(( eeg.expinfo.condition > 3*(h-1) ) & ( eeg.expinfo.condition <= 3*h )));
	%Find the component that best matches the RP template and flip signs if necessary
	[respdata, lockindex, badtrials] = epochsubset(eeg.data(:,:,handtrials), eeg.expinfo.rt(handtrials)+1280);

	respdata = respdata((lockindex-1279):(lockindex+256),:,:); %Lock each window to 1280 samples (1250 ms) before and 256 samples (250 ms) after response


	 % Carry out the filtering
	fprintf('Low pass filtering response-locked ERP data at 4 Hz...\n');
	filteredresp = filtfilthd(Hd2,respdata);

	%Find the RP waveform
	resperpdata = filteredresp;
	resperpdata(:, setdiff(1:nchans,goodchans), :) = 0; %zero out bad channels
	rpbaseline = (size(resperpdata,1)-255):size(resperpdata,1); %at response to 250 ms after response
	resperpdata = eegBaseline(resperpdata,rpbaseline);
	resperp = mean(resperpdata,3); %Calculate the response-locked ERP
	if h==1,
		eeg.resplocked.erp = NaN(1536,129,2);
		eeg.resplocked.erpwindow = linspace(-1249,250,1536); %-1250 to 250 ms in sample space
		eeg.resplocked.weights = NaN(129,2);
		eeg.resplocked.timecourse = NaN(1536,6);
	end
	eeg.resplocked.erp(:,:,h) = resperp;

	%Use singular value decomposition
	[Uh,Sh,Vh] = svd(resperp((size(resperp,1)-1279):(size(resperp,1)-256),:),'econ');

	%Find the component that best matches the RP template (and flip signs if necessary)
	RPtemplatematch = corr(RPwavetemplate(subwindow)',[Uh(:,1:10) -Uh(:,1:10)]);

	handcomp = mod(argmax(RPtemplatematch),10); %mod() ensures a sign flip to match the template
	fprintf('SVD Component %d of the %s hand ERP best matched the RP template! \n',handcomp, hands{h});
	if argmax(RPtemplatematch) > 10,
		% Flip weights in the components
		Vh(:, handcomp) = -Vh(:, handcomp);
		Uh(:, handcomp) = -Uh(:, handcomp);
	end

	handsvd{h}.U = Uh;
	handsvd{h}.S = Sh;
	handsvd{h}.V = Vh;
	handsvd{h}.comp = handcomp;

	eeg.resplocked.weights(:,h) = Vh(:, handcomp);
end


for c=1:6, %per condition (i.e. block)	
	condtrials = find(eeg.expinfo.condition == c);
	
	thesetrain = intersect(condtrials,traintrials);
	erptrials = intersect(goodtrials,thesetrain);

	%Find the N200 waveform for every cutoff condition
	dataforerp = filtered10(:,:,erptrials);
	dataforerp(:, setdiff(1:nchans,goodchans), :) = 0; %zero out bad channels
	dataforerp = eegBaseline(dataforerp,baseline); %baseline the data
	erp = mean(dataforerp,3); %Calculate the ERP

	%Find the N200 waveform for every cutoff condition
	eeg.stimlocked.timecourse(:,c) = erp*((1/S(comp,comp))*V(:,comp)')'; % Derived from identity: M = USV => U = M*inv(SV) => U = MV'*inv(S) => Ustar = MstarV'*inv(S)

	
	%Find the P300 waveform for every cutoff condition
	dataforerp = filtered4(:,:,erptrials);
	dataforerp(:, setdiff(1:nchans,goodchans), :) = 0; %zero out bad channels
	dataforerp = eegBaseline(dataforerp,baseline); %baseline the data
	P300erp = mean(dataforerp,3); %Calculate the ERP
	
	eeg.stimlocked.P300timecourse(:,c) = P300erp*((1/Sp300(compP300,compP300))*Vp300(:,compP300)')'; % Derived from identity: M = USV => U = M*inv(SV) => U = MV'*inv(S) => Ustar = MstarV'*inv(S)
	

	% Calculate beta desynchronization time course
	[respdata, lockindex, badtrials] = epochsubset(eeg.data(:,:,erptrials), eeg.expinfo.rt(erptrials)+1280);
    
    respdata = respdata((lockindex-1279):(lockindex+256),:,:); %Lock each window to 1280 samples (1250 ms) before and 256 samples (250 ms) after response
	fftwindowsize = ceil(.25*eeg.sr);
    fftstep = round(.005*eeg.sr);
    maxfreq = 35;
    windowfn = gausswin(fftwindowsize,3);
    fprintf('Calculating windowed fft...\n');
    tic;
    [fftdata, freqs, ftime] = spectrogram_eff(respdata,eeg.sr,maxfreq,fftwindowsize,fftstep,windowfn);
    mintook = toc/60;
    fprintf('Windowed fft took %.3f minutes!\n',mintook);

    eeg.resplocked.powdata(:,:,:,c) = var(fftdata,[],4);
    % eeg.resplocked.powdata(:,:,:,c) = sgolayfilt(var(fftdata,[],4),2,61,[],1);
    
    eeg.resplocked.lockindex(c) = lockindex;
    eeg.resplocked.freqs = freqs;
    eeg.resplocked.window = (ftime/1.024)-1280; %time before response in ms

    
    % Carry out the filtering
	fprintf('Low pass filtering response-locked ERP data at 4 Hz...\n');
    filteredresp = filtfilthd(Hd2,respdata);

    %Find the RP waveform
    resperpdata = filteredresp;
    resperpdata(:, setdiff(1:nchans,goodchans), :) = 0; %zero out bad channels
    rpbaseline = (size(resperpdata,1)-255):size(resperpdata,1); %at response to 250 ms after response
    resperpdata = eegBaseline(resperpdata,rpbaseline);
    resperp = mean(resperpdata,3); %Calculate the response-locked ERP

    %Which hand?
    if c<4,
    	h=1;
    else
    	h=2;
    end

	eeg.resplocked.timecourse(:,c) = resperp*((1/handsvd{h}.S(handsvd{h}.comp,handsvd{h}.comp))*handsvd{h}.V(:,handsvd{h}.comp)')'; % Derived from identity: M = USV => U = M*inv(SV) => U = MV'*inv(S) => Ustar = MstarV'*inv(S)

end

%Save cleaned data
fprintf('Overwriting the cleaned data...\n');
save(sprintf('%s_%s_sfinal.mat',sub,ses),'-struct','eeg');

%% ----------------------------------------------------------------------%%
	function [newdata lockindex badtrials] = epochsubset(data,newindex,varargin)
	%EPOCHSUBSET - Reepochs each epoch of EEG data by timelocking each epoch
	%              "i" to "newindex(i)" with maximum window size available
	%
	%Useage:  [newdata lockindex badtrials] = epochsubset(data,newindex,lockindex);
	%    OR:  [newdata lockindex badtrials] = epochsubset(data,newindex);
	%
	%Inputs:  data - sample*channel*trial EEG data
	%         newindex - Vector of length "trial"
	%
	%Optional inputs
	%          lockindex - Sample in which newdata is timelocked
	%                      Default: nanmin(newindex)       
	%
	%Outputs:  newdata - Re-timelocked EEG data
	%          lockindex - Sample in which newdata is timelocked
	%          badtrials - Index of trials where newindex contained nans

		%% Code

		if ~isvector(newindex)
		    error('newindex must be a vector!');
		end

		if size(data,3) ~= length(newindex)
		    error('size(data,3) must be equal to length(newindex)');
		end

		if ~isempty(varargin) & isnumeric(varargin{3})
			lockindex = round(varargin{3});
		else
			lockindex = round(nanmin(newindex));
		end
		fprintf('Using lockindex %d \n',lockindex);

		windsize = (size(data,1) - nanmax(newindex)) + lockindex;

		newdata = zeros([windsize, size(data,2), size(data,3)]);

		for t=1:length(newindex)
		    if ~isnan(newindex(t))
		        begin_index = newindex(t)-lockindex + 1;
		        end_index = windsize + begin_index -1;
		        newdata(:,:,t) = data(begin_index:end_index,:,t);
		    end
		end

		badtrials = find(isnan(newindex));
	end %end of epochsubset

	function [fftdata, freqs, ftime] = spectrogram_eff(data,sr,maxfreq,fftwindowsize,fftstep,windowfn);
	%AN EFFICIENT VERSION FOR MULTICHANNEL SPECTROGRAM
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

		ntime = size(data,1);
		nchan = size(data,2);
		nepoch = size(data,3);

		nwin = floor((ntime-fftwindowsize)/fftstep);
		nfreq = ceil(maxfreq*fftwindowsize/sr)+1;
		freqs = [0:nfreq-1]*sr/fftwindowsize;
		ftime = [0:nwin-1]*fftstep+fftwindowsize;
		fftdata = zeros(nwin,nfreq,nchan,nepoch);
		for j = 1:nwin
		timepts = [j-1]*fftstep + [1:fftwindowsize];
		ddata = ndetrend(data(timepts,:,:),1);
		for k = 1:nepoch
		ddata(:,:,k) = (windowfn*ones(1,nchan)).*squeeze(ddata(:,:,k));
		end;
		fdata = fft(ddata,[],1);
		fftdata(j,:,:,:) = fdata(1:nfreq,:,:);
		end;
	end %end of spectogram_eff


end %end of main function