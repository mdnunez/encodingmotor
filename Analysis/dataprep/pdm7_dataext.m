function pdmdata = pdm7_dataext(varargin)
%PDM7_DATAEXT - Extracts xdf data and cleans photocells from PDM experiment7
%
%Usage: pdmdata = pdm7_dataext(nblocks, MaxTrials, MinTrials, xdfdirectory);
%
%
%Inputs:
%  nblock: Number of blocks
%  MaxTrials: Maximum number of trials in each block
%  MinTrials: Minimum number of trials in each block
%  xdfdirectory: Directory with only 'nblocks' + 2 .xdf files corresponding to 
%'nblocks' blocks
%
%Outputs:
% pdmdata: A structure to be used in Cort's functions
%
%To save output (example):
%save s00.mat -struct pdmdata
%
%To load output (example):
%pdmdata = load('s00.mat')
%
%
% Copyright (C) 2017 Michael D. Nunez, <mdnunez1@uci.edu>, Cort Horton, Bill Winter
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
%  7/11/17        Michael Nunez             Adapted from pdm7_dataext
%  7/20/17        Michael Nunez            For use with one photocell

%% Initial
if nargin < 1
    nblocks = 6;
    MaxTrials = 80;
    MinTrials = 80;
    xdfdir = pwd;
elseif nargin == 1
    nblocks = varargin{1};
    MaxTrials = 80;
    MinTrials = 80;
    xdfdir = pwd;
elseif nargin == 2
    nblocks = varargin{1};
    MaxTrials = varargin{2};
    MinTrials = 80;
    xdfdir = pwd;
elseif nargin == 3
    nblocks = varargin{1};
    MaxTrials = varargin{2};
    MinTrials = varargin{3};
    xdfdir = pwd;
elseif nargin == 3
    nblocks = varargin{1};
    MaxTrials = varargin{2};
    MinTrials = varargin{3};
    xdfdir = varargin{4};
else
    error('Too many inputs.');
end

xdf = dir([xdfdir,'/*.xdf']);

if numel(xdf) == nblocks+2
    eyes = 1;
elseif numel(xdf) == nblocks
    eyes = 0;
else
    error('Incorrect number of xdf files found! Please move files.')
end

fprintf('File load order:\n');
for bnum = 1:nblocks
    fprintf('%s ',xdf(bnum).name);
end
fprintf('\n');

behavdir = dir([xdfdir,'/behav/*.mat']);
behavfile = behavdir.name;
fprintf('Behavior file name: %s\n',behavfile);

pdmdata.expinfo = load([xdfdir,'/behav/',behavfile]);
flickerhz = pdmdata.expinfo.flickerhz;
noisehz = pdmdata.expinfo.noisehz;

%% Code
badchans = [125:136];

try
    for bnum = 1:nblocks
        disp(['Cleaning block ',num2str(bnum)]);

        streams = load_xdf([xdfdir,'/',xdf(bnum).name]);

        clear raw srate marker
        %Find the stream containing openvibe data
        for s=1:length(streams),
            if strcmp(streams{s}.info.name,'openvibeSignal')
                raw = streams{s}.time_series;
                eegstamps = streams{s}.time_stamps;
                srate = str2num(streams{s}.info.nominal_srate);
                break;
            end
        end

        for s=1:length(streams),
            if strcmp(streams{s}.info.name,'openvibeMarkers')
                marker = streams{s};               
                break;
            end
        end

        %Filtering step ([1 100] bandpass and detrending-baselining) using Cort's filtereeg function
        [tempdata, ~] = rmchannel(raw',badchans); %Remove photocell channels and average reference the data
        block = pdm7_filtereeg(tempdata,srate);
        [block, ~] = rmchannel(block,[]); %Average reference the remaining data
        
        nbounds = round(srate/(noisehz/2) +[-10 10]);
        sbounds = round(srate/flickerhz +[-10 10]); 

        lengthsamps = floor(srate*1.5);

        spikesignal = getspike(raw(134,:),sbounds(1)); %Photocell of 30 Hz Signal
        psignal = zeros(size(spikesignal));
        psignal(spikesignal > sbounds(1) & spikesignal < sbounds(2)) = ...
            spikesignal(spikesignal > sbounds(1) & spikesignal < sbounds(2));
        
        nonzero = (psignal > 0);
        respstart = strfind(nonzero,[zeros(1,lengthsamps) 1]) + lengthsamps;
        fprintf('%i response intervals were found in block %i \n',length(respstart),bnum);

        %Find where light detector remains on to find feedback interval
        lighton = (abs(raw(134,:)) > std(abs(raw(134,:))));
        feedstart = strfind(lighton, [zeros(1,ceil(.2*srate)), ones(1,ceil(.4*srate))]) + ceil(.2*srate);
        fprintf('%i feedback intervals were found in block %i \n',length(feedstart),bnum);

        ntrials = length(respstart);
        
        missing(bnum) = 0;
        if ntrials < MinTrials
            missing(bnum) = MinTrials - ntrials;
            respstart = [nan(1,missing(bnum)) respstart];
        end
        if ntrials > MaxTrials
            error(sprintf('Number of trials found is larger the number of trials desired in block %i \n',bnum));
        end

        %Note that in marker.time_series: 
        %33028 denotes button push
        %33024 denotes button release
        %33030 denotes right button light detector on (attached to LED indicating button push on Cedrus box)
        %33026 denotes right button light detector off (attached to LED indicating button push on Cedrus box)
        %33029 denotes left button light detector on (attached to LED indicating button push on Cedrus box)
        %33025 denotes left button light detector off (attached to LED indicating button push on Cedrus box)
        leftbutton = zeros(1,size(raw,2));
        rightbutton = zeros(1,size(raw,2));

        buttonpush_times = marker.time_stamps(marker.time_series == 33028);
        for p=1:(length(marker.time_series)-1),
            if marker.time_series(p) == 33028, % if the marker denotes a button push...
                if marker.time_series(p+1) == 33029, %if the next marker after a button push denotes a left button light detector on then...
                    [~,whatsamp] = min(abs(eegstamps - marker.time_stamps(p))); %find closest sample in EEG time
                    leftbutton(whatsamp) = 1; %record as true in leftbutton time series
                elseif marker.time_series(p+1) == 33030 %if the next marker after a button push denotes a right button light detector on then...
                    [~,whatsamp] = min(abs(eegstamps - marker.time_stamps(p))); %find closest sample in EEG time
                    rightbutton(whatsamp) = 1; %record as true in rightbutton time series
                end
            end
        end
        
        %Put data into equal trial windows (4000 sample points to cover all possible window lengths of
        %.25+1+2+.5+.25 seconds) locked on first SIGNAL photocell and approximately on the last
        %photocell (exact photocell placement differences by one or two samples
        %between trials)
        %Note that there may be some EEG overlap
        resplen = ceil(2.750*srate);
        sampsbefore = ceil(1.25*srate);
        data{bnum} = zeros(resplen+sampsbefore,129,ntrials);
        photo{bnum} = zeros(resplen+sampsbefore,ntrials);
        timestamps{bnum} = zeros(resplen+sampsbefore,1,ntrials);
        % button{bnum} = zeros(resplen+sampsbefore,2,ntrials); %Extract behavioral data
        correct{bnum} = nan(1,ntrials);
        rt{bnum} = nan(1,ntrials);
        feedon{bnum} = nan(1,ntrials);
        for i=(missing(bnum)+1):(ntrials+missing(bnum))
            data{bnum}(:,:,i) = block((respstart(i)-sampsbefore):(respstart(i)+resplen-1),1:129); %Segment 250 ms before trial start
            photo{bnum}(:,i) = [spikesignal((respstart(i)-sampsbefore):(respstart(i)+resplen-1))'];
            timestamps{bnum}(:,:,i) = eegstamps((respstart(i)-sampsbefore):(respstart(i)+resplen-1));

            button{bnum}(:,:,i) = [leftbutton((respstart(i)-sampsbefore):(respstart(i)+resplen-1))' , ...
                rightbutton((respstart(i)-sampsbefore):(respstart(i)+resplen-1))'];

            %%Get RT and accuracy
            %Correct, incorrect or answered inappropriately?
            thistrial = i + MinTrials*(bnum-1) + missing(bnum);
            if pdmdata.expinfo.highlow(thistrial) == 2 && any(button{bnum}(:,2,i)) && ~any(button{bnum}(:,1,i))
                correct{bnum}(i) = 1;
                tempsamps = find(button{bnum}(:,2,i));
            elseif pdmdata.expinfo.highlow(thistrial) == 1 && any(button{bnum}(:,1,i)) && ~any(button{bnum}(:,2,i))
                correct{bnum}(i) = 1;
                tempsamps = find(button{bnum}(:,1,i));
            elseif pdmdata.expinfo.highlow(thistrial) == 2 && any(button{bnum}(:,1,i)) && ~any(button{bnum}(:,2,i))
                correct{bnum}(i) = 0;
                tempsamps = find(button{bnum}(:,1,i));
            elseif pdmdata.expinfo.highlow(thistrial) == 1 && any(button{bnum}(:,2,i)) && ~any(button{bnum}(:,1,i))
                correct{bnum}(i) = 0;
                tempsamps = find(button{bnum}(:,2,i));
            else
                tempsamps = NaN;
            end
            %%Check this
            %rt{bnum}(i) = tempsamps(1) - respstart(i) + 1;
            rt{bnum}(i) = tempsamps(1) - sampsbefore;
        end
        feedon{bnum} = sampsbefore + (feedstart-respstart);
        feedon{bnum} = feedon{bnum}((missing(bnum)+1):(ntrials+missing(bnum)));

    end
catch me
    rethrow(me)
    if bnum>1
        disp('Saving all other data...');
        alldata = [];
        allphoto = [];
        allbutton = [];
        allcorrect = [];
        allrt = [];
        allfeedon = [];
        for i=1:(bnum-1)
            alldata = cat(3,alldata,data{i});
            allphoto = cat(3,allphoto,photo{i});
            allbutton = cat(3,allbutton,button{i});
            allcorrect = [allcorrect,correct{i}];
            allrt = [allrt,rt{i}];
            allfeedon = [allfeedon,feedon{i}];
        end
        pdmdata.data = alldata(:,1:129,:);
        pdmdata.photo = allphoto;
        pdmdata.button = allbutton;
        pdmdata.feedon = allfeedon;
        pdmdata.expinfo.correct = allcorrect;
        pdmdata.expinfo.rt = allrt;
        pdmdata.sr = srate;
        pdmdata = addhm(pdmdata,'eginn128');
        pdmdata.error = me;
        pdmdata.missing = missing;
    end
    return
end

alldata = [];
allphoto = [];
allbutton = [];
allcorrect = [];
allrt = [];
allfeedon = [];
for i=1:nblocks
    alldata = cat(3,alldata,data{i});
    allphoto = cat(3,allphoto,photo{i});
    allbutton = cat(3,allbutton,button{i});
    allcorrect = [allcorrect,correct{i}];
    allrt = [allrt,rt{i}];
    allfeedon = [allfeedon,feedon{i}];
end

pdmdata.data = alldata(:,1:129,:);
pdmdata.photo = allphoto;
pdmdata.button = allbutton;
pdmdata.feedon = allfeedon;
pdmdata.expinfo.correct = allcorrect;
pdmdata.expinfo.rt = allrt;
pdmdata.sr = srate;
pdmdata = addhm(pdmdata,'eginn128');
pdmdata.missing = missing;


%% ----------------------------------------------------------------------%%
    function [dataout, badchans] = rmchannel(data,markbad,zeroassign)
        % RMCHANNEL - Removes specified channels from data (either zeros them out
        % or assigns them NANs), average references data
        %
        % Usage: [dataout badchans] = rmchannel(data,badchans)
        %
        % Required Input:
        %    data - sample*chan or sample*chan*trial data, already bad channels are 
        %           NAN or 0
        %    markbad - list of channel indicies for "bad channel" assignment
        % 
        % Optional Input:
        %    zeroassign - Zeros out bad channels (and does not include them in
        %                 average reference calculation) instead of assigning 
        %                 channels NANs (default: 1)
        %
        % Output:
        %    dataout - average (good channels) referenced sample*chan(*trial) data
        %    badchans - new bad channel list (now either NaN or 0's in dataout)

        %% Record of Revisions
        %   Date           Programmers               Description of change
        %   ====        =================            =====================
        %  8/12/14        Michael Nunez                  Original code
        %  8/13/14        Michael Nunez              Finds bad channels near zero
        %  4/06/15        Michael Nunez            Default now zero assigned

        %% Initial
        if nargin < 1
            error('Not enough arguments! See help file')
        elseif nargin > 3
            error('Too many arguments! See help file')
        elseif nargin == 1
            markbad = [];
            zeroassign = 1;
        elseif nargin == 2
            zeroassign = 1;
        end

        %% Code

        %Marks channels that are NaN or near zero as bad
        badindex = all(all(isnan(data),1),3) | all(all((abs(data) < 10^(-9)),1),3);

        %New bad channel list is the union of found list and user defined list
        badchans = union(find(badindex),markbad);
        goodchans = setdiff(1:size(data,2),badchans);

        %Average reference the data
        datamean = nanmean(data(:,goodchans,:),2);
        datamean = repmat(datamean,[1 length(goodchans)]);
        dataout = data;
        dataout(:,goodchans,:) = data(:,goodchans,:) - datamean;

        %Assign bad channels as NaN or 0
        if zeroassign
            dataout(:,badchans,:) = zeros(size(data,1),length(badchans),size(data,3));
        else
            dataout(:,badchans,:) = nan(size(data,1),length(badchans),size(data,3));
        end

    end %end of rmchannel()


    function data=pdm7_filtereeg(data,sr)
        % function fdata=pdm7_filtereeg(data,sr)
        %
        %Filter specifications for experiment type exp5
        %
        %   data - time by chan (by trial)
        %   sr - sampling rate

        %% Record of Revisions
        %   Date           Programmers               Description of change
        %   ====        =================            =====================
        %  5/03/13   Cort Horton, mod: Michael N.   Modified From Cort's for PDM
        %  1/12/15     Michael Nunez                 Low pass stopband: 100
        %  7/10/17     Michael Nunez                  Addition of filtfilthd()
        %
        %% Code

        % Detrend the data
        data=ndetrend(data);

        % Highpass the data
        Fstop = 0.25;        % Stopband Frequency
        Fpass = 1;           % Passband Frequency
        Astop = 10;          % Stopband Attenuation (dB)
        Apass = 1;           % Passband Ripple (dB)
        match = 'passband';  % Band to match exactly

        % Construct an FDESIGN object and call its BUTTER method.
        h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, sr);
        Hd = design(h, 'butter', 'MatchExactly', match);

        % Carry out the filtering
        data=filtfilthd(Hd,data);

        % Notch at 60 Hz
        Fpass1 = 59;          % First Passband Frequency
        Fstop1 = 59.9;        % First Stopband Frequency
        Fstop2 = 60.1;        % Second Stopband Frequency
        Fpass2 = 61;          % Second Passband Frequency
        Apass1 = 1;           % First Passband Ripple (dB)
        Astop  = 10;          % Stopband Attenuation (dB)
        Apass2 = 1;           % Second Passband Ripple (dB)
        match  = 'stopband';  % Band to match exactly

        % Construct an FDESIGN object and call its BUTTER method.
        h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
                              Apass2, sr);
        Hd = design(h, 'butter', 'MatchExactly', match);

        % Carry out the filtering
        data=filtfilthd(Hd,data);

        % Lowpass the data
        Fpass = 40;          % Passband Frequency
        Fstop = 50;          % Stopband Frequency
        Apass = 1;           % Passband Ripple (dB)
        Astop = 10;          % Stopband Attenuation (dB)
        match = 'passband';  % Band to match exactly

        % Construct an FDESIGN object and call its BUTTER method.
        h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, sr);
        Hd = design(h, 'butter', 'MatchExactly', match);

        % Carry out the filtering
        data=filtfilthd(Hd,data);

    end %end of pdm7_filtereeg()


    function [data w] = getspike(data,gap,npeak)
        % cleanstim = getspike(stim,gap,npeak)
        %
        % GETSPIKE processes the data in the stimulus array into a spike train.
        % Spikes less than 'gap' away from the prior spike are dropped.  The result
        % has a power of 1 per period.  GETSPIKE returns the leading edge for
        % positively-oriented spike trains and the lagging edge for negatively
        % oriented spike trains.
        %
        %   stim: channel x sample array of stimulus data
        %    gap: refractory gap
        %  npeak: maximum number of cycles to return
        %
        % Written by Bill Winter, March 2006
        [nChan nFrame] = size(data);
        if nargin < 2 || isempty(gap), gap = 0; end
        if nargin < 3 || isempty(npeak), npeak = size(data,2); end
        w = zeros(nChan,2);
        for k = 1:nChan
            p = getgap(getzeros(ndetrend(data(k,:),2,0)),gap);
            p(npeak+1:end) = [];
            data(k,:) = 0;
            data(k,p) = [diff(p) nFrame-p(end)+1];
            w(k,:) = [p(1) p(end)-1];
        end
    end %end of getspike()

end %end of main function