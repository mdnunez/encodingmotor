function pdmexp7(trainflag)
%% Record of revisions:
%   Date           Programmers               Description of change
%   ====        =================            =====================
%  06/15/17    Michael Nunez, Michelle Cheung, Sean O'Reilly-Jones pdmexp6
%  06/19/17      Michael Nunez              Flip unsure and sad
%  06/22/17      Michael Nunez             Default refresh rate: 144 Hz
%  07/05/17      Michael Nunez 
%       Default refresh rate: 120 Hz, gamma corrrection, photocell position
%  07/07/17      Michael Nunez                 Use of only two photocells
%  07/12/17      Michael Nunez           Fixation text, simplifying changes
%  07/13/17      Michael Nunez                 Display block order
%  07/14/17      Michael Nunez      Caffeine question, remove example images
%                                     Automatic block order config
%  07/17/17                          output() fix, blocktype fix
%  07/18/17      Michael Nunez        Fix number of blocks to subject input
%  07/19/17      Michael Nunez         Record condition per trial
%                                      Change given_feedback to [2 1 0]
%                         If number of blocks is less than 6, skip resting state
%                                      Track noise flicker
%  08/02/17      Michael Nunez     Fix session 1 blocktype generation error

%% Initial PTB3 Code

PsychJavaTrouble;   %make GetChar work (hack fix; needed?)

AssertOpenGL; %Issue warning if PTB3 with non-openGL used

% if strcmp(whichcmp,'k')
%     Screen('Preference', 'SkipSyncTests', 1);
% end

if ~IsLinux
    error('This program was written to run on Ubuntu Linux.');
end

%Find port for reponse box
devs = dir('/dev/ttyUSB*');
if isempty(devs)
    error('The Cedrus Response Box could not be found!');
end
    
if length(devs) > 1
    !dmesg | grep tty
    warning('There are multiple devices that could be the Cedrus Response Box!\n Find "ttyUSB*" in the above output on the same line as "FTDI USB Serial Device"');
end

if nargin > 1
    error('Too many function inputs.');
end
%% Experimenter Prompt

output = struct();

%Inputs Prompt and Output Setup
%Experimenter Prompt
Screenres = get(0,'Screensize');

prompt1={'Subject Number (must begin with letter):','Session Number:','Training session?','Cedrus Port [ttyUSB0 ...]:'};

def1={'SZZ_test','0','1',sprintf('/dev/%s',devs(1).name)};
studytitle='PDM Experiment 7';


lineNo=1;
answer=inputdlg(prompt1,studytitle,lineNo,def1);
%Subject Number
subnum = answer{1};
%ExpSession Number
sesnum = str2num(answer{2});
if isempty(sesnum)
    error('Please enter an appropriate session number (1 or greater)!');
end
output.sesnum = sesnum;
%Training session 'Real session - 0, Training - 1, 
training = str2num(answer{3});

%Window Pointer / 'Home Screen'.  0 - the primary monitor; 1 - the secondary monitor.
whichScreen = 0;
%Screen resolution on the x-axis
xres = Screenres(3);
output.xres = xres;
%Screen resolution on the y-axis
yres = Screenres(4);
output.yres = yres;
%This should be the same as the Refresh Rate shown in the Display
%Properties on the computer.  Always check before running the experiment to
%match flicker frequency.
%This code is currently set up to only handle multiples of 60 fps.
refrate = 120;
realrefrate = Screen(0,'FrameRate');
if refrate ~= Screen(0,'FrameRate')
    error(['The real screen refresh rate is set to ',num2str(realrefrate),...
       'Hz while the proposed screen refresh rate is ',num2str(refrate),'Hz.']);
end
output.refrate = refrate;
%Noise frequency (Hz)
noisehz = 40;
output.noisehz = noisehz;
if round(refrate/noisehz) ~= refrate/noisehz
    error('The noise frequency should be divisor of the refresh rate.');
end

%Gabor flicker frequency (Hz)
flickerhz = 30;
output.flickerhz = flickerhz;
if round(refrate/2/flickerhz) ~= refrate/2/flickerhz
    error('The gabor flicker frequency should be divisor of half the refresh rate.');
end

%Trials per block
if training
    output.tperb = 4;
else
    output.tperb = 80;
end
if round(output.tperb/2) ~= output.tperb/2
    error('There should be an even number of trials per block.');
end

%Initialize block
block = 1;

%SNR
snrs = .6;
output.snrs = snrs;

%signal luminance
lowboundsnr = min(snrs);
slum = lowboundsnr/(1 + lowboundsnr);

%Gabor spatial frequencies (cycles per degree at 57 cm)
gaborspf = [2.4 2.6];

%Noise spatial frequency (cycles per degree at 57 cm)
noisespf = 10;

%Radius of fixation spot (degrees visual angle at 57 cm)
fixrad = .10;

%Cedrus Handle
cport = answer{4};
try
    chandle = CedrusResponseBox('Open',cport);
catch me1
    rethrow(me1);
    fprintf('Cedrus port may need a ''chmod 777'' if you''re getting permission issues.\n');
end

%% Load macro experiment information


%Define block types
temptype{1} = [1 2 3 4 5 6];
temptype{2} = [1 2 3 6 5 4];
temptype{3} = [3 2 1 4 5 6];
temptype{4} = [3 2 1 6 5 4];
temptype{5} = [4 5 6 1 2 3];
temptype{6} = [6 5 4 1 2 3];
temptype{7} = [4 5 6 3 2 1];
temptype{8} = [6 5 4 3 2 1];
randtype = randperm(numel(temptype));

if training
    blocktype = [3 2 1 6 5 4];
else
    blocktype = temptype{randtype(1)};
end


if ~exist([pwd,'/exp7behav'],'dir')
   mkdir('exp7behav');
end

if exist('exp7behav/macroinfo.mat','file') && sesnum ~= 1
    fprintf('Loading macro experiment information to obtain blocktype...\n');
    macroinfo = load('exp7behav/macroinfo.mat');
    subsesfield = sprintf('%s_ses%d',subnum,sesnum-1);
    if isfield(macroinfo,subsesfield)
        lastses = macroinfo.(subsesfield).blocktype;
       if strcmp(num2str(lastses),num2str(temptype{1}))
               blocktype = temptype{8};
       elseif strcmp(num2str(lastses),num2str(temptype{2}))
               blocktype = temptype{7};
       elseif strcmp(num2str(lastses),num2str(temptype{3}))
               blocktype = temptype{6};
       elseif strcmp(num2str(lastses),num2str(temptype{4}))
               blocktype = temptype{5};
       elseif strcmp(num2str(lastses),num2str(temptype{5}))
               blocktype = temptype{4};
       elseif strcmp(num2str(lastses),num2str(temptype{6}))
               blocktype = temptype{3};
       elseif strcmp(num2str(lastses),num2str(temptype{7}))
               blocktype = temptype{2};
       elseif strcmp(num2str(lastses),num2str(temptype{8}))
               blocktype = temptype{1};
       else
       end
    end
end

%% Subject Prompt
prompt2={'Block Order (1r.6 2r.9 3r1.5 4l.6 5l.9 6l1.5)',...
    'What is your gender?',...
    'Age:',...
    'Do you consider yourself right handed, left handed, or both?  (''r'',''l'', or''b'')',...
    'Visual acuity result? (''20/30'' or ''20/35'' or ''20/20'')',...
    'What is your EEG cap size? (''Large'', ''Medium'', or ''Small'')',...
    'Approximate time since last caffeinated beverage in hours:',...
    'Do you have any personal or family history of epilepsy? (''y'' or ''n'')'
    };
def2={sprintf('[%s]',num2str(blocktype)),'','','','','','',''};
demographtitle='Subject Demographics';
lineNo=1;
subdemo=inputdlg(prompt2,demographtitle,lineNo,def2);
switch subdemo{8}
    case 'n'
    otherwise
        CedrusResponseBox('Close',chandle);
        error('You have indicated that you have a personal or family history of epilepsy. This experiment involves a fast flickering image. It is recommended that you NOT participate in this study due to a possible risk of seizure.  Please discuss your options with the experimenters.');
end
output.gender = subdemo{2};
output.age = str2num(subdemo{3});
output.hand = subdemo{4};
output.vision = subdemo{5};
output.capsize = subdemo{6};
output.caffeine = subdemo{7};

blocktype=str2num(subdemo{1});
output.blocktype = blocktype;

%Number of blocks
blocknum = length(blocktype);
if blocknum < 6
    skipeyes = 1;
else
    skipeyes = 0;
end

%Number of Trials
trialnum = blocknum*output.tperb;

%% Code

%Get date and time that the session begins
output.date = date;
output.start_time = clock;
    
%number of rows and columns of image
nCols = 1000;
nRows = 1000;

%Initialize estimated accuracy and speed_cutoff vectors
estcorrect = nan(1,trialnum);
speed_cutoff = nan(1,trialnum);
given_feedback = nan(1,trialnum);
condition = nan(1,trialnum);

%Keyboard keypress variables
advancechar = ' ';
escapechar = 27;

%Colors
txtcolor = round([.5 0 0]*255); %Red
black = [0 0 0];
white = [255 255 255];
gray = 255*[.5 .5 .5];
%darkgray = 255*[.25 .25 .25];
blackwhite{1} = black;
blackwhite{2} = white;

% Load fonts
myfont = '-bitstream-courier 14 pitch-bold-i-normal--0-0-0-0-m-0-ascii-0';
fontsize = 26;

%Define photocell placement
photorect = [0 0 100 100];
% pRect(1,:) = CenterRectOnPoint(photorect, 50, 50);
% pRect(2,:) = CenterRectOnPoint(photorect,xres - 50, 50);
% pRect(3,:) = CenterRectOnPoint(photorect,xres - 50, yres - 50);
% pRect(4,:) = CenterRectOnPoint(photorect, 50, yres - 50);
pRect(1,:) = CenterRectOnPoint(photorect,xres - 50, yres - 50);
pRect(2,:) = CenterRectOnPoint(photorect, 50, yres - 50);

fullscreen = CenterRectOnPoint([0 0 xres yres],round(xres/2),round(yres/2));

%Define button instruction placement
instRow = 374;
instCol = 504;
leftinstruct = CenterRectOnPoint(SetRect(0,0, instCol, instRow), 400, 900);
middleinstruct = CenterRectOnPoint(SetRect(0,0, instCol, instRow), 1000, 900);
rightinstruct = CenterRectOnPoint(SetRect(0,0, instCol, instRow), 1600, 900);

%Flush Cedrus Events
CedrusResponseBox('FlushEvents',chandle);


%The following TRY, CATCH, END statement ends psychtoolbox if an error
%occurs
try
    %Open a fullscreen window on the first screen with black background
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'UseVirtualFramebuffer');
    PsychImaging('AddTask','FinalFormatting','DisplayColorCorrection','SimpleGamma'); %apply power-law based gamma correction: http://docs.psychtoolbox.org/PsychColorCorrection
    [wptr,windowRect] = PsychImaging('OpenWindow', whichScreen,gray);
    PsychGPUControl('FullScreenWindowDisablesCompositor', 1);

    %sets size of gabor field that will be pasted onto Screen
    imageRect=SetRect(0,0,nCols,nRows);
    destRect=CenterRect(imageRect,windowRect);

    %Creates a window of a black screen with gray circle and fixation spot
    fiximage = makefixation([],fixrad);
    fiximage(fiximage == 1) = gray(1); %Gamma correction, gray
    fiximage(isnan(fiximage)) = gray(1);
    fiximage(fiximage == 0) = black(1);
    fixglind = Screen('MakeTexture',wptr,fiximage);
    blackfix = (makefixation([],fixrad) == 0);
    
    %This vector defines the noise frequency for our image
    noiseflic = [];
    for i=1:ceil(4*noisehz)
        noiseflic = [noiseflic 1 zeros(1,(round(refrate/noisehz)- 1))];
    end
    noiseonfind = find(noiseflic);
    
    %This vector defines the Gabor flicker frequency for our image
    gaborflic = [];
    for i=1:ceil(4*flickerhz)
        gaborflic = [gaborflic 2*ones(1,round(refrate/2/flickerhz)) ones(1,round(refrate/2/flickerhz))];
    end
    
    %Set seed based on the time.
    output.seed = round(sum(100*clock));
    % rng('default');  Backwards compatible with older MATLAB versions
    rng(output.seed);

    %Gabor will not be shown for the first 500ms to 1000ms of the trial
    numframes = round(refrate/noisehz);
    minframe = round(.5*refrate/numframes)*numframes;
    maxframe = round(refrate/numframes)*numframes;
    posframes = minframe:numframes:maxframe;
    trialnframes = posframes(randi(length(posframes),1,trialnum));
    output.noisetimes = trialnframes/refrate;
    
    %Inter-trial interval, 1500ms to 2000ms
    output.intertrial = 1.5 + rand(1,trialnum)*.5;
    
    
    %Define SNR vector, ensure even cell counts
    snrvec = [];
    for b=1:blocknum
        blkvec = [];
        for n=1:length(snrs)
            blkvec = [blkvec snrs(n)*ones(1,output.tperb/length(snrs))];
        end
        snrvec = [snrvec blkvec(randperm(numel(blkvec)))];
    end
    output.snrvec = snrvec;

    %%Define 6 block prompts
    % Action selection blocks
    blockprompt{1} = ['In this block use your RIGHT hand to respond within .6 seconds.'];
    blockprompt{2} = ['In this block use your RIGHT hand to respond within .9 seconds.'];
    blockprompt{3} = ['In this block use your RIGHT hand to respond within 1.5 seconds.'];
    blockprompt{4} = ['In this block use your LEFT hand to respond within .6 seconds.'];
    blockprompt{5} = ['In this block use your LEFT hand to respond within .9 seconds.'];
    blockprompt{6} = ['In this block use your LEFT hand to respond within 1.5 seconds.'];

    
    cut = 0; %Counter for ESC
    
    %Calculate the number of frames in a cycle of an image flicker
    numCycleFrames = trialnframes + ceil(refrate*1.5) + ceil(refrate*rand(1,trialnum)*.5);

    %Find the indices of 600, 900, and 1500 ms after the beginning of the noise onset dependent upon the blocktype
    speednframes{1} = trialnframes + ceil(refrate*.6);
    speednframes{2} = trialnframes + ceil(refrate*.9);
    speednframes{3} = trialnframes + ceil(refrate*1.5);
    speednframes{4} = trialnframes + ceil(refrate*.6);
    speednframes{5} = trialnframes + ceil(refrate*.9);
    speednframes{6} = trialnframes + ceil(refrate*1.5);

    %Output stimulus display time in seconds
    output.stimtime = (numCycleFrames/refrate);
    
    %Initialize recording of trialflic
    output.trialflic = cell(1,trialnum);
    output.noiseflic = cell(1,trialnum);
    
    output.nlum = slum./snrvec; %noise luminance
    
    Screen('TextFont',wptr,'Arial');
    Screen('TextSize',wptr,18);
    ShowCursor(0);	% arrow cursor
    sessiontext = 'Loading images...';
    sessiontext2 = sprintf('The experiment will begin soon!');
    sessiontext4 = blockprompt{blocktype(block)};
    if ~training
        sessiontext3 = sprintf('Please keep your eyes fixated on the dot. Session %d has started! Good luck!',sesnum);
    else
        sessiontext3 = sprintf('Please keep your eyes fixated on the dot. A practice session has started! Good luck!',sesnum);
    end
    sessiontext5 = 'Let''s record your brain rhythms with your eyes closed for two minutes!';
    sessiontext6 = 'Let''s record your brain rhythms with your eyes fixated on the center dot for two minutes!';
    trialtext3 = 'Please wait for the experimenter';
    trialtext4 = 'Press the middle button to continue';
    
    HideCursor;

    if training | skipeyes
        Screen('DrawTexture',wptr,fixglind,[],destRect);
        Screen('TextSize', wptr, fontsize);
        Screen('TextFont', wptr, myfont);
        DrawFormattedText(wptr, [sessiontext,'\n\n'], 'center', 'center', txtcolor);
        Screen(wptr,'FillRect',black,pRect'); 
        Screen('Flip',wptr);
    else

        %Display the eyes closed text screen
        Screen('DrawTexture',wptr,fixglind,[],destRect);
        Screen('TextSize', wptr, fontsize);
        Screen('TextFont', wptr, myfont);
        DrawFormattedText(wptr, [sessiontext2,'\n\n',sessiontext5], 'center', 'center', txtcolor);
        Screen(wptr,'FillRect',black,pRect'); 
        Screen('Flip',wptr);

        %Wait for spacebar
        FlushEvents('keyDown');
        [char,when] = GetChar; %Wait for keypress to continue
        notspace=1;
        while notspace
            switch char
                case ' '
                    notspace =0;
                case escapechar %Escape from experiment
                    notspace =0;
                    %RestoreScreen(whichScreen);
                    ShowCursor;
                    Screen('CloseAll');
                    warning on all;
                    CedrusResponseBox('Close',chandle);
                    return
                otherwise
                    [char,when] = GetChar; %Wait for keypress to continue
                    notspace =1;
            end
        end

        %Display the eyes closed text screen, wait for subject to press the middle button
        Screen('DrawTexture',wptr,fixglind,[],destRect);
        Screen('TextSize', wptr, fontsize);
        Screen('TextFont', wptr, myfont);
        DrawFormattedText(wptr, [trialtext4,'\n\n',sessiontext5], 'center', 'center', txtcolor);
        Screen(wptr,'FillRect',black,pRect'); 
        Screen('Flip',wptr);

        CedrusResponseBox('FlushEvents',chandle);
        notmiddle = 1;
        while notmiddle
            evt = CedrusResponseBox('WaitButtonPress',chandle);
            if evt.button == 5
                notmiddle = 0;
            end
        end
        CedrusResponseBox('FlushEvents',chandle);
        
        %Pause while recording eye closed data
        Screen(wptr,'FillRect',black,fullscreen); %black out the screen
        Screen(wptr,'FillRect',black,pRect'); 
        Screen('Flip',wptr);
    end

    %Track generation time
    tic;

    %Load hand position images and resize
%     boxleftleft = imread('example_images/pdmexp6/finger_placement/Box Left Finger on Low Freq.jpg');
%     boxleftleft = imresize(boxleftleft, [instRow instCol]);
%     textleftleft = Screen('MakeTexture',wptr,boxleftleft);
%     boxleftmiddle = imread('example_images/pdmexp6/finger_placement/Box Left Finger on Middle.jpg');
%     boxleftmiddle = imresize(boxleftmiddle, [instRow instCol]);
%     textleftmiddle = Screen('MakeTexture',wptr,boxleftmiddle);
%     boxleftright = imread('example_images/pdmexp6/finger_placement/Box Left Finger on High Freq.jpg');
%     boxleftright = imresize(boxleftright, [instRow instCol]);
%     textleftright = Screen('MakeTexture',wptr,boxleftright);
%     boxrightleft = imread('example_images/pdmexp6/finger_placement/Box Right Finger on Low Freq.jpg');
%     boxrightleft = imresize(boxrightleft, [instRow instCol]);
%     textrightleft = Screen('MakeTexture',wptr,boxrightleft);
%     boxrightmiddle = imread('example_images/pdmexp6/finger_placement/Box Right Finger on Middle.jpg');
%     boxrightmiddle = imresize(boxrightmiddle, [instRow instCol]);
%     textrightmiddle = Screen('MakeTexture',wptr,boxrightmiddle);
%     boxrightright = imread('example_images/pdmexp6/finger_placement/Box Right Finger on High Freq.jpg');
%     boxrightright = imresize(boxrightright, [instRow instCol]);
%     textrightright = Screen('MakeTexture',wptr,boxrightright);

    %Load Feedback Images and resize
    posImage=imread('Happyface.jpeg');
    posImage = posImage(:,:,1);
    unsImage=imread('unsureface.jpg');
    unsImage = unsImage(:,:,1);
    negImage=imread('sadface.jpg');
    negImage = negImage(:,:,1);

    %fix Feedback Images to make adjustable Black and Gray Image
    posImage(posImage<10)=black(1);
    posImage(posImage>10)=gray(1);
    unsImage(unsImage<10)=black(1);
    unsImage(unsImage>10)=gray(1);
    negImage(negImage<10)=black(1);
    negImage(negImage>10)=gray(1);
    
   
    % Resize Feedback Images
    posImage=imresize(posImage,[nRows nCols]);
    unsImage=imresize(unsImage,[nRows nCols]);
    negImage=imresize(negImage,[nRows nCols]); 
    
    %Add fixation cross to Images
    posImage(blackfix)=black(1);
    unsImage(blackfix)=black(1);
    negImage(blackfix)=black(1);
    
    %Write Images as textures
    posTex=Screen('MakeTexture',wptr,posImage);
    unsTex=Screen('MakeTexture',wptr,unsImage);
    negTex=Screen('MakeTexture',wptr,negImage);
    
    %Generate Gabor and noise images for all blocks
    noisecount = noisehz*4;
    %pregenexp5(noisecount,noisespf,trialnum,gaborspf,fixrad);
    [noises, allgabors, allspfs, allrandrot] = genexp5(noisecount,noisespf,trialnum,gaborspf,fixrad);
    
    %Concatenate Gabor patches from high and low distributions
    gabors = cat(3,allgabors{2},allgabors{1});
    spfs = [allspfs{2} allspfs{1}];
    randrots = [allrandrot{2} allrandrot{1}];
    
    %Index of spatial frequencies for Gabor images
    temp = [ones(1,(trialnum/2))*2 ones(1,(trialnum/2))*1];
%     if sesnum <= 1
%         repthisperm = randperm(trialnum,output.tperb); %Create first random permulation of length: trials per block
%         spfperm = [];
%         for b=1:blocknum
%             spfperm = [spfperm repthisperm(randperm(output.tperb))]; %If it is the training session, repeat exact distribution draws
%         end
%     else
        spfperm = randperm(trialnum);
%     end
    output.highlow = temp(spfperm); %Index if it was draw from high or low distribution
    clear temp;
    gabors = gabors(:,:,spfperm);
    output.spfs = spfs(spfperm);
    
    altperm = randperm(trialnum);
    output.randrots = randrots(altperm);
    
    %Locate fixation spot
    blackfix = (makefixation([],fixrad) == 0);
    
    %Locate border spot
    %darkfix = (makefixborder([],fixrad) == 0);
    
    %Change NaN to gray
    noises(isnan(noises)) = 0;
    gabors(isnan(gabors)) = 0;
    
    %Random order of noise images
    whichnoises = nan(trialnum,noisehz*4);
    for t=1:trialnum
        whichnoises(t,:) = randperm(noisehz*4);
    end
    
    %Index for unique SNRs
    [thesesnrs,~,usnrs] = unique(snrvec); 
    
    %Multiply each noise matrix by its luminance
    trialnoises = nan(size(noises,1),size(noises,2),numel(snrs),noisehz*4);
    for b=1:numel(snrs)
        splum = slum./thesesnrs(b); %noise luminance
        trialnoises(:,:,b,:) = splum.*noises;
    end
    thesenoises = 255*(trialnoises/2 + .5);
    %thesenoises(repmat(repmat(darkfix,1,1,b),1,1,1,noisehz*4)) = darkgray(1); %Recreate border
    thesenoises(repmat(repmat(blackfix,1,1,b),1,1,1,noisehz*4)) = black(1); %Recreate fixation spot
    noisescreen = nan(1,size(noises,3));
    
    %Create noise images as OpenGL textures
    for n=1:(noisehz*4)
        for b=1:numel(snrs)
            noisescreen(b,n) = Screen('MakeTexture',wptr,thesenoises(:,:,b,n)); %The square root is in order to account for monitor gamma. That is, the monitor approximately squares the input stimulus color value
        end
    end
    clear thesenoises noises
    
    %Save generation time
    output.gentime = toc;

    if ~training & ~skipeyes
        %Pause while recording eye closed data
        pause(120-output.gentime);

        %Display the eyes open text screen
        Screen(wptr,'FillRect',gray,fullscreen); %black out the screen
        Screen('DrawTexture',wptr,fixglind,[],destRect);
        Screen('TextSize', wptr, fontsize);
        Screen('TextFont', wptr, myfont);
        DrawFormattedText(wptr, [sessiontext2,'\n\n',sessiontext6], 'center', 'center', txtcolor);
        Screen(wptr,'FillRect',black,pRect'); 
        Screen('Flip',wptr);

        %Wait for spacebar
        FlushEvents('keyDown');
        [char,when] = GetChar; %Wait for keypress to continue
        notspace=1;
        while notspace
            switch char
                case ' '
                    notspace =0;
                case escapechar %Escape from experiment
                    notspace =0;
                    %RestoreScreen(whichScreen);
                    ShowCursor;
                    Screen('CloseAll');
                    warning on all;
                    CedrusResponseBox('Close',chandle);
                    return
                otherwise
                    [char,when] = GetChar; %Wait for keypress to continue
                    notspace =1;
            end
        end

        %Display the eyes closed text screen, wait for subject to press the middle button
        Screen('DrawTexture',wptr,fixglind,[],destRect);
        Screen('TextSize', wptr, fontsize);
        Screen('TextFont', wptr, myfont);
        DrawFormattedText(wptr, [trialtext4,'\n\n',sessiontext6], 'center', 'center', txtcolor);
        Screen(wptr,'FillRect',black,pRect'); 
        Screen('Flip',wptr);

        CedrusResponseBox('FlushEvents',chandle);
        notmiddle = 1;
        while notmiddle
            evt = CedrusResponseBox('WaitButtonPress',chandle);
            if evt.button == 5
                notmiddle = 0;
            end
        end
        CedrusResponseBox('FlushEvents',chandle);
        
        %Pause while recording eye closed data
        Screen('DrawTexture',wptr,fixglind,[],destRect);
        Screen(wptr,'FillRect',black,pRect'); 
        Screen('Flip',wptr);
        pause(120);
    end

    %Display second text screen
    Screen('DrawTexture',wptr,fixglind,[],destRect);
    Screen('TextSize', wptr, fontsize);
    Screen('TextFont', wptr, myfont);
    % Screen('DrawText',wptr, sessiontext2,(xres - length(sessiontext2))/2,yres*(5/12),txtcolor);
    % Screen('DrawText',wptr, sessiontext4,(xres - length(sessiontext4))/2,yres*(5/12) + 32,txtcolor);
    DrawFormattedText(wptr, [sessiontext2,'\n\n',sessiontext4], 'center', 'center', txtcolor);
%     if blocktype(block) <= 3
%         Screen('DrawTexture',wptr,textrightleft,[],leftinstruct);
%         Screen('DrawTexture',wptr,textrightmiddle,[],middleinstruct);
%         Screen('DrawTexture',wptr,textrightright,[],rightinstruct);
%     elseif blocktype(block) > 3
%         Screen('DrawTexture',wptr,textleftleft,[],leftinstruct);
%         Screen('DrawTexture',wptr,textleftmiddle,[],middleinstruct);
%         Screen('DrawTexture',wptr,textleftright,[],rightinstruct);
%     end
    Screen(wptr,'FillRect',black,pRect'); 
    Screen('Flip',wptr);
    
    %Wait for spacebar
    FlushEvents('keyDown');
    [char,when] = GetChar; %Wait for keypress to continue
    notspace=1;
    while notspace
        switch char
            case ' '
                notspace =0;
            case escapechar %Escape from experiment
                notspace =0;
                %RestoreScreen(whichScreen);
                ShowCursor;
                Screen('CloseAll');
                warning on all;
                CedrusResponseBox('Close',chandle);
                return
            otherwise
                [char,when] = GetChar; %Wait for keypress to continue
                notspace =1;
        end
    end
    
    %Wait for subject to press the middle button
    Screen('DrawTexture',wptr,fixglind,[],destRect);
    Screen('TextSize', wptr, fontsize);
    Screen('TextFont', wptr, myfont);
    DrawFormattedText(wptr, [trialtext4,'\n\n',sessiontext4], 'center', 'center', txtcolor);
%     if blocktype(block) <= 3
%         Screen('DrawTexture',wptr,textrightleft,[],leftinstruct);
%         Screen('DrawTexture',wptr,textrightmiddle,[],middleinstruct);
%         Screen('DrawTexture',wptr,textrightright,[],rightinstruct);
%     elseif blocktype(block) > 3
%         Screen('DrawTexture',wptr,textleftleft,[],leftinstruct);
%         Screen('DrawTexture',wptr,textleftmiddle,[],middleinstruct);
%         Screen('DrawTexture',wptr,textleftright,[],rightinstruct);
%     end
    Screen(wptr,'FillRect',black,pRect'); 
    Screen('Flip',wptr);
    
    CedrusResponseBox('FlushEvents',chandle);
    notmiddle = 1;
    while notmiddle
        evt = CedrusResponseBox('WaitButtonPress',chandle);
        if evt.button == 5
            notmiddle = 0;
        end
    end
    CedrusResponseBox('FlushEvents',chandle);
    
    %Display third text screen
    Screen('DrawTexture',wptr,fixglind,[],destRect);
    Screen('TextSize', wptr, fontsize);
    Screen('TextFont', wptr, myfont);
    % Screen('DrawText',wptr, sessiontext3,(xres - length(sessiontext3))/2,yres*(5/12),txtcolor);
    DrawFormattedText(wptr, [sessiontext3,'\n\n'], 'center', 'center', txtcolor);
    Screen(wptr,'FillRect',black,pRect'); 
    Screen('Flip',wptr);
    pause(2);
    
    %Pause before beginning of the block, pause for a second after text
    Screen('DrawTexture',wptr,fixglind,[],destRect);
    Screen(wptr,'FillRect',black,pRect'); 
    Screen('Flip',wptr);
    pause(1);
    
    %Initialize timer
    tic;
    for trials = 1:trialnum
      %printiter(trials);
      if ~cut %ESC key track
        %trialinbl = trials - output.tperb*(block-1);
        
        %Create image with specified snr
        trialgabor = slum*gabors(:,:,trials);
        %Shift domain to [0 1], convert to appropriate luminance
        %transformation for gaborimage color values, takes into account
        %monitor gamma
        if trials == 1
            bothscreen = nan(1,noisehz*4);
        end
        for n=1:noisehz*4
            thisimage = 255*( (trialgabor(:,:)+trialnoises(:,:,usnrs(trials),whichnoises(trials,n)) ...
               )/2 + .5);
            %thisimage(darkfix) = darkgray(1);
            thisimage(blackfix) = black(1);
            if trials > 1 %Clear former textures to save memory
                Screen('Close',bothscreen(n));
            end
            bothscreen(n) = Screen('MakeTexture',wptr,thisimage);
        end
        clear trialgabor fixx fixy trialflic
        
        trialflic = [ones(1,trialnframes(trials)) gaborflic];
        output.trialflic{trials} = trialflic;
        output.noiseflic{trials} = noiseflic;
          
        %Wait at least lboundwait seconds between trials
        lboundwait = output.intertrial(trials);
        output.elapsedtime(trials) = toc;
        if output.elapsedtime(trials) < lboundwait
            pause(lboundwait-output.elapsedtime(trials));
        end
        output.fixedtime(trials) = toc;
        
        CedrusResponseBox('FlushEvents',chandle);
        
        %Display rush loops (Rush is apparently obsolete in PTB3, test this)
        Priority(MaxPriority(wptr)); %New to PTB3
        
        %Loop 1: Noise interval for 500ms - 1000ms
        %Flush Cedrus box responses after the noise interval, does this effect the stimulus draw?
        %Response interval for 1000ms - 2000ms, accept responses
        noisenum = 0;
        bwswitch = 1;
        clearresp = 0;
        speed_acc = 0;
        evt = [];
        for i = 1:numCycleFrames(trials)
            if noiseflic(i)
                noisenum = noisenum + 1;
                bwswitch = mod(bwswitch,2) + 1; %Changes 1 to 2 and vica versa
            end
            if trialflic(i) == 2
                Screen('DrawTexture',wptr,bothscreen(noisenum),[],destRect);
                if clearresp == 0
                    clearresp = 1;
                end
            else
                Screen('DrawTexture',wptr,noisescreen(usnrs(trials),whichnoises(trials,noisenum)),[],destRect);
            end
            Screen(wptr,'FillRect',blackwhite{bwswitch},pRect(1,:)');
            Screen(wptr,'FillRect',blackwhite{trialflic(i)},pRect(2,:)');
            Screen('Flip',wptr);
            if clearresp == 1 %Flush Cedrus box responses after the noise interval
                CedrusResponseBox('FlushEvents',chandle);
                clearresp = 2;
            end
            if (i == speednframes{blocktype(block)}(trials))
                evt = CedrusResponseBox('GetButtons',chandle);
            
                if isempty(evt) %If there is no response given in time the subject is given unhappy feedback
                    speed_acc = 0;
                else
                    speed_acc = 1;
                end
                
            end
        end
    
        %Loop 2: Keep displaying black fixation spot (only) for 250ms to collect responses
        for frames = 1:round(refrate/4)
            Screen('DrawTexture',wptr,fixglind,[],destRect);
            Screen(wptr,'FillRect',black,pRect'); 
            Screen('Flip',wptr);
        end

        %Timer to calculate time between the last trial and the next
        tic;
        
        %Calculate trial accuracy, despite subject feedback
        if isempty(evt)
            evt = CedrusResponseBox('GetButtons',chandle);
        end
        if isempty(evt)
            correct = NaN;
        elseif (evt.button == 6 && output.highlow(trials) == 2) || ... %High frequency
                (evt.button == 4 && output.highlow(trials) == 1) %Low frequency
            correct = 1;
        else
            correct = 0;
        end
        estcorrect(trials) = correct;
        speed_cutoff(trials) = speed_acc;
        condition(trials) = blocktype(block);
 
        %Loop 3: Display feedback for 500 ms
        for frames = 1:round(refrate/2)
            if (correct == 1) & (speed_acc ~= 0)
                Screen('DrawTexture',wptr,posTex,[], destRect);
                given_feedback(trials) = 2;
            elseif (speed_acc == 0)
                Screen('DrawTexture',wptr,negTex,[],destRect);
                given_feedback(trials) = 0;
            else
                Screen('DrawTexture',wptr,unsTex,[],destRect);
                given_feedback(trials) = 1;
            end
            Screen('FillRect',wptr,white,pRect');
            Screen('Flip',wptr);
        end
        % Trial Number Display
        Screen('DrawTexture',wptr,fixglind,[],destRect);
        Screen('TextSize', wptr, 18);
        Screen('TextFont', wptr, myfont);
        Screen('DrawText', wptr, sprintf('B%i',block), 10, yres-140, black);
        Screen('DrawText', wptr, sprintf('T%i',trials), 10, yres-120, black);
        Screen(wptr,'FillRect',black,pRect');
        Screen('Flip',wptr);
        pause(.5)
        if ~cut
            if trials == trialnum
                %Show ending screen for 1 second
                peranswered = sum(given_feedback((trials-output.tperb+1):trials) > 0)/output.tperb;
                endtext = ['Done!  ',...
                    num2str(round(peranswered*100)),'% answered this block. Thank you for participating!'];                    
                Screen('DrawTexture',wptr,fixglind,[],destRect);
                Screen('TextSize', wptr, fontsize);
                Screen('TextFont', wptr, myfont);
                % Screen('DrawText',wptr, endtext,(xres - length(endtext))/2,yres*(5/12),txtcolor);
                DrawFormattedText(wptr, [endtext,'\n\n'], 'center', 'center', txtcolor);
                Screen(wptr,'FillRect',black,pRect'); 
                %Pause for 1 second
                pause(1);
                Screen('Flip',wptr);
                %Pause for 1 second
                pause(1);
                %Wait for spacebar to end program
                FlushEvents('keyDown');
                [char,~] = GetChar; %Wait for keypress to continue
                notspace=1;
                while notspace
                    switch char
                        case advancechar
                            notspace =0;
                        otherwise
                            [char,~] = GetChar; %Wait for keypress to continue
                            notspace =1;
                    end
                end
            elseif trials/output.tperb == round(trials/output.tperb)
                 %Take a break every 'output.tperb' trials and show ending Screens
                peranswered = sum(given_feedback((trials-output.tperb+1):trials) > 0)/output.tperb;  
                
                trialtext = sprintf('Block %i complete! %i%% answered this block. You may now take a break!',...
                    block, round(peranswered*100));
                
     
                block = block + 1;
                 

                Screen('DrawTexture',wptr,fixglind,[],destRect);
                Screen('TextSize', wptr, fontsize);
                Screen('TextFont', wptr, myfont);
                trialtext2 = blockprompt{blocktype(block)};
                DrawFormattedText(wptr, [trialtext,'\n\n',trialtext2,'\n\n',trialtext3], 'center', 'center', txtcolor);
%                 if blocktype(block) <= 3
%                     Screen('DrawTexture',wptr,textrightleft,[],leftinstruct);
%                     Screen('DrawTexture',wptr,textrightmiddle,[],middleinstruct);
%                     Screen('DrawTexture',wptr,textrightright,[],rightinstruct);
%                 elseif blocktype(block) > 3
%                     Screen('DrawTexture',wptr,textleftleft,[],leftinstruct);
%                     Screen('DrawTexture',wptr,textleftmiddle,[],middleinstruct);
%                     Screen('DrawTexture',wptr,textleftright,[],rightinstruct);
%                 end
                Screen(wptr,'FillRect',black,pRect'); 
                %Pause for 1 second
                pause(1);
                Screen('Flip',wptr);
                
                %Wait for spacebar
                FlushEvents('keyDown');
                [char,~] = GetChar; %Wait for keypress to continue
                notspace=1;
                while notspace
                    switch char
                        case advancechar
                            notspace =0;

                            %Wait for subject to press the middle button
                            Screen('DrawTexture',wptr,fixglind,[],destRect);
                            Screen('TextSize', wptr, fontsize);
                            Screen('TextFont', wptr, myfont);
                            DrawFormattedText(wptr, [trialtext,'\n\n',trialtext2,'\n\n',trialtext4], 'center', 'center', txtcolor);
%                             if blocktype(block) <= 3
%                                 Screen('DrawTexture',wptr,textrightleft,[],leftinstruct);
%                                 Screen('DrawTexture',wptr,textrightmiddle,[],middleinstruct);
%                                 Screen('DrawTexture',wptr,textrightright,[],rightinstruct);
%                             elseif blocktype(block) > 3
%                                 Screen('DrawTexture',wptr,textleftleft,[],leftinstruct);
%                                 Screen('DrawTexture',wptr,textleftmiddle,[],middleinstruct);
%                                 Screen('DrawTexture',wptr,textleftright,[],rightinstruct);
%                             end
                            Screen(wptr,'FillRect',black,pRect'); 
                            Screen('Flip',wptr);
                            
                            CedrusResponseBox('FlushEvents',chandle);
                            notmiddle = 1;
                            while notmiddle
                                evt = CedrusResponseBox('WaitButtonPress',chandle);
                                if evt.button == 5
                                    notmiddle = 0;
                                end
                            end
                            CedrusResponseBox('FlushEvents',chandle);

                            Screen('DrawTexture',wptr,fixglind,[],destRect);
                            Screen('TextSize', wptr, fontsize);
                            Screen('TextFont', wptr, myfont);
                            % Screen('DrawText',wptr,sprintf('Block %i of the experiment has started! Good luck!',block),(xres - length(sessiontext3))/2,yres*(5/12),txtcolor);
                            DrawFormattedText(wptr, [sprintf('Please keep your eyes fixated on the dot. Block %i has started! Good luck!',block),'\n\n'], 'center', 'center', txtcolor);
                            Screen(wptr,'FillRect',black,pRect'); 
                            Screen('Flip',wptr);
                            %Timer to calculate time between the last trial and the next
                            %Pause for 2 seconds
                            pause(2);
                            tic;
                        case escapechar %Escape from experiment and save current data (for experimenter)
                            notspace =0;
                            %RestoreScreen(whichScreen);
                            ShowCursor;
                            Screen('CloseAll');
                            output.ESC_time = clock;
                            output.estcorrect = estcorrect;
                            output.speed_cutoff = speed_cutoff;
                            output.given_feedback = given_feedback;
                            output.condition = condition;
                            %Organize data and time in a string
                            rightnow = clock;
                            rightnow = num2cell(rightnow)';
                            timestr = sprintf('_%i',rightnow{1:5});
                            eval([subnum,'_ExpSession',num2str(sesnum),'=output;']);
                            if ~strcmp(subnum,'SZZ_test')
                                if ~exist([pwd,'/exp7behav'],'dir')
                                    mkdir('exp7behav');
                                end
                                eval(['save(''exp7behav/',subnum,'_ses',num2str(sesnum),timestr,'.mat'',''-struct'', ''',subnum,'_ExpSession',num2str(sesnum),''');']);
                                if exist('exp7behav/macroinfo.mat','file')
                                    macroinfo = load('exp7behav/macroinfo.mat');
                                    subsesfield = sprintf('%s_ses%d',subnum,sesnum);
                                    macroinfo.(subsesfield) = output;
                                    save('exp7behav/macroinfo.mat','-struct','macroinfo');
                                end
                            end
                            warning on all;
                            CedrusResponseBox('Close',chandle);
                            return
                        otherwise
                            [char,when] = GetChar; %Wait for keypress to continue
                            notspace =1;
                    end
                end
                
                %Pause before beginning of the block, pause for a second after text
                Screen('DrawTexture',wptr,fixglind,[],destRect);
                Screen(wptr,'FillRect',black,pRect'); 
                Screen('Flip',wptr);
                pause(1);
            end
            
        end
      end  
    end
catch me
    
    fprintf('\n');
    %RestoreScreen(whichScreen);
    ShowCursor;
    Screen('CloseAll');
    output.error_time = clock;
    output.estcorrect = estcorrect;
    output.speed_cutoff = speed_cutoff;
    output.given_feedback = given_feedback;
    output.condition = condition;
    %Organize data and time in a string
    rightnow = clock;
    rightnow = num2cell(rightnow)';
    timestr = sprintf('_%i',rightnow{1:5});
    eval([subnum,'_ExpSession',num2str(sesnum),'=output;']);
    if ~strcmp(subnum,'SZZ_test')
        if ~exist([pwd,'/exp7behav'],'dir')
            mkdir('exp7behav');
        end
        eval(['save(''exp7behav/',subnum,'_ses',num2str(sesnum),timestr,'.mat'',''-struct'', ''',subnum,'_ExpSession',num2str(sesnum),''');']);
        if exist('exp7behav/macroinfo.mat','file')
            macroinfo = load('exp7behav/macroinfo.mat');
            subsesfield = sprintf('%s_ses%d',subnum,sesnum);
            macroinfo.(subsesfield) = output;
            save('exp7behav/macroinfo.mat','-struct','macroinfo');
        end
    end
    CedrusResponseBox('Close',chandle);
    %PsychPortAudio('Close', goodhand);
    %PsychPortAudio('Close', badhand);
    rethrow(me); %rethrow reproduces the original error, stored in the object 'me'
end

fprintf('\n');

%RestoreScreen(whichScreen);
ShowCursor;
Screen('CloseAll');

%Output time finished
output.finish_time = clock;

%Estimated accuracy
output.estcorrect = estcorrect;
output.speed_cutoff = speed_cutoff;
output.given_feedback = given_feedback;
output.condition = condition;

%Organize data and time in a string
rightnow = clock;
rightnow = num2cell(rightnow)';
timestr = sprintf('_%i',rightnow{1:5});
eval([subnum,'_ExpSession',num2str(sesnum),'=output;']);
if ~strcmp(subnum,'SZZ_test')
    if ~exist([pwd,'/exp7behav'],'dir')
        mkdir('exp7behav');
    end
    eval(['save(''exp7behav/',subnum,'_ses',num2str(sesnum),timestr,'.mat'',''-struct'', ''',subnum,'_ExpSession',num2str(sesnum),''');']);
    if exist('exp7behav/macroinfo.mat','file')
        macroinfo = load('exp7behav/macroinfo.mat');
        subsesfield = sprintf('%s_ses%d',subnum,sesnum);
        macroinfo.(subsesfield) = output;
        save('exp7behav/macroinfo.mat','-struct','macroinfo');
    end
end

warning on all;

CedrusResponseBox('Close',chandle);

%% ----------------------------------------------------------------------%%
function [noises, allgabors, allspfs, allrandrot] = genexp5(numnoise,noisespfrq,numgabor,gaborspfrq,radius)
%GENEXP% - generates images for pdmexp5
%
%Useage: 
%  >> genexp5(numnoise,noisespfrq,ngabor,gaborspfrq)
%
%Inputs:
%   numnoise - Number of noise images to generate
%
%   noisespfrq - Spatial frequency of pixelated visual noise (cycles per cm)
%
%   numgabor - Number of gabor images to pregenerate
%
%   gaborspfrq - Spatial frequencies of gabors (cycles per cm)
%
%   radius - Radius size of fixation spot (cycles per cm)

%% Code

%Gabor size
gaborsize = 10;

if round(numgabor/length(gaborspfrq)) ~= numgabor/length(gaborspfrq)
    error('Number of Gabor spatial frequencies must be a divisor of the number of Gabor images');
end

if numnoise > 0
    fprintf('Building noise images...\n');
end
noises = nan(1000,1000,numnoise);
% for m=1:numnoise
%     noises(:,:,m) = makepixeled([],radius,noisespfrq,[]);
% end

%%Combine two spatial frequencies to equally mask both high and low signal
%%stimuli
for m=1:numnoise
%     temphn = makepixeled([],radius,3,[]);
%     templn = makepixeled([],radius,2,[]);
%     noises(:,:,m) = (temphn + templn)/2;
      noises(:,:,m) = makenoise([],2,10,radius,[2 3]);
end


%%Make Gabor images
if ~exist([pwd,'/pregen/gabor'],'dir')
    mkdir('pregen/gabor');
end

if numgabor > 0
    fprintf('Building Gabor images...\n');
    for f=1:length(gaborspfrq)
        numimages = numgabor/length(gaborspfrq);
        gabors = nan(1000,1000,numimages);
        spf = nan(numimages,1);
        randrot = nan(numimages,1);
        for m=1:numimages
            %Draw from truncated normal distributions
            if f==1
%                 draw = 100;
%                 while draw > 0 %There is a 2.28% chance of the draw being less than 0
%                     draw = normrnd(-.1,.05);
%                 end
            draw = -.1;
            elseif f==2
%                 draw = -100;
%                 while draw < 0 %There is a 2.28% chance of the draw being more than 0
%                     draw = normrnd(.1,.05);
%                 end
            draw = .1;
            end
            spf(m) = draw + 2.5; %modes = 2.4 | 2.6
            randrot(m) = rand(1)*360;
            gabors(:,:,m) = makegabor([],gaborsize,randrot(m)*360,[],spf(m));
        end
        allgabors{f} = gabors;
        allspfs{f} = spf;
        allrandrot{f} = randrot;
    end
end

timeels = toc;
fprintf('The images took %3.2f minutes to generate and save.\n',timeels/60);