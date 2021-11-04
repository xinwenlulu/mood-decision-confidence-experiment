data.savefile               = 'test.mat'
disp(['Saving to: ', data.savefile]); save(data.savefile, 'data')

%% Set experiment settings

% set cogent settings
addpath 'Cogent2000v1.33\Toolbox' % Update this with where you put Cogent folder
%addpath 'N:\DesktopSettings\Desktop\images'
screenMode              = 0;                 % 0 for small window, 1 for full screen, 2 for second screen if attached
screenRes               = 3;                 % 1=640x480; 2=800x600; 3=1024x768
foreground              = [0.5 0.5 0.5];           % foreground colour 111 is white
data.background         = [0.5 0.5 0.5];           % background colour
fontName                = 'Arial';           % font parameters
fontSize                = 32;
nbits                   = 0;          

config_sound;
config_keyboard;
config_display(screenMode, screenRes, data.background, foreground, fontName, fontSize, 5, nbits);   % open graphics window
start_cogent

% Define experiment settings
data                  = struct;
data.nblocks = 3; % Because 3x4 = 12
conditions = [1 2 3]; % Reward distributions 1 = bad; 2 = goodvL; 3 = goodvH
data.pictures = {'bull.bmp';'chick.bmp';'crab.bmp';'fox.bmp';'hedgehog.bmp';'hippopotamus.bmp';'koala.bmp';'lemur.bmp';'pig.bmp';'tiger.bmp';'whale.bmp';'zebra.bmp'};
%data.pictures = {'flowers.png';'flowers.png';'flowers.png';'flowers.png';'flowers.png';'flowers.png';'flowers.png';'flowers.png';'flowers.png'};
data.background         = [0.5 0.5 0.5];           % background colour
data.ntrials = 6;

picprob = zeros(data.nblocks,6); %stimuli used and associated probabilities
index_side = randperm(length(data.pictures)); % Shuffle stimuli
for n=1:data.nblocks
     picprob(n,:) = [index_side(3*n-2:3*n) conditions];
end
data.picprob = picprob;
loadsound('Cheer2.wav', 15);
loadsound('Ding3.wav',16);
loadsound('down.wav',17);
loadsound('Boo.wav',18);

data.sidelocation             = [-200 200];
data.FontSizeText             = 32;
data.FontSizeTask             = 80;

% task timing
data.iti                      =  800;
data.iti_std                  =  300;
data.iti_instruction          =  800;
data.timeout                  = 1000;
data.timeout_instruction      = 1200;
data.optiontime               =    0;
data.confirmationtime         =  600;
data.outcometime              =  800;
data.time.outcome_delivery    =  600;
data.time.outcome_delivery_practice    = 0;


data.time.instructiontime     = 300;
data.time.delay               = 1000; %delay until outcome revealed
data.time.outcometime         = 700;
data.time.iti                 = 600; %min ITI time
data.time.maxiti              = 800; %max ITI time
data.time.ratewait            = 0; %before can make rating - no time limit

if exist('debugging', 'var')
    data.time.instructiontime     = 300;
    data.time.delay               = 50; %delay until outcome revealed
    data.time.outcometime         = 50;
    data.time.iti                 = 50; %min ITI time
    data.time.maxiti              = 60; %max ITI time
    data.time.ratewait            = 0; %before can make rating - no time limit
end
%end;
% hardware mapping
data.leftrep                  = 1;
data.rightrep                 = 4;

% task parameters
%data.nslots = 2;
data.npracticetrials = 2; % divisible by 2
data.ntrials = 3; % % divisible by 3
%data.nblocks = 4;

disp(['nblocks=,', num2str(data.nblocks), ' ntrials=', num2str(data.ntrials)])

%data.pictures = {'bull.bmp';'chick.bmp';'crab.bmp';'fox.bmp';'hedgehog.bmp';'hippopotamus.bmp';'koala.bmp';'lemur.bmp';'pig.bmp';'tiger.bmp';'whale.bmp';'zebra.bmp'};
data.stimDimensions = [300, 300]; 
data.rewardstim.picture = 'gemstone.bmp';
data.rewardstim.stimDimensions = [98, 93];
data.norewardstim.picture = 'Silhouettegemstone.bmp';
data.norewardstim.stimDimensions = [98, 93];

% Reward
data.endowment = 0; %0 gems to begin with
data.rewardsize = 10;
data.currentscore = data.endowment;
%%
%%%%COMMENT THIS OUT%%%%
blocknum = 1; wof = 0;;%data.ntrials=30; %data.picprob = [3 7 0.3 0.65]; %for testing
%%%%%%%%

cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); cgtext('+',0,0);
cgflip(data.background);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Drop columns with WoF outcomes
picprob_only = data.picprob(:,1:6); 
% Drop columns with WoF outcomes

block_picprob = picprob_only(blocknum, :);
% Create 3 pairs from the 3 stimuli in the current block
pairs = zeros(3,4);
comparison = zeros(3,4);
for n = 1:3 %3 pairs 
    if n ~= 3
        pairs(n,:) = [block_picprob(1,n:n+1) block_picprob(1,n+3:n+4)];
    else
        pair_index = [1 3 4 6];
        pairs(n,:) = block_picprob(pair_index);
    end
end

%%% Decide when to flip left / right
ntrial_perpair = data.ntrials/3; % 3 pairs
trialmtx = ones(data.ntrials,4);
for n = 1:ntrial_perpair
    trialmtx(3*n-2:3*n,:) = pairs;
end

for n = 1:data.ntrials
    if rand>0.5
        trialmtx(n,:) = trialmtx(n,[2 1 4 3]);
    end
end

%shuffle the trials in the trial matrix
trial_index = randperm(data.ntrials);
unshuffledmtx = trialmtx;
for n = 1:data.ntrials
    trialmtx(n,:) = unshuffledmtx(trial_index(n),:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%risktrials = 5:5:30;

%[~,~,~,mp] = cgmouse;
tstart = time;                      %get time in ms from first time call
nrating = 1; %nratingrisk = 1;


wait(data.time.iti);

%columns: 1 trial num, 2-3 stim identities (0 forced choice, safe value if
%>1), 4-5 stim probs, 6 button pressed (1 or 2), 7 is stimulus chosen, 8
%outcome (in points), 9 is rt, 10-12 happiness rating, 13 is WoF (0 if
%none, otherwise outcome), 14-18 trial timing info, 19-22 happy rating
%info, 23-36 risky trial info
experiment_data = zeros(data.ntrials,23);
experiment_data(:,1) = 1:data.ntrials;
experiment_data(:,13) = wof; %wheel of fortune outcome %%%%NEW%%%%
experiment_data(:,2:5) = trialmtx; %zero anything that is empty - 0 stimulus and 0 prob
%%%%%%%%%%%%%%%%%%%%%%%%

%Reward distributions 1 = bad; 2 = goodvL; 3 = goodvH
bad = [0 2 2 2 3 3 3 4 4 4 6];
goodvL = [4 5 5 6 6 6 6 6 7 7 8];
goodvH = [0 3 4 5 5 6 7 7 8 9 14];

%predetermine rewards for each option
rewardmL = repmat(bad,1,ceil(data.ntrials*2/3/size(bad,2)));
rewardmHvL = repmat(goodvL,1,ceil(data.ntrials*2/3/size(goodvL,2)));
rewardmHvH = repmat(goodvH,1,ceil(data.ntrials*2/3/size(goodvH,2)));

% shuffle the rewards in each trial
rewardmL = rewardmL(randperm(length(rewardmL)));
rewardmHvL = rewardmHvL(randperm(length(rewardmHvL)));
rewardmHvH = rewardmHvH(randperm(length(rewardmHvH)));

chosemL = 1; chosemHvL = 1; chosemHvH = 1; %how many of each chosen so far %%%%NEW%%%%

v = belief_probe(1, getkeymap)

for itrial = 1:data.ntrials
    disp(experiment_data(:,1:14));
    allowleftclick = -1; allowrightclick = -1; % by default, no mouse presses are defined as acceptable
    stimDimensions = data.stimDimensions;
    
    % check if free choice, forced left, or forced right trial
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); cgtext('+',0,0);
    if experiment_data(itrial,2)>0,
        allowleftclick=1; % allow click left
        pic = data.pictures{experiment_data(itrial,2)};
        cgloadbmp(1,pic,stimDimensions(1),stimDimensions(2));
        cgpencol(0,0,0);
        cgrect(420,350,200,100);
        cgpencol(1,1,1)
        cgfont('Arial',data.FontSizeText)
        cgtext('Your Gems',420,365);
        cgtext(num2str(data.currentscore+experiment_data(itrial,23)),420,330);
        cgrect(-200,0,stimDimensions(1),stimDimensions(2));
        cgdrawsprite(1,-200,0);
    end
    if experiment_data(itrial,3)>0,
        allowrightclick=4; % allow click right
        pic = data.pictures{experiment_data(itrial,3)};
        cgloadbmp(2,pic,stimDimensions(1),stimDimensions(2));
        cgpencol(0,0,0);
        cgrect(420,350,200,100);
        cgpencol(1,1,1)
        cgfont('Arial',data.FontSizeText)
        cgtext('Your Gems',420,365);
        cgtext(num2str(data.currentscore),420,330);
        cgrect(200,0,stimDimensions(1),stimDimensions(2));
        cgdrawsprite(2,200,0)
    end
    cgflip(data.background);
    
    % get the choice
    tstart2 = time; experiment_data(itrial,14) = (tstart2 - tstart)/1000; %time since start of this task - %%%%NEW%%%%edited
    [~,~,~,mp]=cgmouse; mp=0; %no mouse press yet
    while (mp ~=allowleftclick) && (mp ~=allowrightclick), %1 from leftmouse, 4 for rightmouse - if a free choice trial
        [~,~,~, mp] = cgmouse;
    end;
    rt = (time - tstart2)/1000; experiment_data(itrial,15) = (time - tstart)/1000; %rt in s %%%%NEW%%%%edited
    choicemade = sqrt(mp);
    experiment_data(itrial,6) = choicemade; %which button pressed (1 left, 2 right) or 0 for none
    
    
    experiment_data(itrial,7) = experiment_data(itrial,1+choicemade); %%%%NEW%%%%
    if experiment_data(itrial,3+choicemade) == 1 %if chose the badmL option %%%%NEW%%%%
        experiment_data(itrial,8) = rewardmL(chosemL); chosemL=chosemL+1; %%%%NEW%%%%
    elseif experiment_data(itrial,3+choicemade) == 2 %if chose the goodmL option
        experiment_data(itrial,8) = rewardmHvL(chosemHvL); chosemHvL = chosemHvL + 1;
        %%%%NEW%%%%
    else %condition of the picture chosen == 3 (goodmH)
        experiment_data(itrial,8) = rewardmHvH(chosemHvH); chosemHvH=chosemHvH+1; %%%%NEW%%%%
    end; %%%%NEW%%%%
    %    experiment_data(itrial,8)=data.rewardsize*(experiment_data(itrial,3+choicemade)>rand); %determine outcome %%%%NEW%%%%
    %Add variable in col 23 cumulatively summing rewards
     experiment_data(itrial,23) = data.currentscore + experiment_data(itrial,8);
    
    experiment_data(itrial,9) = rt;   %rt in s
    if choicemade==1
         cgpencol(0,0,0);
         cgrect(420,350,200,100);
         cgpencol(1,1,1)
         cgfont('Arial',data.FontSizeText)
         cgtext('Your Gems',420,365);
         cgtext(num2str(data.currentscore),420,330);
        cgpencol(1,.8,0)
        cgrect(-200,0,stimDimensions(1)+20,stimDimensions(2)+20);
        if experiment_data(itrial,3)>0; cgdrawsprite(2,200,0); end %Other stimulus does not disappear
        cgdrawsprite(1,-200,0);
    elseif choicemade==2
         cgpencol(0,0,0);
         cgrect(420,350,200,100);
         cgpencol(1,1,1)
         cgfont('Arial',data.FontSizeText)
         cgtext('Your Gems',420,365);
         cgtext(num2str(data.currentscore),420,330);
         cgpencol(1,.8,0)
        if experiment_data(itrial,2)>0; cgdrawsprite(1,-200,0); end %Other stimulus does not disappear
        cgrect(200,0,stimDimensions(1)+20,stimDimensions(2)+20);
        cgdrawsprite(2,200,0)
    end
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); cgtext('+',0,0);
    cgflip(data.background);
    wait(data.time.delay);
    
    data.currentscore = experiment_data(itrial,23);
    
    %%% display the outcome
    cgpencol(1,1,1); cgtext('+',0,0);
    img = data.rewardstim.picture; % Load reward stimuli (dimensions of which will be used for background of loss stimuli too)
    rewardstimDimensions = data.rewardstim.stimDimensions;
    cgloadbmp(3,img,rewardstimDimensions(1),rewardstimDimensions(2))
    cgfont('Arial',data.FontSizeTask);
    cgpencol(1,1,1); % White text
    if choicemade==1 %left pic
        cgpencol(1,.8,0)
        cgrect(-200,0,stimDimensions(1)+20,stimDimensions(2)+20);
        cgdrawsprite(1,-200,0)
        if experiment_data(itrial,3)>0; cgdrawsprite(2,200,0); end %Other stimulus does not disappear
            cgdrawsprite(3,-200,0)
            cgpencol(1,1,1)
            cgtext(num2str(experiment_data(itrial,8)), -200, 0);
            cgpencol(0,0,0);
            cgrect(420,350,200,100);
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Your Gems',420,365);
            cgtext(num2str(data.currentscore),420,330);
    elseif choicemade==2
            cgpencol(1,.8,0)
            cgrect(200,0,stimDimensions(1)+20,stimDimensions(2)+20);
            if experiment_data(itrial,2)>0; cgdrawsprite(1,-200,0); end %Other stimulus does not disappear
            cgdrawsprite(2,200,0)
            cgdrawsprite(3,200,0)
            cgpencol(1,1,1)
            cgtext(num2str(experiment_data(itrial,8)), 200, 0);
            cgpencol(0,0,0);
            cgrect(420,350,200,100);
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Your Gems',420,365);
            cgtext(num2str(data.currentscore),420,330);
    end
        
    
    cgflip(data.background); experiment_data(itrial,16) = (time - tstart)/1000; %outcome revealed %%%%NEW%%%%edited
    wait(data.time.outcometime);
    cgfont('Arial',data.FontSizeTask);
    cgpencol(1,1,1); cgtext('+', 0, 0);
    cgpencol(0,0,0);
    cgrect(420,350,200,100);
    cgpencol(1,1,1)
    cgfont('Arial',data.FontSizeText)
    cgtext('Your Gems',420,365);
    cgtext(num2str(data.currentscore),420,330);
    cgflip(data.background); experiment_data(itrial,17) = (time - tstart)/1000; %iti starts %%%%NEW%%%%edited
    wait(data.time.iti+round((data.time.maxiti-data.time.iti)*rand)); %jitter between min and max iti
    experiment_data(itrial,18) = (time - tstart)/1000; %trials end - could jitter %%%%NEW%%%%edited
end
wait(data.time.iti);

%% Run WoF
data.initials = 'test'

% Decide WoF draws and shuffle
wofdraws = [1 -1];
data.index_wofdraws = randperm( length(wofdraws) );
for i=1:length(wofdraws)
    wofdraws_run1(i) = wofdraws(data.index_wofdraws(i));
end
data.wofdraws=wofdraws_run1;

% Load wheel positions into buffer to save time later
cgpencol(1,1,1); cgfont('Arial',data.FontSizeText); cgtext(['Please take a rest while the next experiment loads'], 0, 50)
cgflip(data.background);
cgloadbmp(10,'Test_WoF_new.bmp'); w=500; h=500; %%have changed from 'WoF_01' --> Also need to change run_WoF_2
disp('Generating rotated sprites...')
tic
for ang=1:2:359
    buff=10+ang;
    cgmakesprite(buff, w,h, [0 0 0])
    cgsetsprite(buff)
    cgrotatesprite(10,0,0,-ang)
end
toc
cgsetsprite(0)
outcome = run_WoF_2(data, 1)