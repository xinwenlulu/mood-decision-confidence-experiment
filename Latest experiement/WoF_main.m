% %% This needs to go in top level meta-script! Subject ID has to be
% the same on both calls (and no 'subject 1-2', 'subject 1-2', etc, needs to be integer) otherwise
% it will crash / they won't do both task variants
% 

% LMblocks = 6; % Because 2x3
% LMconditions = [20 60]; % Reward probabilities
% LMpictures = {'bull.bmp';'chick.bmp';'crab.bmp';'fox.bmp';'hedgehog.bmp';'hippopotamus.bmp';'koala.bmp';'lemur.bmp';'pig.bmp';'tiger.bmp';'whale.bmp';'zebra.bmp'};
% LMpicprob = zeros(LMblocks,4); %stimuli used and associated probabilities
% LMindex_side = randperm(length(LMpictures)); % Shuffle stimuli
% for n=1:LMblocks,
%    LMindex_prob = randperm(length(LMconditions)); % Shuffle reward contingencies
%    LMpicprob(n,:) = [LMindex_side(2*n-1:2*n) LMconditions(LMindex_prob)/100];
%      LMpicprob(n,:) = [LMindex_side(2*n-1:2*n) LMconditions/100];
% end;



function [data] = WoF_main(subjectid, subjinitials, picprob, currentrun)

%% GENERAL STARTUP
%clear
% subjectid
%subjectid = lower(input('Subject number? ','s'));

%addpath(genpath('D:\Data\Cogent2000v1.33'))
addpath(genpath('C:\Program Files\MATLAB\toolbox'))
rand('seed',sum(100*clock));



data                  = struct;
data.nblocks = 6; % Because 2x3
conditions = [20 60]; % Reward probabilities
data.pictures = {'bull.bmp';'chick.bmp';'crab.bmp';'fox.bmp';'hedgehog.bmp';'hippopotamus.bmp';'koala.bmp';'lemur.bmp';'pig.bmp';'tiger.bmp';'whale.bmp';'zebra.bmp'};

%%% Comment this out if picprob passed into function
%%% 
% picprob = zeros(data.nblocks,4); %stimuli used and associated probabilities
% index_side = randperm(length(data.pictures)); % Shuffle stimuli
% for n=1:data.nblocks,
% %    index_prob = randperm(length(conditions)); % Shuffle reward contingencies
% %    picprob(n,:) = [index_side(2*n-1:2*n) conditions(index_prob)/100];
%      picprob(n,:) = [index_side(2*n-1:2*n) conditions/100];
% end;
%%% 
%%% Comment this out if picprob passed into function

%%% Depending on if doing one or two runs of experiment
%data.picprob = picprob; %stimuli and probabilities
picprob_run1=picprob(1:data.nblocks/2,:);
picprob_run2=picprob((data.nblocks/2)+1:end,:);
%%%

% Decide WoF draws and shuffle
wofdraws = [1 -1];
data.index_wofdraws = randperm( length(wofdraws) );
for i=1:length(wofdraws)
    wofdraws_run1(i) = wofdraws(data.index_wofdraws(i));
end;
picprob_run1(1,5:6) = [0, 0];
picprob_run1(2,5:6) = [wofdraws_run1(1), wofdraws_run1(1)];
picprob_run1(3,5:6) = [wofdraws_run1(2), wofdraws_run1(2)]


%%% Second run
data.picprob = picprob_run2;
%%% Manually reshuffle WoF draws 
data.index_wofdraws = randperm( length(wofdraws) );
for i=1:length(wofdraws)
    wofdraws_run2(i) = wofdraws(data.index_wofdraws(i));
end;
picprob_run2(1,5:6) = [0, 0];
picprob_run2(2,5:6) = [wofdraws_run2(1), wofdraws_run2(1)];
picprob_run2(3,5:6) = [wofdraws_run2(2), wofdraws_run2(2)]



%%% initialise block structure


if strcmp(subjectid,'0') || subjectid == 0
    data.initials = 'debugging';
    data.dob = '01/01/01';
    debugging = 1;
else
    %data.initials = lower(input('Subject initials? ','s'));
    data.initials = subjinitials;
    %data.dob  = input('Date of birth (dd/mm/yyyy)? ','s');

end

%data.id               = ['Smartphone_v2_',num2str(subjectid),'_01'];
data.starttimenum     = now;
data.starttime        = datestr(data.starttimenum,'HH:MM:SS');
data.date             = datestr(date,'yyyy-mm-dd');
%data.dir              = [data.id '_' data.date];
%data.filename         = sprintf('%s_sl7mri_%s_%s_001,data.id,data.starttime(1:2),data.starttime(4:5)');
data.netdir           = '\\abba\RutledgeLab\rawdata\LM\Smartphone_Data';
%eval(sprintf('save %s data',fullfile(data.dir,data.filename))); %save data locally

%make local data directory
%if ~exist(data.dir),
    %mkdir(data.dir);
%end;
% if ~isempty(dir(fullfile(data.dir,'*_MoodSample_*.mat'))), %check if existing data file
%     fprintf(1,'Date file for this task for %s already exists.\n',data.id);
%     overwrite         = input('Continue and make a new data file (yes/no)? ','s');
%     if length(overwrite) && strmatch(overwrite,'yes'), %cant just press enter has to be y or yes
%     else
%         return;
%     end;
% end;

% if exist(data.netdir), %can access network directory
%     if ~exist(fullfile(data.netdir,data.dir)),
%         %mkdir(fullfile(data.netdir,data.dir));
%     end;
%     eval(sprintf('save %s data',fullfile(data.netdir,data.dir,data.filename))); %save data
%     fprintf(1,'data saved on network\n');
% else
%     fprintf(1,'data NOT saved on network\n');
% end;



%%% Cogent configuration
% screenMode              = 0;                 % 0 for small window, 1 for full screen, 2 for second screen if attached
% screenRes               = 3;                 % 1=640x480; 2=800x600; 3=1024x768
foreground              = [0.5 0.5 0.5];           % foreground colour 111 is white
 data.background         = [0.5 0.5 0.5];           % background colour
fontName                = 'Arial';           % font parameters
fontSize                = 32;
nbits                   = 0;                 % 0 selects the maximum possible bits per pixel

% config_display(screenMode, screenRes, data.background, foreground, fontName, fontSize, 5, nbits);   % open graphics window
% config_sound;
% config_keyboard;
% start_cogent;

% Sounds for WoF
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
data.npracticetrials = 10; % divisible by 2
data.ntrials = 22; % % divisible by 2
%data.nblocks = 4;

disp(['nblocks=,', num2str(data.nblocks), ' ntrials=', num2str(data.ntrials)])

data.pictures = {'bull.bmp';'chick.bmp';'crab.bmp';'fox.bmp';'hedgehog.bmp';'hippopotamus.bmp';'koala.bmp';'lemur.bmp';'pig.bmp';'tiger.bmp';'whale.bmp';'zebra.bmp'};
data.stimDimensions = [300, 300];
data.rewardstim.picture = 'gemstone.bmp';
data.rewardstim.stimDimensions = [98, 93];
data.norewardstim.picture = 'Silhouettegemstone.bmp';
data.norewardstim.stimDimensions = [98, 93];

% Reward
data.endowment = 100; %100 gems to begin with
data.rewardsize = 10;
data.currentscore = data.endowment;


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



%% Not needed for group battery

% WELCOME AND OVERVIEW
% cgflip(data.background);
% cgpencol(1,1,1);
% cgtext('Welcome to this experiment.',0,300);
% cgtext('Before you begin, please provide the following ratings when prompted.',0,250)
% cgtext('Click when you are ready to proceed.',0,100); 
% cgflip(data.background);
% wait(data.timeout_instruction );
% [~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;

% % Lifetime happiness rating
%[data.happybaseline(1,1), data.happybaseline(1,2), data.happybaseline(1,3)] = happy_rating('lifehappy', data.time.ratewait+3000, data);

% %Pause before starting questionnaires
% cgflip(data.background);
% cgpencol(1,1,1);
% cgfont('Arial',data.FontSizeText);
% cgtext('You will now be presented with several questionnaires.',0,250);
% cgtext('Please follow the instructions on each page.',0,200);
% cgtext('Click when you are ready to begin the questionnaires.',0,100);
% cgflip(data.background);
% wait(data.timeout_instruction );
% [~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;

% % QUESTIONNAIRES
%  data = STAI(data,0);
%  data = PHQ(data,0);
%  data = BDI(data,0);
%  data = HPS(data,0);
%  data = RR(data,0);
 
%  
% save data
% eval(sprintf('save %s data',fullfile(data.dir,data.filename)));


%% MAIN TASK AND INSTRUCTIONS
cgflip(data.background);
cgpencol(1,1,1); cgfont('Arial',data.FontSizeText);
cgtext('Welcome to Safari!',0,300);
cgtext('Use the left and right mouse buttons to collect gems from the animals',0,150);
cgtext('You will earn more gems if you figure out which animals have more',0,50);
cgtext(['You start with ', num2str(data.endowment), ' gems. Collect as many as you can!'],0,-50);
cgtext('Click to continue',0,-200);
cgflip(data.background);
wait(data.time.instructiontime);
[~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;

%Life happiness rating previously went here
%wait(1500);
%[data.happybaseline(2,1), data.happybaseline(2,2), data.happybaseline(2,3)] = happy_rating('nowhappy', data.time.ratewait+3000, data);


%%% Call tasks
%%% Task version: 0 = [20 60], 1 = [25 75]

if mod(subjectid, 2) == 0 % Subj number even
    taskorder = [0 1];
    taskver = {'20_60'; '25_75'};
else % Subject number odd
    taskorder = [1 0];
    taskver = {'25_75'; '20_60'};
end 
    
%%% Check which run we are in
if currentrun == 1
    data.picprob = picprob_run1;
    data.wofdraws = wofdraws_run1;
elseif currentrun ==2
    data.picprob = picprob_run2;
    data.wofdraws = wofdraws_run2;
end

data.savefile               = ['Smartphone_v2_subj',num2str(subjectid),'_task', taskver{currentrun}, '_date', data.date, '_time', data.starttime(1:2),'-', data.starttime(4:5)]; % Change save file
disp(['Saving to: ', data.savefile]); save(data.savefile, 'data')

data.maintask{1}                = MoodSampling_task_LM(data,1,0, taskorder(currentrun)); save(data.savefile, 'data');    
data.currentscore               = data.maintask{1}(data.ntrials,23);
outcome = run_WoF_2(data, 1); data.currentscore = data.currentscore + outcome;

data.maintask{2} = MoodSampling_task_LM(data,2,data.wofdraws(1), taskorder(currentrun)); save(data.savefile, 'data');
data.currentscore               = data.maintask{2}(data.ntrials,23);
data.test{1} = testblock(data, 2); save(data.savefile, 'data');
data.currentscore = data.currentscore + 50; % %For now, just estimate gems collected in test block (average reward rate is .4)
outcome = run_WoF_2(data, 2); data.currentscore = data.currentscore + outcome;

data.maintask{3} = MoodSampling_task_LM(data,3,data.wofdraws(2), taskorder(currentrun)); save(data.savefile, 'data');
data.test{2} = testblock(data, 3); save(data.savefile, 'data');
data.currentscore = data.currentscore + 80; % %For now, just estimate gems collected in test block (average reward rate is .4)

%%% End
%%% CURRENTLY FEEDBACK DOES NOT TAKE ACCOUNT OF TEST BLOCK

cgpencol(1,1,1); cgfont('Arial',data.FontSizeText);
cgtext('Congratulations!',0,300);
cgtext('This task is over',0,250);
%cgtext('You can call the experimenter.',0,200);
cgtext(['You found: ', num2str(data.currentscore), ' gems'], 0, 100)
cgtext(['Well done!'], 0, 50)
cgtext('Click to continue',0,-50);
cgflip(data.background);
[~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;

[data.happyendtask(2,1), data.happyendtask(2,2), data.happyendtask(2,3)] = happy_rating('happy', data.time.ratewait+3000, data);
[~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;

%stop_cogent