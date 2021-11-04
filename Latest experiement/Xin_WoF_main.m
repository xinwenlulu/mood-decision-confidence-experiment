% %% This needs to go in top level meta-script! Subject ID has to be
% the same on both calls (and no 'subject 1-2', 'subject 1-2', etc, needs to be integer) otherwise
% it will crash / they won't do both task variants
% 

% blocks = 4; 
% conditions = [1 2 3]; %Reward distributions 1 = bad; 2 = goodvL; 3 = goodvH
% pictures = {'bull.bmp';'chick.bmp';'crab.bmp';'fox.bmp';'hedgehog.bmp';'hippopotamus.bmp';'koala.bmp';'lemur.bmp';'pig.bmp';'tiger.bmp';'whale.bmp';'zebra.bmp'};
% picprob = zeros(blocks,6); %stimuli used and associated probabilities
% index_side = randperm(length(pictures)); % Shuffle stimuli
% for n=1:blocks,
%index_prob = randperm(length(conditions)); % Shuffle reward contingencies
   %picprob(n,:) = [index_side(3*n-2:3*n) conditions(index_prob)];  
% end;



function [data] = Xin_WoF_main(screen,subjectid, subjinitials, dob, gender, picprob, wof,wofscore,previousearning)

%% GENERAL STARTUP

%addpath(genpath('D:\Data\Cogent2000v1.33'))
%addpath(genpath('C:\Program Files\MATLAB\toolbox'))
rand('seed',sum(100*clock));

data                  = struct;
data.nblocks = 2; % Because 2x3

if rand > 0.5
    data.pictures = {'2.bmp';'3.bmp';'4.bmp';'7.bmp';'9.bmp';'10.bmp'};
else
    data.pictures = {'3.bmp';'6.bmp';'7.bmp';'8.bmp';'9.bmp';'10.bmp'};
end

data.practice_pic = {'p1.bmp';'p2.bmp'};
data.confidence = {'left.bmp';'right.bmp'};
picprob_used=picprob(1:2,:); %use 2 blocks of pictures

% append wof to the array
picprob_used(1,7) = [0];
picprob_used(2,7) = [wof];

%%% initialise block structure

if  subjectid == '0'
    data.subjectid = subjectid;
    data.initials = 'debugging';
    data.dob = '01/01/01';
    debugging = 1;
else
    %data.initials = lower(input('Subject initials? ','s'));
    data.subjectid = subjectid;
    data.initials = subjinitials;
    data.gender = gender;
    data.dob = dob;
    %data.dob  = input('Date of birth (dd/mm/yyyy)? ','s');
end


%data.id               = ['Smartphone_v2_',num2str(subjectid),'_01'];
data.starttimenum     = now;
data.starttime        = datestr(data.starttimenum,'HH:MM:SS');
data.date             = datestr(date,'yyyy-mm-dd');


%%% Cogent configuration
screenMode              = screen;                 % 0 for small window, 1 for full screen, 2 for second screen if attached
screenRes               = 3;                 % 1=640x480; 2=800x600; 3=1024x768
foreground              = [0.5 0.5 0.5];           % foreground colour 111 is white
data.background         = [0.5 0.5 0.5];           % background colour
fontName                = 'Arial';           % font parameters
fontSize                = 32;
nbits                   = 0;                 % 0 selects the maximum possible bits per pixel

 config_display(screenMode, screenRes, data.background, foreground, fontName, fontSize, 5, nbits);   % open graphics window
 config_sound;
 config_keyboard;
 start_cogent;

% Sounds for WoF
loadsound('Cheer2.wav', 15);
loadsound('Ding3.wav',16);
loadsound('down.wav',17);
loadsound('Boo.wav',18);

data.sidelocation             = [-200 200];
data.FontSizeText             = 32;
data.FontSizeTask             = 80;

% task timing
data.iti                      =  1000;
data.iti_std                  =  300;
data.iti_instruction          =  800;
data.timeout                  = 1000;
data.timeout_instruction      = 1200;
data.optiontime               =    0;
data.confirmationtime         =  600;
data.outcometime              =  1500;
data.time.outcome_delivery    =  600;
data.time.outcome_delivery_practice    = 0;


data.time.instructiontime     = 1500;
data.time.delay               = 800; %delay until outcome revealed
data.time.outcometime         = 1500;
data.time.iti                 = 600; %min ITI time
data.time.maxiti              = 1200; %max ITI time
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
data.npracticetrials = 10; 
data.ntrials = 44; % % divisible by 8 + 4 (2learning paris * 4 section) one mHcomparison per pair
%data.nblocks = 4;

disp(['nblocks=,', num2str(data.nblocks), ' ntrials=', num2str(data.ntrials)])

data.stimDimensions = [300, 300];
data.rewardstim.picture = 'gemstone.bmp';
data.rewardstim.stimDimensions = [100, 100];
data.norewardstim.picture = 'Silhouettegemstone.bmp';
data.norewardstim.stimDimensions = [98, 93];
data.winwofstim.picture = 'goldcoin.bmp';
data.winwofstim.stimDimensions = [98, 93];
data.losewofstim.picture = 'Greycoin.bmp';
data.losewofstim.stimDimensions = [98, 93];


% Reward
data.currentscore = 0;
data.wofscore = wofscore;
data.currentearning = 0;

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


%% MAIN TASK AND INSTRUCTIONS
cgpencol(1,1,1); cgfont('Arial',data.FontSizeText);
if data.wofscore == 0
    cgtext('Welcome to this experiment',0,250);
    cgtext('In this experiment, you will play two games.',0,150);
    cgtext('In each game, your task is to earn as many gems as possible.',0,50);
    cgtext('At the end of the experiment, your gems will be ',0,-50);
    cgtext('converted to £, and added to your payout.',0,-150);
    cgflip(data.background);
    wait(data.time.instructiontime);
    cgtext('Welcome to this experiment',0,250);
    cgtext('In this experiment, you will play two games.',0,150);
    cgtext('In each game, your task is to earn as many gems as possible.',0,50);
    cgtext('At the end of the experiment, your gems will be ',0,-50);
    cgtext('converted to £, and added to your payout.',0,-150);
    cgtext('Click to continue.',0,-300);
    cgflip(data.background);
    [~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;
end    

% % Lifetime happiness rating
[data.happybaseline(1,1), data.happybaseline(1,2), data.happybaseline(1,3)] = happy_rating('lifehappy', data.time.ratewait+3000, data);

cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10);
cgtext('Welcome to Gem Garden!',0,250);
cgfont('Arial',data.FontSizeText+5);
cgtext('Move the mouse and use left-click to collect gems from the flowers.',0,100);
cgtext('The gems count towards your payout.',0,0);
cgtext('You will earn more gems if you figure out which flower has more.',0,-100);
cgflip(data.background);
wait(data.time.instructiontime);
cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10);
cgtext('Welcome to Gem Garden!',0,250);
cgfont('Arial',data.FontSizeText+5);
cgtext('Move the mouse and use left-click to collect gems from the flowers.',0,100);
cgtext('The gems count towards your payout.',0,0);
cgtext('You will earn more gems if you figure out which flower has more.',0,-100);
cgtext('Click to continue',0,-250);
cgflip(data.background);
[~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;

cgpencol(1,1,1); cgfont('Arial',data.FontSizeText);
cgtext('With each choice, you will also indicate your level of confidence',0,250);
cgtext('that the chosen flower will have more gems than the other one',0,150);
cgtext('Use the slider to indicate your confidence right at that moment',0,50);
cgtext('Please use the whole range of the scale',0,-50);
cgtext('Your confidence rating does not affect the number of gems you get',0,-200);
cgflip(data.background);
wait(data.time.instructiontime);
cgtext('With each choice, you will also indicate your level of confidence',0,250);
cgtext('that the chosen flower will have more gems than the other one',0,150);
cgtext('Use the slider to indicate your confidence right at that moment',0,50);
cgtext('Please use the whole range of the scale',0,-50);
cgtext('Your confidence rating does not affect the number of gems you get',0,-200);
cgtext('Click to continue',0,-300);
cgflip(data.background);
[~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;


cgtext('Are you ready for some practice?',0,200);
cgtext('Click to continue',0,-100);
cgflip(data.background);
wait(600);
[~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;

data.picprob = picprob_used;
data.wofdraws = wof;

data.savefile               = ['Mood&Confidence_subj',num2str(subjectid),'_date', data.date, '_time', data.starttime(1:2),'-', data.starttime(4:5)]; % Change save file
disp(['Saving to: ', data.savefile]); save(data.savefile, 'data')

%practice block here
[data.practice] = Xin_PracticeBlock(data); save(data.savefile, 'data');
data.currentscore = 0;

cgtext('Are you ready to enter the gem garden?',0,200);
cgtext('Click to continue',0,-100);
cgflip(data.background);
wait(600);
[~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;

%Block 1
 data.maintask{1}                = Xin_MoodSampling_task(data,1,0); save(data.savefile, 'data');    
 data.currentscore               = data.maintask{1}(data.ntrials,23);


cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+15);
cgtext(['You will see these flowers again!'],0,100);
cgflip(data.background);
wait(data.time.instructiontime);
cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+15);
cgtext(['You will see these flowers again!'],0,100);
cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+5);
cgtext('Click to continue',0,-300);
cgflip(data.background);
[~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;

% text before the WoF converting gems to pound amounts
data.divisor = 5000; %5000 gems per £1
data.currentearning = round(data.currentscore/data.divisor,2);
if data.wofscore == 0
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeText);
 
    cgtext(['You so far earned ', num2str(data.currentscore), ' gems.'],0,150);
    cgtext(['Every ', num2str(data.divisor), ' gems in this round were worth £1.00.'],0,50);
    cgtext(['That means you so far earned £',num2str(data.currentearning),' extra on top of your payout.'],0,-50);
else
    data.currentearning = data.currentearning + previousearning;
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeText);
    cgtext(['You so far earned ', num2str(data.currentscore), ' gems in this task.'],0,200);
    cgtext(['Every ', num2str(data.divisor), ' gems in this round were worth £1.00.'],0,100);
    cgtext(['That means together with your earnings from the previous experiment,'],0,0);
    cgtext(['you so far earned £',num2str(data.currentearning), ' in gems '],0,-100);  
    cgtext(['and £', num2str(data.wofscore),' bonus from the first wheel on top of your payout.'],0,-200);  
end
cgflip(data.background);
wait(data.time.instructiontime);
if data.wofscore == 0
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeText);
 
    cgtext(['You so far earned ', num2str(data.currentscore), ' gems.'],0,150);
    cgtext(['Every ', num2str(data.divisor), ' gems in this round were worth £1.00.'],0,50);
    cgtext(['That means you so far earned £',num2str(data.currentearning),' extra on top of your payout.'],0,-50);
else
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeText);
    cgtext(['You so far earned ', num2str(data.currentscore), ' gems in this task.'],0,200);
    cgtext(['Every ', num2str(data.divisor), ' gems in this round were worth £1.00.'],0,100);
    cgtext(['That means together with your earnings from the previous experiment,'],0,0);
    cgtext(['you so far earned £',num2str(data.currentearning), ' in gems '],0,-100);  
    cgtext(['and £', num2str(data.wofscore),' bonus from the first wheel on top of your payout.'],0,-200);  
end
cgtext('Click to continue',0,-300);
cgflip(data.background);
[~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;

outcome = run_WoF_2(data, 1); 
data.wofscore = data.wofscore+outcome;

cgflip(data.background);
cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+5);
cgtext(['Ready to meet more flowers?'],0,100);
cgtext('Click to continue',0,-200);
cgflip(data.background);
wait(data.time.instructiontime);
[~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;

%block 2
data.maintask{2} = Xin_MoodSampling_task(data,2,data.wofdraws); save(data.savefile, 'data');
data.currentscore               = data.maintask{2}(data.ntrials,23);

%test block
data.test{1} = Xin_testblock(data, 2); save(data.savefile, 'data');
data.currentscore = data.currentscore + 1000; % %For now, just estimate gems collected in test block (average reward rate is .4)
cgfont('Arial',data.FontSizeText+5); cgpencol(1,1,1); 
cgtext('Finally, we would like to ask you how consistent each flower is',0,200);
cgtext('(how much did the number of gems vary for each flower ?)',0,100);
cgtext('It did not vary at all means you always got the same number of gems',0,0);
cgtext('It varied a lot means you got a very different number of gems each time',0,-100);
cgflip(data.background);
wait(1000);
cgtext('Finally, we would like to ask you how consistent each flower is',0,200);
cgtext('(how much did the number of gems vary for each flower ?)',0,100);
cgtext('It did not vary at all means you always got the same number of gems',0,0);
cgtext('It varied a lot means you got a very different number of gems each time',0,-100);
cgtext('Click to continue',0,-250);
cgflip(data.background)
[~,~,~,mp]=cgmouse; mp=0; %no mouse press yet
while (mp ~=1) && (mp ~=4), %1 from leftmouse, 4 for rightmouse, 2 for middle mouse
    [~,~,~, mp] = cgmouse;
end;

data.variance{1}                = Xin_VarianceEstimate(data,1); save(data.savefile,'data');
data.variance{2}                = Xin_VarianceEstimate(data,2); save(data.savefile,'data');

%%% End
%%% CURRENTLY FEEDBACK DOES NOT TAKE ACCOUNT OF TEST BLOCK
data.currentearning = round(data.currentscore/data.divisor,2) + previousearning; %all gem money
if data.wofdraws == 1
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeText);
    %cgtext('Congratulations!',0,300);
    cgtext('This game is over',0,250);
    %cgtext('You can call the experimenter.',0,200);
    %cgtext(['You found: ', num2str(data.currentscore), ' gems'], 0, 150)
    cgtext(['The gems you earned have been banked'], 0, 150)
    %cgtext(['Well done! You earned £',num2str(data.currentearning),'altogether'], 0, 100)
    cgtext('Click to continue',0,-50);
else
    data.currentearning = data.currentearning + data.wofscore;
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeText);
    %cgtext('Congratulations!',0,300);
    cgtext('This game is over',0,250);
    %cgtext('You can call the experimenter.',0,200);
    %cgtext(['You found: ', num2str(data.currentscore), ' gems in this task'], 0, 150)
    if data.currentearning >= 0
        cgtext(['You earned £',num2str(data.currentearning),' in total from the 2 tasks'], 0, 100)
    else
        cgtext(['You earned - £',num2str(abs(data.currentearning)),' in total from the 2 tasks'], 0, 100)
    end 
    cgtext('Click to continue',0,-50);
end    
cgflip(data.background);
[~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;

save(data.savefile, 'data');
%[data.happyendtask(2,1), data.happyendtask(2,2), data.happyendtask(2,3)] = happy_rating('happy', data.time.ratewait+3000, data);
%[~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;

stop_cogent