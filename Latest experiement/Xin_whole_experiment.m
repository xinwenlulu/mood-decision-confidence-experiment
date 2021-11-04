
clear
% subjectid
subjectid = lower(input('Subject number? ','s'));
initials = lower(input('Subject initials? ','s'));
dob  = input('Date of birth (dd/mm/yyyy)? ','s');
gender = lower(input('Gender: ', 's'));

%% Call entire experiment
addpath 'Cogent2000v1.33\Toolbox'
addpath(genpath('C:\Program Files\MATLAB\toolbox'))

blocks = 2; % before vs after WoF;
conditions = [1 2 3]; % Reward distributions 1 = bad; 2 = goodvL; 3 = goodvH
%pictures = {'bull.bmp';'chick.bmp';'crab.bmp';'fox.bmp';'hedgehog.bmp';'hippopotamus.bmp';'koala.bmp';'lemur.bmp';'pig.bmp';'tiger.bmp';'whale.bmp';'zebra.bmp'};
picprob = zeros(blocks,6); %stimuli used and associated probabilities
index_side = randperm(6); % Shuffle stimuli
for n=1:blocks,
   index_prob = randperm(length(conditions)); % Shuffle reward contingencies
   picprob(n,:) = [index_side(3*n-2:3*n) conditions];
end;

screen              = 0;                 % 0 for small window, 1 for full screen, 2 for second screen if attached

if mod(str2num(subjectid), 2) == 0
    taskorder = [0 1]; %0=Xin 1= Sascha
else
   taskorder = [1 0];
end 

wof = [1,-1];
score = [0,5];

 data = Xin_WoF_main(screen,subjectid,initials,dob,gender,picprob,-1,score(1),0)


if taskorder(1) == 0 % taskorder = [0 1]
    data = Xin_WoF_main(screen, subjectid,initials,dob,gender, picprob,wof(1),score(1),0) %Xin's experiment

    pause
    starttime        = datestr(now,'HH:MM:SS');
    date             = datestr(date,'yyyy-mm-dd');
    
    questionnaires = struct;
    questionnaires = call_questionnaires(questionnaires,screen,0);
    questionnaires.savefile = ['Questionnaires_subj',num2str(subjectid),'_date', date, '_time',starttime(1:2),'-',starttime(4:5)]; % Change save file
    disp(['Saving to: ', questionnaires.savefile]); save(questionnaires.savefile, 'questionnaires');

    pause

    WoF_main_Sascha(screen,subjectid, initials, dob, gender, wof(2), score(2),data.currentearning) %Sascha's experiment
else % taskorder = [1 0]
    data = WoF_main_Sascha(screen,subjectid, initials, dob, gender, wof(1),score(1),0) %Sascha's experiment
   
    pause
    starttime        = datestr(now,'HH:MM:SS');
    date             = datestr(date,'yyyy-mm-dd');
    
    questionnaires = struct;
    questionnaires = call_questionnaires(questionnaires,screen,0);
    questionnaires.savefile = ['Questionnaires_subj',num2str(subjectid),'_date', date, '_time',starttime(1:2),'-',starttime(4:5)]; % Change save file
    disp(['Saving to: ', questionnaires.savefile]); save(questionnaires.savefile, 'questionnaires');
    pause   
    
    Xin_WoF_main(screen, subjectid,initials,dob,gender, picprob, wof(2),score(2),data.currentearning) %Xin's experiment
end

