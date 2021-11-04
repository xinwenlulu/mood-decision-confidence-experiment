% this is the Altman Self-Rating Mania Scale. data is entered using numeric keys.

function player_struct = Altman(screen,player_struct,debugging)
%% General setup and start Cogent

foreground              = [0.5 0.5 0.5];           % foreground colour 111 is white
background              = [0.5 0.5 0.5];           % background colour
fontName                = 'Arial';           % font parameters
fontSize                = 32;
waitTime                = 4000;

%% Instructions 

 cgflip(background);
 cgpencol(1,1,1);
 cgfont('Arial',32);
 cgtext('You will now be presented with the second questionnaire.',0,250);
 cgtext('There are 5 groups of statements in this questionnaire.',0,150);
 cgtext('Read each group of statements carefully.',0,100);
 cgtext('You should choose the statement in each group that best',0,0);
 cgtext('describes the way you have been feeling for the past week.',0,-50);
 cgflip(background);
 wait(waitTime);
 cgflip(background);
 cgpencol(1,1,1);
 cgfont('Arial',32);
 cgtext('You will now be presented with a brief questionnaire.',0,250);
 cgtext('There are 5 groups of statements in this questionnaire.',0,150);
 cgtext('Read each group of statements carefully.',0,100);
 cgtext('You should choose the statement in each group that best',0,0);
 cgtext('describes the way you have been feeling for the past week.',0,-50);
 cgtext('Click when you are ready to begin the questionnaire.',0,-250);
 cgflip(background);
 [~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;
 
 cgflip(background);
 cgpencol(1,1,1);
 cgfont('Arial',32);
 cgtext('Please note:',0,150);
 cgtext('The word OCCASIONALLY when used here means once or twice,',0,50);
 cgtext('OFTEN means several times or more, and',0,0);
 cgtext('FREQUENTLY means most of the time.',0,-50);
 cgflip(background);
 wait(waitTime);
 cgflip(background);
 cgpencol(1,1,1);
 cgfont('Arial',32);
 cgtext('Please note:',0,150);
 cgtext('The word OCCASIONALLY when used here means once or twice,',0,50);
 cgtext('OFTEN means several times or more, and',0,0);
 cgtext('FREQUENTLY means most of the time.',0,-50);
 cgtext('Click when you are ready to begin the questionnaire.',0,-250);
 cgflip(background);
 [~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Its a set up, Its a set up, Its a set up %%

instructionwait = 1000; %minimum wait in ms after instructions screen
player_struct.lagtime = 500; %lag time in ms so subject can see number typed

if debugging == 1;
    instructionwait = 0;
    player_struct.lagtime = 100;
end

keys = getkeymap;
player_struct.abortkey = [keys.F12]; %no quitting allowed in LOTR actually
player_struct.controlkey = [keys.Return];
player_struct.respkeys = [keys.K1 keys.K2 keys.K3 keys.K4 keys.K5]; %for buttons 1-5
player_struct.subjectkey = [keys.Space];

%%----------------------------------------------------------------------
%%run the Altman scale
%%----------------------------------------------------------------------

% HPS - 20 statements

% HPS_items = {'I get into moods where I feel very speeded up and irritable.'; 
%     'I think that my moods don’t change more than most people’s do.';
%     'I tend to feel happy and irritable at the same time.';
%     'I can slow myself down when I want to.';
%     'I am a person whose moods go up and down easily.';
%     'I find that my thoughts are racing.';
%     'I am usually in an average sort of mood, not too high and not too low.';
%     'I am often so restless that it is impossible for me to sit still.';
%     'I get so happy or energetic that I am almost giddy.';
%     'I feel emotions with extreme intensity.';
%     'I am considered to be kind of eccentric.';
%     'I understand the reasons when I feel very excited or happy.';
%     };
% 
% player_struct.HPS.HPS_items = HPS_items;

Altman_items = {'Positive Mood';
    'Self-Confidence';
    'Sleep Patterns';
    'Speech';
    'Activity Levels';
    };
player_struct.Altman.Altman_items = Altman_items;
Altman_raw = zeros(length(Altman_items),2); %Creating matrix for answers

%Item 1 - Positive Mood

    clearpict(3);
    settextstyle('Arial',30);
    preparestring('Please choose the statement that best describes the way',3,0,250);
    preparestring('you have been feeling for the past week.',3,0,200);
    settextstyle('Arial',36); setforecolour(0.5,0,0.5); setforecolour(1,1,1);
    settextstyle('Arial',30);
    
    preparestring('1 = I do not feel happier or more cheerful than usual.',3,0,80);
    preparestring('2 = I occasionally feel happier or more cheerful than usual.',3,0,50);
    preparestring('3 = I often feel happier or more cheerful than usual.',3,0,20);
    preparestring('4 = I feel happier or more cheerful than usual most of the time.',3,0,-10);
    preparestring('5 = I feel happier or more cheerful than usual all of the time.',3,0,-40); %start
    
    preparestring('Your answer: ',3,-30,-110);
    drawpict(3);
    wait(player_struct.lagtime); %wait before response allowed
    clearkeys;
    [k, t, npress] = waitkeydown(inf,player_struct.respkeys);
    key = find(player_struct.respkeys == k(1));
    preparestring(num2str(key),3,70,-110);
    drawpict(3);
    wait(player_struct.lagtime); %wait so that subjects can see their answer
    Altman_raw(1,1) = key-1; % save answers in matrix
    
%Item 2 - Self-Confidence
    
    clearpict(3);
    settextstyle('Arial',30);
    preparestring('Please choose the statement that best describes the way',3,0,250);
    preparestring('you have been feeling for the past week.',3,0,200);
    settextstyle('Arial',36); setforecolour(0.5,0,0.5); setforecolour(1,1,1);
    settextstyle('Arial',30);
    
    preparestring('1 = I do not feel more self-confident than usual.',3,0,80);
    preparestring('2 = I occasionally feel more self-confident than usual.',3,0,50);
    preparestring('3 = I often feel more self-confident than usual.',3,0,20);
    preparestring('4 = I feel more self-confident than usual.',3,0,-10);
    preparestring('5 = I feel extremely self-confident all of the time.',3,0,-40); %start
    
    preparestring('Your answer: ',3,-30,-110);
    drawpict(3);
    wait(player_struct.lagtime); %wait before response allowed
    clearkeys;
    [k, t, npress] = waitkeydown(inf,player_struct.respkeys);
    key = find(player_struct.respkeys == k(1));
    preparestring(num2str(key),3,70,-110);
    drawpict(3);
    wait(player_struct.lagtime); %wait so that subjects can see their answer
    Altman_raw(2,1) = key-1; % save answers in matrix
    
%Item 3 - Sleep Patterns
    
    clearpict(3);
    settextstyle('Arial',30);
    preparestring('Please choose the statement that best describes the way',3,0,250);
    preparestring('you have been feeling for the past week.',3,0,200);
    settextstyle('Arial',36); setforecolour(0.5,0,0.5); setforecolour(1,1,1);
    settextstyle('Arial',30);
    
    preparestring('1 = I do not need less sleep than usual.',3,0,80);
    preparestring('2 = I occasionally need less sleep than usual.',3,0,50);
    preparestring('3 = I often need less sleep than usual.',3,0,20);
    preparestring('4 = I frequently need less sleep than usual.',3,0,-10);
    preparestring('5 = I can go all day and night without any sleep',3,0,-40);
    preparestring('and still not feel tired.',3,0,-70); %start
    
    preparestring('Your answer: ',3,-30,-110);
    drawpict(3);
    wait(player_struct.lagtime); %wait before response allowed
    clearkeys;
    [k, t, npress] = waitkeydown(inf,player_struct.respkeys);
    key = find(player_struct.respkeys == k(1));
    preparestring(num2str(key),3,70,-110);
    drawpict(3);
    wait(player_struct.lagtime); %wait so that subjects can see their answer
    Altman_raw(3,1) = key-1; % save answers in matrix

%Item 4 - Speech
    
    clearpict(3);
    settextstyle('Arial',30);
    preparestring('Please choose the statement that best describes the way',3,0,250);
    preparestring('you have been feeling for the past week.',3,0,200);
    settextstyle('Arial',36); setforecolour(0.5,0,0.5); setforecolour(1,1,1);
    settextstyle('Arial',30);
    
    preparestring('1 = I do not talk more than usual.',3,0,80);
    preparestring('2 = I occasionally talk more than usual.',3,0,50);
    preparestring('3 = I often talk more than usual.',3,0,20);
    preparestring('4 = I frequently talk more than usual.',3,0,-10);
    preparestring('5 = I talk constantly and cannot be interrupted.',3,0,-40); %start
    
    preparestring('Your answer: ',3,-30,-110);
    drawpict(3);
    wait(player_struct.lagtime); %wait before response allowed
    clearkeys;
    [k, t, npress] = waitkeydown(inf,player_struct.respkeys);
    key = find(player_struct.respkeys == k(1));
    preparestring(num2str(key),3,70,-110);
    drawpict(3);
    wait(player_struct.lagtime); %wait so that subjects can see their answer
    Altman_raw(4,1) = key-1; % save answers in matrix

%Item 5 - Activity Level
    
    clearpict(3);
    settextstyle('Arial',30);
    preparestring('Please choose the statement that best describes the way',3,0,250);
    preparestring('you have been feeling for the past week.',3,0,200);
    settextstyle('Arial',36); setforecolour(0.5,0,0.5); setforecolour(1,1,1);
    settextstyle('Arial',30);
    
    preparestring('1 = I have not been more active (either socially, sexually,',3,0,80);
    preparestring('at work, home or school) than usual.',3,0,50);
    preparestring('2 = I have occasionally been more active than usual.',3,0,20);
    preparestring('3 = I have often been more active than usual.',3,0,-10);
    preparestring('4 = I have frequently been more active than usual.',3,0,-40);
    preparestring('5 = I am constantly active or on the go all the time.',3,0,-70); %start
    
    preparestring('Your answer: ',3,-30,-110);
    drawpict(3);
    wait(player_struct.lagtime); %wait before response allowed
    clearkeys;
    [k, t, npress] = waitkeydown(inf,player_struct.respkeys);
    key = find(player_struct.respkeys == k(1));
    preparestring(num2str(key),3,70,-110);
    drawpict(3);
    wait(player_struct.lagtime); %wait so that subjects can see their answer
    Altman_raw(5,1) = key-1; % save answers in matrix
    
 cgflip(background);
 cgpencol(1,1,1);
 cgfont('Arial',32);
 cgtext('Thank you for completing this questionnaire.',0,100);
 cgtext('Press any mouse button to continue.',0,-150);
 cgflip(background);
 [~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;
    
% reverse coding for items

Altman_score = Altman_raw(1:5,1);
Altman_score = Altman_score';

Altmanpos = [1 2 3 4 5];
Altmanneg = [];

%Altman_score_mean = (mean([Altman_score(Altmanpos)
%5-Altman_score(Altmanneg)]) - 1) / 4; % not sure what this does. Why minus
%1 and divide by 4?
Altman_score = sum([Altman_score(Altmanpos) 5-Altman_score(Altmanneg)]);

player_struct.Altman.Altman_raw = Altman_raw(:,1);
%player_struct.Altman.Altman_score_mean = Altman_score_mean;
player_struct.Altman.Altman_score = Altman_score;
settextstyle('Arial',32);

end