% this is Affective Instability test. data is entered using numeric keys.

function player_struct= HPS(screen, player_struct,debugging)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.FontSizeText             = 32;
data.timeout_instruction      = 1200;

foreground              = [0.5 0.5 0.5];           % foreground colour 111 is white
data.background         = [0.5 0.5 0.5];           % background colour
fontName                = 'Arial';           % font parameters
fontSize                = 32;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

instructionwait = 4000; %minimum wait in ms after instructions screen
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
%%run the HPS
%%----------------------------------------------------------------------
% Instructions
clearpict(3); %clear buffer 3
settextstyle('Arial',30);
preparestring('This questionnaire consists of 12 statements.',3,0,300);
preparestring('Each of the following statements describes people’s behaviour',3,0,200);
preparestring('Please indicate how accurately each statement describes you.',3,0,150);
%preparestring('Describe yourself as you generally are now,',3,0,0);
%preparestring('NOT as you wish to be in the future.',3,0,-50);
%preparestring('Describe yourself as you honestly see yourself,',3,0,-100);
%preparestring('in relation to other people you know of the same sex as',3,0,-150);
%preparestring('you are, and roughly your same age',3,0,-200);
preparestring('Press the space bar when ready.',3,0,-300);

drawpict(3);
%wait(instructionwait); %wait 3s at least
%waitkeydown(inf,player_struct.subjectkey); %wait for space bar to be pressed

% HPS - 20 statements

HPS_items = {'I get into moods where I feel very speeded up and irritable.'; 
    'I think that my moods don’t change more than most people’s do.';
    'I tend to feel happy and irritable at the same time.';
    'I can slow myself down when I want to.';
    'I am a person whose moods go up and down easily.';
    'I find that my thoughts are racing.';
    'I am usually in an average sort of mood, not too high and not too low.';
    'I am often so restless that it is impossible for me to sit still.';
    'I get so happy or energetic that I am almost giddy.';
    'I feel emotions with extreme intensity.';
    'I am considered to be kind of eccentric.';
    'I understand the reasons when I feel very excited or happy.';
    };

player_struct.HPS.HPS_items = HPS_items;
HPS_raw = zeros(length(HPS_items),2); %Creating matrix for answers

for n = 1:length(HPS_items) % loop over all 26 MASQ items
    
    clearpict(3);
    settextstyle('Arial',30);
    preparestring('How accurately does the following statement describe you?',3,0,250);
    settextstyle('Arial',36); setforecolour(0.5,0,0.5);preparestring(HPS_items{n},3,0,150); setforecolour(1,1,1);
    settextstyle('Arial',30);
    
    preparestring('1 = Very inaccurate',3,0,80);
    preparestring('2 = Moderately inaccurate',3,0,50);
    preparestring('3 = Neither accurate nor inaccurate',3,0,20);
    preparestring('4 = Moderately accurate',3,0,-10);
    preparestring('5 = Very accurate',3,0,-40); %start
    
    preparestring('Your answer: ',3,-30,-110);
    drawpict(3);
    wait(player_struct.lagtime); %wait before response allowed
    clearkeys;
    [k, t, npress] = waitkeydown(inf,player_struct.respkeys);
    key = find(player_struct.respkeys == k(1));
    preparestring(num2str(key),3,70,-110);
    drawpict(3);
    wait(player_struct.lagtime); %wait so that subjects can see their answer
    HPS_raw(n,1) = key; % save answers in matrix
    
end

% reverse coding for items

HPS_score = HPS_raw(1:12,1);
HPS_score = HPS_score';

HPSpos = [1 3 5 6 8:11];
HPSneg = [2 4 7 12];

HPS_score_mean = (mean([HPS_score(HPSpos) 6-HPS_score(HPSneg)]) - 1) / 4;
HPS_score = sum([HPS_score(HPSpos) 6-HPS_score(HPSneg)]);            
            
player_struct.HPS.HPS_raw = HPS_raw(:,1);
player_struct.HPS.HPS_score_mean = HPS_score_mean;
player_struct.HPS.HPS_score = HPS_score;
settextstyle('Arial',32);

end