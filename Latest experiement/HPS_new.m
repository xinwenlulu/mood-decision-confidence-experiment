% this is an Patient Health Questionnaire. data is entered using numeric keys.
function player_struct= HPS_new(player_struct,debugging)

instructionwait = 4000; %minimum wait in ms after instructions screen
player_struct.lagtime = 500; %lag time in ms so subject can see number typed

if debugging == 1;
    instructionwait = 0;
    player_struct.lagtime = 100;
end

keys = getkeymap;
player_struct.abortkey = [keys.F12]; %no quitting allowed in LOTR actually
player_struct.controlkey = [keys.Return];
player_struct.respkeys = [keys.K1 keys.K2 keys.K3 keys.K4 keys.K5]; %for buttons 1-4
player_struct.subjectkey = [keys.Space];


%%----------------------------------------------------------------------
%%run the PHQ
%%----------------------------------------------------------------------
% Instructions
cgflip(0.5,0.5,0.5);
cgfont('Arial',30);
cgtext('This questionnaire consists of 12 statements',0,200);
%cgtext('Each of the following statements describe people’s behaviour',0,200);
cgtext('Please indicate how accurately each statement describes you.',0,150);

cgtext('Press the space bar when ready.',0,-300);


cgflip(0.5,0.5,0.5);
wait(instructionwait); %wait 3s at least
waitkeydown(inf,player_struct.subjectkey); %wait for space bar to be pressed

% HPS - 20 statements

HPS_items = {'I frequently get into moods where I feel very speeded up and irritable.';
    'I think that my moods don’t change more than most people’s do.';
    'I have often felt happy and irritable at the same time.';
    'I can slow myself down when I want to.';
    'I am a person whose moods go up and down easily.';
    'I frequently find that my thoughts are racing.';
    'I am usually in an average sort of mood, not too high and not too low.';
    'I am often so restless that it is impossible for me to sit still.';
    'I get so happy or energetic that I am almost giddy.';
    'I feel emotions with extreme intensity.';
    'I am considered to be kind of eccentric.';
    'When I feel very excited and happy, I almost always know the reason.';
    };

player_struct.HPS.HPS_items = HPS_items;
HPS_raw = zeros(length(HPS_items),2); %Creating matrix for answers


for n = 1:length(HPS_items) % loop over all 26 MASQ items  
    
    cgflip(0.5,0.5,0.5);
    cgfont('Arial',30);
    cgtext('How accurately does the following statement describe you?',0,300);
    cgfont('Arial',30); setforecolour(0.5,0,0.5);
  %  if n<6
        cgtext(HPS_items{n},0,150); 
%     elseif n==6 
%         cgtext(HPS_items{n},0,180);
%         cgtext('failure or have let yourself or your family down',0,150);        
%     elseif n==7 
%         cgtext(HPS_items{n},0,180);
%         cgtext('reading the newspaper or watching television',0,150);        
%     elseif n==8 
%         cgtext(HPS_items{n},0,200);
%         cgtext('people could have noticed. Or the opposite -',0,170);
%         cgtext('being so fidgety or restless that you have',0,140); 
%         cgtext('been moving around a lot more than usual',0,110);           
%     elseif n==9 
%         cgtext(HPS_items{n},0,180);
%         cgtext('or of hurting yourself in some way',0,150);        
%     end
        
        
    cgfont('Arial',24);
    cgtext('1 = Very inaccurate',0,80);
    cgtext('2 = Moderately inaccurate',0,50);
    cgtext('3 = Neither accurate nor inaccurate',0,20);
    cgtext('4 = Moderately accurate',0,-10);
    cgtext('5 = Very accurate',0,-40);
    
    cgtext('Your answer: ',-30,-110);
    cgflip(0.5,0.5,0.5);
    wait(player_struct.lagtime); %wait before response allowed
    clearkeys;
    [k, t, npress] = waitkeydown(inf,player_struct.respkeys);
    key = find(player_struct.respkeys == k(1));
    cgtext(num2str(key),70,-110);
    
    cgflip(0.5,0.5,0.5);
    cgfont('Arial',30);
    cgtext('How accurately does the following statement describe you?',0,300);
    
    cgfont('Arial',30); setforecolour(0.5,0,0.5);
    
    cgtext(HPS_items{n},0,150); 
       
    cgfont('Arial',24);
    cgtext('1 = Very inaccurate',0,80);
    cgtext('2 = Moderately inaccurate',0,50);
    cgtext('3 = Neither accurate nor inaccurate',0,20);
    cgtext('4 = Moderately accurate',0,-10);
    cgtext('5 = Very accurate',0,-40);
    
    
    cgtext('Your answer: ',-30,-110);
    
    cgflip(0.5,0.5,0.5);
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

cgfont('Arial',30);

end