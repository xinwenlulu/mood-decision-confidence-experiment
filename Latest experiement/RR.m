% This is a Reward responsiveness test. Data is entered using numeric keys.

function player_struct= RR(player_struct,debugging)

instructionwait = 4000; %minimum wait in ms after instructions screen
player_struct.lagtime = 500; %lag time in ms so subject can see number typed

if debugging == 1;
    instructionwait = 0;
    player_struct.lagtime = 100;
end

keys = getkeymap;
player_struct.abortkey = [keys.F12]; %no quitting allowed in LOTR actually
player_struct.controlkey = [keys.Return];
player_struct.respkeys = [keys.K1 keys.K2 keys.K3 keys.K4]; %for buttons 1-4
player_struct.subjectkey = [keys.Space];


%%----------------------------------------------------------------------
%%run the RR
%%----------------------------------------------------------------------
% Instructions
clearpict(3); %clear buffer 3
settextstyle('Arial',30);
preparestring('This questionnaire consists of 8 statements.',3,0,300);
preparestring('Each item of this questionnaire is a statement that',3,0,200);
preparestring('a person may either agree with or disagree with.',3,0,175);
preparestring('For each item, indicate how much you agree or disagree with what the item says.',3,0,100);
preparestring('Please respond to all items; do not leave any blank.',3,0,0);
preparestring('Choose only one response to each statement.',3,0,-50);
preparestring('Please be as accurate and honest as you can be.',3,0,-100);
%preparestring('Respond to each item as if it were the only item.',3,0,-150);
%preparestring('Do not worry about being "consistent" in your responses',3,0,-200);
preparestring('Press the space bar when ready.',3,0,-300);

drawpict(3);
wait(instructionwait); %wait 3s at least
waitkeydown(inf,player_struct.subjectkey); %wait for space bar to be pressed

% RR - 8 statements

RR_items = {'I am someone who goes all-out.';
    'If I discover something new I like, I usually continue doing it for a while.';
    'I would do anything to achieve my goals';
    'When I am successful at something, I continue doing it.';
    'When I go after something I use a "no holds barred" approach.';
    'When I see an opportunity for something I like, I get excited right away.';
    'When I am doing well at something, I love to keep at it.';
    'If I see a chance of something I want, I move on it right away.';
    };

player_struct.RR.RR_items = RR_items;
RR_raw = zeros(length(RR_items),2); %Creating matrix for answers

for n = 1:length(RR_items) % loop over all 8 RR items
    
    clearpict(3);
    settextstyle('Arial',30);
    preparestring('How accurately does the following statement describe you?',3,0,250); %remove this? Re-word maybe?
    settextstyle('Arial',36); setforecolour(0.5,0,0.5);preparestring(RR_items{n},3,0,150); setforecolour(1,1,1);
    settextstyle('Arial',30);
    
    preparestring('1 = Very true for me',3,0,80);
    preparestring('2 = Somewhat true for me',3,0,50);
    preparestring('3 = Somewhat false for me',3,0,20);
    preparestring('4 = Very false for me',3,0,-10); %start
    
    preparestring('Your answer: ',3,-30,-110);
    drawpict(3);
    wait(player_struct.lagtime); %wait before response allowed
    clearkeys;
    [k, t, npress] = waitkeydown(inf,player_struct.respkeys);
    key = find(player_struct.respkeys == k(1));
    preparestring(num2str(key),3,70,-110);
    drawpict(3);
    wait(player_struct.lagtime); %wait so that subjects can see their answer
    RR_raw(n,1) = key; % save answers in matrix
    
end

RR_score = RR_raw(1:8,1);
RR_score = sum(RR_score');

RR_score_mean = (RR_score)/8;

player_struct.RR.RR_raw = RR_raw(:,1);
player_struct.RR.RR_score_mean = RR_score_mean;
player_struct.RR.RR_score = RR_score;
settextstyle('Arial',32);

%{
% reverse coding for items

HPSpos = [1 3 5 6 8:11];
HPSneg = [2 4 7 12];


(mean([RR_score(HPSpos) 6-RR_score(HPSneg)]) - 1) / 4;
RR_score = sum([RR_score(HPSpos) 6-RR_score(HPSneg)]);            
         
%}

end