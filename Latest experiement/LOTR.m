%% lossAversion_LOTR %%
% LOT-R optimism questionnaire for expt. One item per screen, response
% needs to be bound to num keys.
% Scheier et al. (1994) http://local.psy.miami.edu/faculty/ccarver/sclLOT-R.html

function player_struct = LOTR(player_struct,debugging)
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

% if debugging == 1;
%     instructionwait = 0;
%     player_struct.lagtime = 100;
% end

keys = getkeymap;
player_struct.abortkey = [keys.F12]; %no quitting allowed in LOTR actually
player_struct.controlkey = [keys.Return];
player_struct.respkeys = [keys.K1 keys.K2 keys.K3 keys.K4 keys.K5]; %for buttons 1-5
player_struct.subjectkey = [keys.Space];


% %% Instructions
% 
% background = [0.5 0.5 0.5];  
%  cgflip(background);
%  cgpencol(1,1,1);
%  cgfont('Arial',32);
%  cgtext('You will now be presented with a brief questionnaire.',0,250);
%  cgtext('This questionnaire consists of 10 questions.',0,100);
%  cgtext('about the past week, including today.',0,50);
%  cgflip(background);
%  wait(4000);
%  cgflip(background);
%  cgpencol(1,1,1);
%  cgfont('Arial',32);
%  cgtext('You will now be presented with a brief questionnaire.',0,250);
%  cgtext('This questionnaire consists of 10 questions.',0,100);
%  cgtext('about the past week, including today.',0,50);
%  cgtext('Click when you are ready to begin the questionnaire.',0,-150);
%  cgflip(background);
%  [~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Items from the LOTR %%

% LOTR Items
LOTRitems = {'In uncertain times, I usually expect the best.';
              'It is easy for me to relax.';
              'If something can go wrong for me, it will.';
              'I am always optimistic about my future.';
              'I enjoy my friends a lot.';
              'It is important for me to keep busy.';
              'I hardly ever expect things to go my way.';
              'I do not get upset too easily.';
              'I rarely count on good things happening to me.';
              'Overall, I expect more good things to happen to me than bad.'};
          
% Save these items into the player_struct
player_struct.LOTR.LOTRitems = LOTRitems; 
LOTRraw = zeros(length(LOTRitems),2);

%% Items start here %%

for item = 1:length(LOTRitems)

cgflip(data.background);
cgpencol(1,1,1); cgfont('Arial',30);
cgtext('To what extent have you experienced the following',0,300);
cgtext('over the past week, including today?',0,250);
cgtext(char(LOTRitems(item)),0,100);
cgtext('1 = Not at All',0,40);
cgtext('2 = Slightly',0,10);
cgtext('3 = Moderately',0,-20);
cgtext('4 = Very',0,-50);
cgtext('5 = Extremely',0,-80);
cgtext('Your answer: ',0,-150);
cgflip(data.background);
   
    wait(player_struct.lagtime);
    clearkeys;
    [k, t, npress] = waitkeydown(inf,player_struct.respkeys);
    key = find(player_struct.respkeys == k(1));
    
cgpencol(1,1,1); cgfont('Arial',30);
cgtext('To what extent have you experienced the following',0,300);
cgtext('over the past week, including today?',0,250);
cgtext('1 = Not at All',0,40);
cgtext('2 = Slightly',0,10);
cgtext('3 = Moderately',0,-20);
cgtext('4 = Very',0,-50);
cgtext('5 = Extremely',0,-80);
cgtext('Your answer: ',0,-150);
cgtext(num2str(key),0,-180);
cgflip(data.background);
    
    wait(player_struct.lagtime); %wait so that subjects can see their answer
    LOTRraw(item,1) = key; % save answers in matrix
   
end

%% Scored it %%

% Items 2, 5, 6 and 8 are 'fillers' so need to be omitted from actual score
fillerItems = [2 5 6 8];
realItems = [1 3 4 7 9 10];
player_struct.LOTR.raw = LOTRraw;
player_struct.LOTR.score = sum(LOTRraw(realItems,:));
player_struct.LOTR.meanScore = player_struct.LOTR.score/22;

end

