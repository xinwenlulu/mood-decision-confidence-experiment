%function [alltest_t1, alltest_t2] = testblock(data, nBlocks)
function trialmtx = Xin_testblock(data, nBlocks)

%error('still need to add in WoF outcome preceding each stim!')

%%% Build list of all possible test comparisons from stimuli seen so far
%allstim = [data.picprob(:,[1 3]); data.picprob(:,[2 4])];
allstim = [data.picprob(:,[1 4 7]); data.picprob(:,[2 5 7]); data.picprob(:,[3 6 7])];
ind_allpairs = nchoosek(1:size(allstim,1),2); % Generate indices for all unique stimuli pairs
allpairs = allstim(ind_allpairs); % Retrieve the stimuli for all unique pairs

stimseen_probs=allstim(:,2); % Do again for probabilities of the stimuli
allprobs = stimseen_probs(ind_allpairs);
stimseen_wofs=allstim(:,3); % Do again for WoF outcomes
allwofs = stimseen_wofs(ind_allpairs);
alltest = [allpairs, allprobs, allwofs];

%%% Build list of test comparisons seen at test 1
%stimseen_t1 = data.picprob(1:2,:);
%stimseen_t1 = [stimseen_t1(:,[1 3]); stimseen_t1(:,[2 4])]; % With stim and prob, no WoF
%stimseen_t1 = [stimseen_t1(:,[1 3 5]); stimseen_t1(:,[2 4 6])]; % With WoF

%ind_allpairs_t1 = nchoosek(1:size(stimseen_t1,1),2); % Generate indices for all unique stimuli pairs
%allpairs_t1 = stimseen_t1(ind_allpairs_t1); % Retrieve the stimuli for all unique pairs
%stimseen_probs_t1=stimseen_t1(:,2); % Do again for probabilities of the stimuli
%allprobs_t1 = stimseen_probs_t1(ind_allpairs_t1);

%stimseen_wofs_t1=stimseen_t1(:,3); % Do again for probabilities of the stimuli
%allwofs_t1 = stimseen_wofs_t1(ind_allpairs_t1);
%alltest_t1 = [allpairs_t1, allprobs_t1, allwofs_t1];

%%% Build list of test comparisons for test 2 (dropping those seen at test 1) 
%alltest_t2 = alltest(find(~ismember(alltest, alltest_t1, 'rows')==1),:);


% Fill 2nd to 5th colums with shuffled test stimuli
%if nBlocks == 2
    %trialmtx = nan(size(alltest_t1,1),14);
    %alltest_t1 = alltest_t1(randperm(size(alltest_t1,1)),:); % Shuffle test stimuli
    %trialmtx(:,2:7) = alltest_t1; 
%elseif nBlocks == 3
    %trialmtx = nan(size(alltest_t2,1),14);
    %alltest_t2 = alltest_t2(randperm(size(alltest_t2,1)),:); % Shuffle test stimuli
    %trialmtx(:,2:7) = alltest_t2; 
%end

keyComparisons = [];
for n = 1:size(alltest,1)
    if alltest(n,3) == alltest(n,4) || alltest(n,5) == alltest(n,6)
        keyComparisons(end+1,:) = alltest(n,:);
    end
end

trialmtx = nan(size(keyComparisons,1),14);
keyComparisons = keyComparisons(randperm(size(keyComparisons,1)),:);
trialmtx(:,2:7) = keyComparisons;

flippedtrials = [ceil(size(trialmtx,1)/2)+1 : size(trialmtx,1)]; % Flip halfway to reduce bias of high prob to one side
trialmtx(flippedtrials,2:7) = [fliplr(trialmtx(flippedtrials,2:3)) fliplr(trialmtx(flippedtrials,4:5)) fliplr(trialmtx(flippedtrials,6:7))]; %flip lr the second half
trialmtx = trialmtx(randperm(size(trialmtx,1)),:); % Shuffle test stimuli again

trialmtx =[trialmtx; trialmtx]; % They do each comparison twice, once binary and once rating

starttime=time;
cgfont('Arial',data.FontSizeText); cgpencol(1,1,1); 
cgtext('Now we will show you those flowers in pairs',0,250);
cgtext('Pick the flower you think will give you the most gems',0,200);
cgtext('Watch out, you will see each pair only once.',0,150);
cgtext('So think well before making your choice.',0,100);

cgtext('In this round you can earn extra gems, if you do well.',0,-50);
cgtext('Outcomes will be added to your total earnings at the end.',0,-100);
cgflip(data.background);
wait(1000);
cgtext('Now we will show you those flowers in pairs',0,250);
cgtext('Pick the flower you think will give you the most gems',0,200);
cgtext('Watch out, you will see each pair only once.',0,150);
cgtext('So think well before making your choice.',0,100);

cgtext('In this round you can earn extra gems, if you do well.',0,-50);
cgtext('Outcomes will be added to your total earnings at the end.',0,-100);
cgtext('Click to continue',0,-200);
cgflip(data.background)
[~,~,~,mp]=cgmouse; mp=0; %no mouse press yet
while (mp ~=1) && (mp ~=4), %1 from leftmouse, 4 for rightmouse, 2 for middle mouse
    [~,~,~, mp] = cgmouse;
end;

cgfont('Arial',data.FontSizeTask); cgtext('+',0,0); cgflip(data.background);wait(data.time.delay);
stimDimensions = data.stimDimensions;

for itrial=1:size(trialmtx,1)
   
    disp(trialmtx)
    
    cgfont('Arial',data.FontSizeText+10); cgtext('Outcomes not revealed until the end.',0,300);
    cgfont('Arial',data.FontSizeTask); %cgtext('+',0,0);    
    
    % Draw left stim
    leftstim = data.pictures{trialmtx(itrial,2)};
    cgpencol(0,0,0); cgrect(-200,0,stimDimensions(1)+10,stimDimensions(2)+10)
    cgloadbmp(1,leftstim,stimDimensions(1),stimDimensions(2));cgdrawsprite(1,-200,0)
    
    % Draw right stim
    rightstim = data.pictures{trialmtx(itrial,3)};
    cgpencol(0,0,0); cgrect(200,0,stimDimensions(1)+10,stimDimensions(2)+10)
    cgloadbmp(2,rightstim,stimDimensions(1),stimDimensions(2)); cgdrawsprite(2,200,0)

    cgflip(data.background);
    
    choicemade = nan; stimchosen =nan; rating = nan; rt = nan; rating_rt = nan; %in case early exit
    tstart=time;
    %alltime = tstart;

    %% Get line rating instead of L/R click if more than half way through trials

    if itrial>ceil( size(trialmtx,1)/2 ) 
    
        originy = -200; 
        sx = [-200-stimDimensions(1)/2 200+stimDimensions(1)/2]; 
        sy = [originy, originy]; % for lefthandside
        timelimit=Inf;
        %randstart = 0.5;               %skip the randstart and just start from the middle
        mp2 = 0; 
        x = -280+560*rand; y = 0; %initialize mouse to center of screen at moment line appears
        %cgmouse((sx(2)-sx(1))*randstart + sx(1), 0); %mouse jumps to random location on line
        cgmouse(x, 0); %mouse jumps to random location on line

        while (mp2 ~= 1) && time < tstart + timelimit;

            cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); cgtext('Outcomes not revealed until the end.',0,300); 
            cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); cgtext('No preference',0,originy-30);
            cgtext('Strongly confident left',-350, originy-30); %cgtext('|',-200-stimDimensions(1)/2, originy);
            cgtext('Strongly confident right',350, originy-30); %cgtext('|',200+stimDimensions(1)/2,originy);

            % Keep drawing stim
            cgpencol(0,0,0); cgrect(-200,0,stimDimensions(1)+10,stimDimensions(2)+10)
            cgdrawsprite(1,-200,0) % Left stim

            cgrect(200,0,stimDimensions(1)+10,stimDimensions(2)+10)
            cgdrawsprite(2,200,0) % Right stim

            % Draw rating line
            cgpencol(0,0,0); cgdraw(sx(1),sy(1),sx(2),sy(2)); 
            cgpencol(1,1,1); cgellipse(sx(1),sy(1),10,10,'f');
            cgellipse(sx(2),sy(2),10,10,'f');

            % Update mouse position
            [x, ~, ~, mp2] = cgmouse; y = originy;
            if x < sx(1), x = sx(1); end
            if x > sx(2), x = sx(2); end;
            cgpencol(1,1,0); cgellipse(x,y,10,10,'f');
            cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); %cgtext('+',0,0); 
            cgflip(data.background);
            
            
        end
        
        % After choice leave up rating (no stim)
        cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); cgtext('Outcomes not revealed until the end.',0,300); 
            cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); cgtext('No preference',0,originy-30);
            cgtext('Strongly confident left',-350, originy-30); %cgtext('|',-200-stimDimensions(1)/2, originy);
            cgtext('Strongly confident right',350, originy-30); %cgtext('|',200+stimDimensions(1)/2,originy);
            
        cgpencol(0,0,0); cgdraw(sx(1),sy(1),sx(2),sy(2)); 
        cgpencol(1,1,1); cgellipse(sx(1),sy(1),10,10,'f');
        cgellipse(sx(2),sy(2),10,10,'f');
        
        % Update mouse position
            [x, ~, ~, mp2] = cgmouse; y = originy;
            if x < sx(1), x = sx(1); end
            if x > sx(2), x = sx(2); end;
            cgpencol(1,1,0); cgellipse(x,y,10,10,'f');
            cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); %cgtext('+',0,0);
        
        rating = (x - sx(1)) / (sx(2) - sx(1)); %0 to 1
        rating_rt = (time - tstart)/1000; %in s since end of waittime
        
        if rating < .5, choicemade = 1;
        elseif rating > .5, choicemade = 2;
        else choicemade = 0;
        end
        
         cgpencol(1,.8,0); 
        if choicemade == 1
            cgrect(-200,0,stimDimensions(1)+20,stimDimensions(2)+20);
            cgdrawsprite(1,-200,0)
        elseif choicemade == 2
            cgrect(200,0,stimDimensions(1)+20,stimDimensions(2)+20)
            cgdrawsprite(2,200,0)
        end
        
         cgflip(data.background);
         wait(data.time.outcometime)
         
         cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); %cgtext('Outcomes not revealed until the end.',0,300); 
         cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); 
         cgtext('No preference',0,originy-30);
         cgtext('Strongly confident left',-350, originy-30); %cgtext('|',-200-stimDimensions(1)/2, originy);
         cgtext('Strongly confident right',350, originy-30); %cgtext('|',200+stimDimensions(1)/2,originy);
      
         % Draw rating line
         cgpencol(0,0,0); cgdraw(sx(1),sy(1),sx(2),sy(2)); 
         cgpencol(1,1,1); cgellipse(sx(1),sy(1),10,10,'f');
         cgellipse(sx(2),sy(2),10,10,'f');
    
    else % Trial does not need a line rating

        originy = -180; sx = [-200-stimDimensions(1)/2 200+stimDimensions(1)/2]; sy = [originy, originy]; % for lefthandside
        timelimit=Inf;
        %randstart = 0.5;               %skip the randstart and just start from the middle
        mp2 = 0; x = -280+560*rand; y = 0; %initialize mouse to center of screen at moment line appears
        %cgmouse((sx(2)-sx(1))*randstart + sx(1), 0); %mouse jumps to random location on line
        cgmouse(x, 0); %mouse jumps to random location on line

        while (mp2 ~= 1) && time < tstart + timelimit;

            cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); cgtext('Outcomes not revealed until the end.',0,300); 
            %cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); cgtext('No preference',0,originy-30);
            %cgtext('Strongly confident left',-350, originy-30); %cgtext('|',-200-stimDimensions(1)/2, originy);
            %cgtext('Strongly confident right',350, originy-30); %cgtext('|',200+stimDimensions(1)/2,originy);

            % Keep drawing stim
            cgpencol(0,0,0); cgrect(-200,0,stimDimensions(1)+10,stimDimensions(2)+10)
            cgdrawsprite(1,-200,0) % Left stim

            cgrect(200,0,stimDimensions(1)+10,stimDimensions(2)+10)
            cgdrawsprite(2,200,0) % Right stim

            % Draw rating line
            %cgpencol(0,0,0); cgdraw(sx(1),sy(1),sx(2),sy(2)); 
            %cgpencol(1,1,1); cgellipse(sx(1),sy(1),10,10,'f');
            %cgellipse(sx(2),sy(2),10,10,'f');

            % Update mouse position
            [x, ~, ~, mp2] = cgmouse; y = 0;
            if x < sx(1), x = sx(1); end
            if x > sx(2), x = sx(2); end;
            cgpencol(1,1,0); cgellipse(x,y,10,10,'f');
            %cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); %cgtext('+',0,0); 
            cgflip(data.background);
            
            
        end
        
        % After choice leave up rating (no stim)
        cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); cgtext('Outcomes not revealed until the end.',0,300); 
        
        % Update mouse position
            [x, ~, ~, mp2] = cgmouse; y = originy;
            if x < sx(1), x = sx(1); end
            if x > sx(2), x = sx(2); end;
            
        rating = (x - sx(1)) / (sx(2) - sx(1)); %0 to 1
        rating_rt = (time - tstart)/1000; %in s since end of waittime
        
        if rating < .5, choicemade = 1;
        elseif rating > .5, choicemade = 2;
        else choicemade = 0;
        end
        
        %[~,~,~,mp]=cgmouse; mp=0; %no mouse press yet
        %while (mp ~=1) && (mp ~=4), %1 from leftmouse, 4 for rightmouse, 2 for middle mouse
            %[~,~,~, mp] = cgmouse;
        %end;
        %choicemade = sqrt(mp);
        %rt = (time - tstart)/1000; %rt in s

        cgpencol(1,.8,0); 
        if choicemade == 1
            cgrect(-200,0,stimDimensions(1)+20,stimDimensions(2)+20);
            cgdrawsprite(1,-200,0)
        elseif choicemade == 2
            cgrect(200,0,stimDimensions(1)+20,stimDimensions(2)+20)
            cgdrawsprite(2,200,0)
        end
        %cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); cgtext('Outcomes not revealed until the end.',0,300);
        %cgfont('Arial',data.FontSizeTask); %cgtext('+',0,0);
        cgflip(data.background);
        wait(data.time.outcometime);            
    end
    
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); %cgtext('+',0,0);
    cgfont('Arial',data.FontSizeText+10); cgtext('Outcomes not revealed until the end.',0,300);
    cgflip(data.background); 
    wait(data.iti);
        
    % Logging
    if choicemade > 0, stimchosen = trialmtx(itrial,1+choicemade);
    elseif choicemade == 0, stimchosen = 0;
    else stimchosen = nan;
    end
    
    trialmtx(itrial,1) = itrial;
   
    if isnan(rating)
        trialmtx(itrial,9) = stimchosen; % Don't look up stim chosen if rating was exactly inbetween.
        trialmtx(itrial,8) = choicemade;
    else
        trialmtx(itrial,12) = stimchosen; % Don't look up stim chosen if rating was exactly inbetween.
    end
    
    trialmtx(itrial,10) = rt;
    
    trialmtx(itrial,11) = rating;
    trialmtx(itrial,13) = rating_rt;
    trialmtx(itrial,14) = (time-starttime)/1000;
    save(data.savefile, 'data');
    
end

cgfont('Arial',data.FontSizeTask); cgpencol(1,1,1); %cgtext('+',0,0);
cgflip(data.background);
wait(data.time.delay);