function [experiment_data] = ConfidenceSampling_task_sascha(data, blocknum, wof, showpounds)

%%%%COMMENT THIS OUT%%%%
%blocknum = 1; wof = 0; data.ntrials=30; %data.picprob = [3 7 0.3 0.65]; %for testing
%%%%%%%%

cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); cgtext('+',0,0);
cgflip(data.background);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%current design - 15 trials, 5 ratings %indexratings = 5:5:30;
% NEW current design - 20 trials, 5 ratings %indexratings = 5:5:30;

%t=[3 3 3; 3 4 4]; for n=1:size(t,2), t(:,n)=t(randperm(2),n); end; indexratings = cumsum(t(:)); %trials with ratings after them
%t = [3 3 4 4 4]; indexratings = cumsum(t(randperm(length(t)))); % Five happiness ratings
%t = [3 3 3 4 4 4]; indexratings = cumsum(t(randperm(length(t)))); % Six happiness ratings
%t = [3 3 4 4 5 5 5]; indexratings = cumsum(t(randperm(length(t)))); % Seven happiness ratings
%picprob is the stimulus numbers and their probabilities - shuffle mtx

%indexratings2 = [2 data.ntrials/2 data.ntrials-2]; % 3 happiness ratings
indexratings2 = [14 28]; % 2 happiness ratings, plus the one before and one after the block (see below)
%t = [3 3 4 4 5 5 5]; indexratings = cumsum(t(randperm(length(t)))); % Seven confidence ratings
%t = [3 3 4 4 4 4 4 5]; indexratings = cumsum(t(randperm(length(t)))); % Eight confidence ratings
%t = [3 3 3 4 4 4 4 5 5]; indexratings = cumsum(t(randperm(length(t)))); % Nine confidence ratings
t = [3 3 4 4 4 4 4 4 5 5]; indexratings = cumsum(t(randperm(length(t)))); % Ten confidence ratings
%indexratings = [2 5 9 14 17 22 26 31 36 41]; % Ten confidence ratings

%%% Build in forced choice and shuffle
ind = reshape(1:data.ntrials-6,6,(data.ntrials-6)/6); % Subtract six trials here to ensure first two and last five are free choice (below)
for n=1:size(ind,2), ind(:,n)=ind(randperm(6),n); end;
%done = 0; while ~done, if ind(end)>12, done=1; else ind(:,end)=ind(randperm(5),end); end; end; %constrain so ends on a choice
t = repmat([1 0; 0 1; 1 1; 1 1; 1 1; 1 1],6,1); t = t(ind,:); %shuffle forced trials locally
%t=[t;ones(2)]; % Add on two extra free choice at end
%t(1:2,:) = ones(2,2); % Make first two trials always free choice
t2=[[1 1; 1 1; 1 0; 0 1]]; t2(randsample(1:length(t2),length(t2)),:); %add in two more forced choice trials
t=[[1 1];t2;t;[1 1]]; % Add on one free choice trials at the beginning and one at the end each to get back to 42 trials


% Drop columns with WoF outcomes
picprob_only = data.picprob(:,1:4); 
% Drop columns with WoF outcomes

%%% Decide when to flip left / right
trialmtx = ones(data.ntrials,1)*picprob_only(blocknum,:).*[t t];
if rand>0.5, trialmtx = trialmtx(:,[2 1 4 3]); end; %flip lr everything half of the time so rich could start on either side
%flippedtrials = [indexratings(round(length(indexratings)/2))+1:indexratings(end)];
%flippedtrials = [indexratings(2)+1:indexratings(4)]; % Flip after 2nd and 4th happiness ratings (but ?imbalance of trials on one side...)
% flippedtrials = [round(data.ntrials/2)+1 : data.ntrials]; % Flip exactly halfway regardless of happiness rating
flippedtrials = [indexratings(2):indexratings(4), indexratings(6):data.ntrials];length(flippedtrials) % Flip three times
trialmtx(flippedtrials,:) = [fliplr(trialmtx(flippedtrials,1:2)) fliplr(trialmtx(flippedtrials,3:4))]; %flip lr the second half

%shuffle = [0;0;1;1]; shuffled = repmat(shuffle,1,ceil(data.ntrials/size(shuffle,1)) );
%for n=1:size(shuffled,2), shuffled(:,n)=shuffled(randperm(size(shuffle,1)),n);end; %%%%NEW%%%%
%shuffled = shuffled(:); shuffled=shuffled(1:data.ntrials);

side=[0 1]; shuffled=[];
for i=1:ceil(data.ntrials/2)
    %shuffled(i,:) = side(randperm(2) );
    shuffled = [shuffled, side(randperm(2) )];
end
trialmtx(shuffled==1,:) = [fliplr(trialmtx(shuffled==1,1:2)) fliplr(trialmtx(shuffled==1,3:4))]; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%risktrials = 5:5:30;

%[~,~,~,mp] = cgmouse;
tstart = time;                      %get time in ms from first time call
nrating = 1; %nratingrisk = 1;
nrating2 = 1; %nratingrisk = 1;


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

% reward probabilities [60 20]
himax = [1;1;1;0;0]
lomax = [0;0;0;0;1]

% reward probabilities [66 33]
% himax1 = [1;1;0]
% himax2 = [1;1;1;1;0;0]
% hitotal= cat(1,himax1, himax2)
% lomax1 = [0;0;1]
% lomax2 = [0;0;0;0;1;1]
% lototal=

%predetermine rewards for each option
%hireward = repmat(himax,1,data.ntrials/size(himax,1)); loreward = repmat(lomax,1,data.ntrials/size(himax,1)); %%%%NEW%%%%
hireward = repmat(himax,1,ceil(data.ntrials/size(himax,1)) ); loreward = repmat(lomax,1,ceil( data.ntrials/size(himax,1) )); % NEWER: rounds up

for n=1:size(hireward,2), hireward(:,n)=hireward(randperm(size(himax,1)),n); loreward(:,n)=loreward(randperm(size(himax,1)),n); end; %%%%NEW%%%%
hireward = hireward(:); loreward = loreward(:); 
%hireward = hireward(1:data.ntrials)'; loreward = loreward(1:data.ntrials)'; % Matrix is too big now due to round up, so just take enough for ntrials
%length(find(hireward==1))/22

chosehi = 1; choselo = 1; %how many of each chosen so far %%%%NEW%%%%



%%%%%%%

for itrial = 1:data.ntrials
    disp(experiment_data(:,1:23));
    
    if itrial == 1 % only for first trial
        [experiment_data(itrial,19), experiment_data(itrial,20), experiment_data(itrial,21), alltime] = happy_rating_sascha('nowhappy',data.time.ratewait,data); %respkeys, waittime, resptime, startlocation, pixperkey);
         experiment_data(itrial,22) = NaN;
        cgpencol(1,1,1); cgtext('+', 0, 0); cgflip(data.background); 
    else
        experiment_data(itrial,[19:22]) = NaN;
    end    
    
    allowleftclick = -1; allowrightclick = -1; % by default, no mouse presses are defined as acceptable
    stimDimensions = data.stimDimensions;
    
    % check if free choice, forced left, or forced right trial
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); cgtext('+',0,0);
    if experiment_data(itrial,2)>0,
        allowleftclick=1; % allow click left
        pic = data.pictures{experiment_data(itrial,2)};
        cgloadbmp(1,pic,stimDimensions(1),stimDimensions(2));
        
        % show Bonus box in second learning block
        if showpounds == 1   
           cgpencol(0,0,0);
           cgrect(420,350,200,100);
           if data.wofscore >= 0 
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('  £',400,330);       
            cgalign('l','c')
            cgtext(num2str(data.wofscore),428,330);
            cgalign('c','c');
           else
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('- £',400,330);       
            cgalign('l','c')
            cgtext(num2str(abs(data.wofscore)),428,330);
            cgalign('c','c');
          end
        else
        end
        
        cgrect(-200,0,stimDimensions(1),stimDimensions(2));
        cgdrawsprite(1,-200,0);
    end
    if experiment_data(itrial,3)>0,
        allowrightclick=4; % allow click right
        pic = data.pictures{experiment_data(itrial,3)};
        cgloadbmp(2,pic,stimDimensions(1),stimDimensions(2));

                
        % show Bonus box in second learning block
        if showpounds == 1  
            cgpencol(0,0,0);
            cgrect(420,350,200,100);
          if data.wofscore >= 0 
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('  £',400,330);       
            cgalign('l','c')
            cgtext(num2str(data.wofscore),428,330);
            cgalign('c','c');
           else
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('- £',400,330);       
            cgalign('l','c')
            cgtext(num2str(abs(data.wofscore)),428,330);
            cgalign('c','c');
          end
        else
        end
        
        cgrect(200,0,stimDimensions(1),stimDimensions(2));
        cgdrawsprite(2,200,0)
    end
    cgflip(data.background);
    
    % get the choice
    tstart2 = time; experiment_data(itrial,14) = NaN; %experiment_data(itrial,14) = (tstart2 - tstart)/1000; %time since start of this task - %%%%NEW%%%%edited
    [~,~,~,mp]=cgmouse; mp=0; %no mouse press yet
    while (mp ~=allowleftclick) && (mp ~=allowrightclick), %1 from leftmouse, 4 for rightmouse - if a free choice trial
        [~,~,~, mp] = cgmouse;
    end;
    rt = (time - tstart2)/1000; experiment_data(itrial,15) = NaN; %experiment_data(itrial,15) = (time - tstart)/1000; %rt in s %%%%NEW%%%%edited
    choicemade = sqrt(mp);
    experiment_data(itrial,6) = choicemade; %which button pressed (1 left, 2 right) or 0 for none
    
    
    experiment_data(itrial,7) = experiment_data(itrial,1+choicemade); %%%%NEW%%%%
    if experiment_data(itrial,3+choicemade)>0.5, %if chose hireward option %%%%NEW%%%%
        experiment_data(itrial,8) = data.rewardsize*hireward(chosehi); chosehi=chosehi+1; %%%%NEW%%%%
    else %%%%NEW%%%%
        experiment_data(itrial,8) = data.rewardsize*loreward(choselo); choselo=choselo+1; %%%%NEW%%%%
    end; %%%%NEW%%%%
    %    experiment_data(itrial,8)=data.rewardsize*(experiment_data(itrial,3+choicemade)>rand); %determine outcome %%%%NEW%%%%
    %Add variable in col 23 cumulatively summing rewards
     experiment_data(itrial,23) = data.currentscore + experiment_data(itrial,8);
   
    experiment_data(itrial,9) = rt;   %rt in s
    if choicemade==1
              
        % show Bonus box in second learning block
        if showpounds == 1  
           cgpencol(0,0,0);
           cgrect(420,350,200,100);
           if data.wofscore >= 0 
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('  £',400,330);       
            cgalign('l','c')
            cgtext(num2str(data.wofscore),428,330);
            cgalign('c','c');
           else
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('- £',400,330);       
            cgalign('l','c')
            cgtext(num2str(abs(data.wofscore)),428,330);
            cgalign('c','c');
          end
        else
        end
         
        cgpencol(1,.8,0)
        cgrect(-200,0,stimDimensions(1)+20,stimDimensions(2)+20);
        if experiment_data(itrial,3)>0; cgdrawsprite(2,200,0); end %Other stimulus does not disappear
        cgdrawsprite(1,-200,0);
    elseif choicemade==2
        % show Bonus box in second learning block
        if showpounds == 1  
            cgpencol(0,0,0);
            cgrect(420,350,200,100);
           if data.wofscore >= 0 
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('  £',400,330);       
            cgalign('l','c')
            cgtext(num2str(data.wofscore),428,330);
            cgalign('c','c');
          else

            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('- £',400,330);       
            cgalign('l','c')
            cgtext(num2str(abs(data.wofscore)),428,330);
            cgalign('c','c');
          end
        else
        end
         
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
    if experiment_data(itrial,8)>0 % Win
        cgloadbmp(3,img,rewardstimDimensions(1),rewardstimDimensions(2))
        cgfont('Arial',data.FontSizeTask);
        cgpencol(1,1,1); % White text
        if choicemade==1
             cgpencol(1,.8,0)
             cgrect(-200,0,stimDimensions(1)+20,stimDimensions(2)+20);
            cgdrawsprite(1,-200,0)
            if experiment_data(itrial,3)>0; cgdrawsprite(2,200,0); end %Other stimulus does not disappear
            cgdrawsprite(3,-200,0)
            cgpencol(1,1,1)
            cgtext(num2str(data.rewardsize), -200, 0);

            
        % show Bonus box in second learning block
        if showpounds == 1   
            cgpencol(0,0,0);
            cgrect(420,350,200,100);
           if data.wofscore >= 0 
       
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('  £',400,330);       
            cgalign('l','c')
            cgtext(num2str(data.wofscore),428,330);
            cgalign('c','c');
           else
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('- £',400,330);       
            cgalign('l','c')
            cgtext(num2str(abs(data.wofscore)),428,330);
            cgalign('c','c');
          end
        else
        end
            
        elseif choicemade==2
             cgpencol(1,.8,0)
             cgrect(200,0,stimDimensions(1)+20,stimDimensions(2)+20);
            if experiment_data(itrial,2)>0; cgdrawsprite(1,-200,0); end %Other stimulus does not disappear
            cgdrawsprite(2,200,0)
            cgdrawsprite(3,200,0)
            cgpencol(1,1,1)
            cgtext(num2str(data.rewardsize), 200, 0);
     % show Bonus box in second learning block
        if showpounds == 1   
            cgpencol(0,0,0);
            cgrect(420,350,200,100);
           if data.wofscore >= 0 
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('  £',400,330);       
            cgalign('l','c')
            cgtext(num2str(data.wofscore),428,330);
            cgalign('c','c');
           else
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('- £',400,330);       
            cgalign('l','c')
            cgtext(num2str(abs(data.wofscore)),428,330);
            cgalign('c','c');
          end
        else
        end
            
        end
    elseif experiment_data(itrial,8)==0 % No win
        cgfont('Arial',data.FontSizeTask);
        if choicemade==1
            cgpencol(1,.8,0)
            cgrect(-200,0,stimDimensions(1)+20,stimDimensions(2)+20);
            cgdrawsprite(1,-200,0)
            if  experiment_data(itrial,3)>0; cgdrawsprite(2,200,0); end %Other stimulus does not disappear
            cgpencol(.2,.2,.2); % Grey background
            cgrect(-200,0,rewardstimDimensions(1),rewardstimDimensions(2))
            cgpencol(1,0,0); % Red cross
            cgtext('0', -200, 0);
  
        % show Bonus box in second learning block
        if showpounds == 1 
            cgpencol(0,0,0);
            cgrect(420,350,200,100);
           if data.wofscore >= 0 
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('  £',400,330);       
            cgalign('l','c')
            cgtext(num2str(data.wofscore),428,330);
            cgalign('c','c');
          else

            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('- £',400,330);       
            cgalign('l','c')
            cgtext(num2str(abs(data.wofscore)),428,330);
            cgalign('c','c');
          end
        else
        end
            
        elseif choicemade==2
            cgpencol(1,.8,0)
            cgrect(200,0,stimDimensions(1)+20,stimDimensions(2)+20);
            if experiment_data(itrial,2)>0; cgdrawsprite(1,-200,0); end %Other stimulus does not disappear
            cgdrawsprite(2,200,0)
            cgpencol(.2,.2,.2); % Grey background
            cgrect(200,0,rewardstimDimensions(1),rewardstimDimensions(2))
            cgpencol(1,0,0); % Red cross
            cgtext('0', 200, 0);
            
        % show Bonus box in second learning block
        if showpounds == 1   
            cgpencol(0,0,0);
            cgrect(420,350,200,100);
           if data.wofscore >= 0 
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('  £',400,330);       
            cgalign('l','c')
            cgtext(num2str(data.wofscore),428,330);
            cgalign('c','c');
           else
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('- £',400,330);       
            cgalign('l','c')
            cgtext(num2str(abs(data.wofscore)),428,330);
            cgalign('c','c');
          end
        else
        end
            
        end
    end
    cgflip(data.background); experiment_data(itrial,16) = NaN; %experiment_data(itrial,16) = (time - tstart)/1000; %outcome revealed %%%%NEW%%%%edited
    wait(data.time.outcometime);
    cgfont('Arial',data.FontSizeTask);
    cgpencol(1,1,1); cgtext('+', 0, 0);

        % show Bonus box in second learning block
        if showpounds == 1
            cgpencol(0,0,0);
            cgrect(420,350,200,100);
           if data.wofscore >= 0 
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('  £',400,330);       
            cgalign('l','c')
            cgtext(num2str(data.wofscore),428,330);
            cgalign('c','c');
           else
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);
            cgtext('- £',400,330);       
            cgalign('l','c')
            cgtext(num2str(abs(data.wofscore)),428,330);
            cgalign('c','c');
          end
        else
        end
    
    cgflip(data.background); experiment_data(itrial,17) = NaN; %experiment_data(itrial,17) = (time - tstart)/1000; %iti starts %%%%NEW%%%%edited
    wait(data.time.iti+round((data.time.maxiti-data.time.iti)*rand)); %jitter between min and max iti
    experiment_data(itrial,18) = NaN; %experiment_data(itrial,18) = (time - tstart)/1000; %trials end - could jitter %%%%NEW%%%%edited
    
    % if needed ask for confidence - 700ms ITI after 700ms of rating shown
    if itrial == indexratings(nrating), % || itrial == data.ntrials
        if nrating<length(indexratings)
            nrating = nrating + 1;
        end
        [experiment_data(itrial,10), experiment_data(itrial,11), experiment_data(itrial,12), alltime] = confidence_rating_sascha('choiceconfidence',data.time.ratewait,data); %respkeys, waittime, resptime, startlocation, pixperkey);
        cgpencol(1,1,1); cgtext('+', 0, 0); cgflip(data.background);   
    else
        experiment_data(itrial,[10:12]) = NaN;
    end

    % if needed ask for happiness - 700ms ITI after 700ms of rating shown
    if itrial == indexratings2(nrating2), % || itrial == data.ntrials
        if nrating2<length(indexratings2)
            nrating2 = nrating2 + 1;
        end
        [experiment_data(itrial,19), experiment_data(itrial,20), experiment_data(itrial,21), alltime] = happy_rating_sascha('nowhappy',data.time.ratewait,data); %respkeys, waittime, resptime, startlocation, pixperkey);
        experiment_data(itrial,22) = NaN;
        cgpencol(1,1,1); cgtext('+', 0, 0); cgflip(data.background);
%         experiment_data(itrial,19) = (alltime(1) - tstart)/1000; %happiness question onset time
%         experiment_data(itrial,20) = (alltime(2) - tstart)/1000; %happiness line onset
%         experiment_data(itrial,21) = experiment_data(itrial,12) + experiment_data(itrial,19); %happines rating made
%         experiment_data(itrial,22) = (time - tstart)/1000;       %ITI (0.7s) and 1s of looking at rating already built in  
    else
    end    
    
    if itrial == 42,
        cgflip(data.background);
        cgpencol(1,1,1); cgfont('Arial',data.FontSizeText);
        cgtext('This round is now over.',0,0);
        cgflip(data.background);
        wait(data.time.instructiontime);
        cgtext('This round is now over.',0,0);
        cgtext('Click to continue',0,-250);
        cgflip(data.background);
        [~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;
 
        [experiment_data(itrial,19), experiment_data(itrial,20), experiment_data(itrial,21), alltime] = happy_rating_sascha('nowhappy',data.time.ratewait,data); %respkeys, waittime, resptime, startlocation, pixperkey);
         experiment_data(itrial,22) = NaN;
        cgpencol(1,1,1); cgtext('+', 0, 0); cgflip(data.background);
        disp(experiment_data(:,1:23));
    else
    end
    data.maintask{blocknum} = experiment_data;
    save(data.savefile, 'data');
    
    % if needed, offer risky gamble
    %     if itrial == risktrials(nratingrisk) || itrial == data.ntrials
    %         if nratingrisk<length(risktrials)
    %             nratingrisk = nratingrisk + 1;
    %         end
    
    
    %  %%%%%%%%%%%%%%%%%%%%NEEDS TO BE FIXED%%%%%%%%%%%%%%%%%%%%%%%%%
    %  %THIS SHOULD ALL GO IN COLUMNS 23-36 IF 7x2 COLUMNS FOR RISK
    %         %next 5 lines here have been updated and do_risktrial
    %         %liam will provide 7 numbers for each trial
    %         %needs to have the pictures for the correct stimulus numbers like:         picleft = data.pictures{experiment_data(itrial,2)};
    %         %with the corresponding probability experiment_data(itrial,3+choicemade)
    %         %but these flip left-right halfway so really need to do this based on the stimulus and not the location
    %         [experiment_data(itrial,23), experiment_data(itrial,24), experiment_data(itrial,25), experiment_data(itrial,26), alltime] = do_risktrial('risky',data.time.ratewait,data); %respkeys, waittime, resptime, startlocation, pixperkey);
    %         cgpencol(1,1,1); cgtext('+', 0, 0); Do you know which animal has the most gems?(data.background);
    %         experiment_data(itrial,27) = (alltime(1) - tstart)/1000; %risky choice made
    %         experiment_data(itrial,28) = (alltime(2) - tstart)/1000; %not sure what any of this is
    %         experiment_data(itrial,29) = experiment_data(itrial,25) + experiment_data(itrial,27); %choice made
    %         experiment_data(itrial,30) = (time - tstart)/1000;
    % [experiment_data(itrial, 23), experiment_data(itrial, 24),experiment_data(itrial, 25),experiment_data(itrial, 26),experiment_data(itrial, 27),experiment_data(itrial, 28),experiment_data(itrial, 29),experiment_data(itrial, 30),experiment_data(itrial, 31),experiment_data(itrial, 32),experiment_data(itrial, 33),experiment_data(itrial, 34), alltime] = do_risktrial(data, experiment_data, blocknum, itrial, 3); % Risk trial with fixed gamble of 3
    %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %         data.maintask{blocknum} = experiment_data;
    %         %save data after every happy rating
    %         eval(sprintf('save %s data',fullfile(data.dir,data.filename))); %save data
    %     else
    %         experiment_data(itrial,[23:36]) = NaN;
    %     end
end
wait(data.time.iti);