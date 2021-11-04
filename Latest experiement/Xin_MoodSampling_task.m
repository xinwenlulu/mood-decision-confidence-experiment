function [experiment_data]   = Xin_MoodSampling_task(data, blocknum, wof)

%%%%COMMENT THIS OUT%%%%
%blocknum = 1; wof = 0; data.ntrials=30; %data.picprob = [3 7 0.3 0.65]; %for testing
%%%%%%%%

cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); cgtext('+',0,0);
cgflip(data.background);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

happy1 = 1;
happy2 = randi([8 12],1);
happy3 = randi([ceil(data.ntrials*2/3)-2 ceil(data.ntrials*2/3)+2],1);
happy4 = randi([data.ntrials-3 data.ntrials-2],1);
indexratings = [happy1 happy2 happy3 happy4]; %5th happy rating at the end


% Drop columns with WoF outcomes

picprob_only = data.picprob(:,1:6);

block_picprob = picprob_only(blocknum, :);
% Create 3 pairs from the 3 stimuli in the current block
pairs = zeros(3,4);
for n = 1:3 %3 pairs 
    if n ~= 3
        pairs(n,:) = [block_picprob(1,n:n+1) block_picprob(1,n+3:n+4)];
    else
        pair_index = [1 3 4 6];
        pairs(n,:) = block_picprob(pair_index);
    end
end

learningPairs = pairs;

for n = 1:3
    if (pairs(n,3) == 2 && pairs(n,4) == 3) || (pairs(n,3) == 3 && pairs(n,4) == 2)
        learningPairs(n,:) = [];
        mHpair = pairs(n,:);
    end
    
end

trialmtx = ones(data.ntrials/4,4);
ntrial_persection = (data.ntrials)/4;%4 sections, one mHpair per section
npair_persection = (ntrial_persection-1)/2;

for n = 1:npair_persection
    trialmtx(2*n-1:2*n,:) = learningPairs;
end

trial_index = randperm(ntrial_persection-1);
unshuffledmtx = trialmtx;
for n = 1:ntrial_persection-1
    trialmtx(n,:) = unshuffledmtx(trial_index(n),:);
end
trialmtx(end,:) = mHpair;
trialmtx = repmat(trialmtx,4,1);
data.subjectid = '10'
subjectid = str2num(data.subjectid);

% randomise the location of the bad option in the first trial (counterbalancing the order of two experiments)
    if blocknum == 1
        if  mod(subjectid, 4) == 0
            flippedtrials = [ceil(data.ntrials/8): ceil(data.ntrials/4) ceil(data.ntrials*3/8): ceil(data.ntrials/2) ceil(data.ntrials*5/5):ceil(data.ntrials*3/4) ceil(data.ntrials*7/8):data.ntrials]; 
        elseif mod(subjectid, 4) == 1
            flippedtrials = [1 ceil(data.ntrials/8): ceil(data.ntrials/4) ceil(data.ntrials*3/8): ceil(data.ntrials/2) ceil(data.ntrials*5/5):ceil(data.ntrials*3/4) ceil(data.ntrials*7/8):data.ntrials]; 
        elseif mod(subjectid, 4) == 2
            flippedtrials = [1 ceil(data.ntrials/8): ceil(data.ntrials/4) ceil(data.ntrials*3/8): ceil(data.ntrials/2) ceil(data.ntrials*5/5):ceil(data.ntrials*3/4) ceil(data.ntrials*7/8):data.ntrials]; 
        else %mod(id,4) ==3
            flippedtrials = [ceil(data.ntrials/8): ceil(data.ntrials/4) ceil(data.ntrials*3/8): ceil(data.ntrials/2) ceil(data.ntrials*5/5):ceil(data.ntrials*3/4) ceil(data.ntrials*7/8):data.ntrials]; 
        end
    else %second block
        if  mod(subjectid, 4) == 0
            flippedtrials = [1 ceil(data.ntrials/8): ceil(data.ntrials/4) ceil(data.ntrials*3/8): ceil(data.ntrials/2) ceil(data.ntrials*5/5):ceil(data.ntrials*3/4) ceil(data.ntrials*7/8):data.ntrials]; 
        elseif mod(subjectid, 4) == 1
            flippedtrials = [1 ceil(data.ntrials/8): ceil(data.ntrials/4) ceil(data.ntrials*3/8): ceil(data.ntrials/2) ceil(data.ntrials*5/5):ceil(data.ntrials*3/4) ceil(data.ntrials*7/8):data.ntrials]; 
        elseif mod(subjectid, 4) == 2
            flippedtrials = [ceil(data.ntrials/8): ceil(data.ntrials/4) ceil(data.ntrials*3/8): ceil(data.ntrials/2) ceil(data.ntrials*5/5):ceil(data.ntrials*3/4) ceil(data.ntrials*7/8):data.ntrials]; 
        else %mod(id,4) ==3
            flippedtrials = [ceil(data.ntrials/8): ceil(data.ntrials/4) ceil(data.ntrials*3/8): ceil(data.ntrials/2) ceil(data.ntrials*5/5):ceil(data.ntrials*3/4) ceil(data.ntrials*7/8):data.ntrials];  
        end
    end


%%% Decide when to flip left / right
%flippedtrials = [round(data.ntrials/4)+1 : round(data.ntrials/2) round(data.ntrials*3/4)+1:data.ntrials]; % Flip half the time
trialmtx(flippedtrials,:) = [fliplr(trialmtx(flippedtrials,1:2)) fliplr(trialmtx(flippedtrials,3:4))];

%for n = 1:data.ntrials
    %if rand>0.5
        %trialmtx(n,:) = trialmtx(n,[2 1 4 3]);
    %end
%end


%shuffle the trials in the trial matrix
%trial_index1 = randperm((data.ntrials-2)/2);
%trial_index2 = randperm((data.ntrials-2)/2);
%unshuffledmtx = trialmtx;
%for n = 1:data.ntrials-1
%    if n < data.ntrials/2
%        trialmtx(n,:) = unshuffledmtx(trial_index1(n),:);
%    elseif n > data.ntrials/2
%        trialmtx(n,:) = unshuffledmtx(trial_index2(n-data.ntrials/2),:);
%    end
%end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%risktrials = 5:5:30;

%[~,~,~,mp] = cgmouse;
tstart = time;                      %get time in ms from first time call
nrating = 1; %nratingrisk = 1;


wait(data.time.iti);

%columns: 1 trial num, 2-3 stim identities (0 forced choice, safe value if
%>1), 4-5 stim probs, 6 button pressed (1 or 2), 7 is stimulus chosen, 8
%outcome (in points), 9 is rt, 10-12 happiness rating, 13 is WoF (0 if
%none, otherwise outcome), 14-18 trial timing info, 19-22 happy rating
%info, 23-36 risky trial info
experiment_data = zeros(data.ntrials,24);
experiment_data(:,1) = 1:data.ntrials;
experiment_data(:,13) = wof; %wheel of fortune outcome %%%%NEW%%%%
experiment_data(:,2:5) = trialmtx; %zero anything that is empty - 0 stimulus and 0 prob
%%%%%%%%%%%%%%%%%%%%%%%%

%Reward distributions 1 = bad; 2 = goodvL; 3 = goodvH
bad = reward_distributions(40,10,15);
goodvL = reward_distributions(60,7,15);
goodvH = reward_distributions(60,25,15);

data.bad = bad;
data.goodvL = goodvL;
data.goodvH = goodvH;


%predetermine rewards for each option
rewardmL = repmat(bad,1,ceil(data.ntrials/size(bad,2)));
rewardmHvL = repmat(goodvL,1,ceil(data.ntrials/2/size(goodvL,2)));
rewardmHvH = repmat(goodvH,1,ceil(data.ntrials/2/size(goodvH,2)));


for n=1:ceil(data.ntrials/2/size(goodvL,2))
    rewardmHvL(15*n-14:15*n)=rewardmHvL(randperm(length(goodvL))); 
    rewardmHvH(15*n-14:15*n)=rewardmHvH(randperm(length(goodvH))); 
end

for n=1:ceil(data.ntrials/size(bad,2))
    rewardmL(15*n-14:15*n)=rewardmL(randperm(length(bad))); 
end

chosemL = 1; chosemHvL = 1; chosemHvH = 1; %how many of each chosen so far %%%%NEW%%%%

%%%%%%%

for itrial = 1:data.ntrials
    
       % if needed ask for happiness - 700ms ITI after 700ms of rating shown
    if itrial == indexratings(nrating), % || itrial == data.ntrials
        if nrating<length(indexratings)
            nrating = nrating + 1;
        end
        [experiment_data(itrial,10), experiment_data(itrial,11), experiment_data(itrial,12), alltime] = happy_rating('nowhappy',data.time.ratewait,data); %respkeys, waittime, resptime, startlocation, pixperkey);
        cgpencol(1,1,1); %cgtext('+', 0, 0); 
        cgflip(data.background);
        experiment_data(itrial,19) = (alltime(1) - tstart)/1000; %happy question onset time
        experiment_data(itrial,20) = (alltime(2) - tstart)/1000; %happy line onset
        experiment_data(itrial,21) = experiment_data(itrial,12) + experiment_data(itrial,19); %happy rating made
        experiment_data(itrial,22) = (time - tstart)/1000;       %ITI (0.7s) and 1s of looking at rating already built in
        data.maintask{blocknum} = experiment_data;
        save(data.savefile, 'data');    
    else
        experiment_data(itrial,[10:12 19:22]) = NaN;
    end
    
    disp(experiment_data(:,1:14));
    allowleftclick = 1; %allowrightclick = -1; 
    stimDimensions = data.stimDimensions;
    
    
   if blocknum == 2
        cgpencol(0,0,0);
        cgrect(420,350,200,100);
        cgpencol(1,1,1)
        cgfont('Arial',data.FontSizeText)
        cgtext('Bonus',420,365);
        if data.wofscore >= 0
            cgtext('  £',400,330);
            cgalign('l','c')
            cgtext(num2str(data.wofscore),428,330)
            cgalign('c','c');
        else
            cgtext('- £',400,330);        
            cgalign('l','c')
            cgtext(num2str(abs(data.wofscore)),428,330);
            cgalign('c','c');
        end    
   end
   
    if experiment_data(itrial,2)>0,
        pic = data.pictures{experiment_data(itrial,2)};
        cgloadbmp(1,pic,stimDimensions(1),stimDimensions(2));
        cgrect(-200,0,stimDimensions(1),stimDimensions(2));
        cgdrawsprite(1,-200,0);
    end
    if experiment_data(itrial,3)>0
        pic = data.pictures{experiment_data(itrial,3)};
        cgloadbmp(2,pic,stimDimensions(1),stimDimensions(2));
        cgrect(200,0,stimDimensions(1),stimDimensions(2));
        cgdrawsprite(2,200,0)
    end
    cgflip(data.background);
    
    originy = -180;
    pic = data.confidence{1};
    cgloadbmp(3,pic,350,50);
    cgrect(-175,originy-25,350,50);
    cgdrawsprite(3,-175,originy-25);
    pic = data.confidence{2};
    cgloadbmp(4,pic,350,50);
    cgrect(175,originy-25,350,50);
    cgdrawsprite(4,175,originy-25);
    
    % get the choice
    tstart2 = time; experiment_data(itrial,14) = (tstart2 - tstart)/1000; %time since start of this task - %%%%NEW%%%%edited
   
    sx = [-200-stimDimensions(1)/2 200+stimDimensions(1)/2]; sy = [originy, originy]; % for lefthandside
        timelimit=Inf;
        randstart = 0.5;               %skip the randstart and just start from the middle
        mp2 = 0; 
        x = -280+560*rand;
        %y = 0; %initialize mouse to center of screen at moment line appears
        cgmouse(x,0); %mouse jumps to random location on line
        rating = 0.5;

          while rating > 0.49 && rating < 0.51
                 while (mp2 ~= 1) && time < tstart + timelimit
                    cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); %cgtext('Outcomes added to total earnings',0,300); 
                    cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); 
                    cgtext('No preference',0,originy-80);
                    cgtext('Strongly confident left',-350, originy-80); %cgtext('|',-200-stimDimensions(1)/2, originy);
                    cgtext('Strongly confident right',350, originy-80); %cgtext('|',200+stimDimensions(1)/2,originy);


                    % Keep drawing stim
                    cgpencol(0,0,0); cgrect(-200,0,stimDimensions(1),stimDimensions(2))
                    cgdrawsprite(1,-200,0) % Left stim

                    cgrect(200,0,stimDimensions(1),stimDimensions(2))
                    cgdrawsprite(2,200,0) % Right stim
                    
                    cgrect(-175,originy-25,350,50)
                    cgdrawsprite(3,-175,originy-25)
                    cgrect(175,originy-25,.350,50)
                    cgdrawsprite(4,175,originy-25) 

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
             
                    if blocknum == 2
                        cgpencol(0,0,0);
                        cgrect(420,350,200,100);
                        cgpencol(1,1,1)
                        cgfont('Arial',data.FontSizeText)
                        cgtext('Bonus',420,365);
                        if data.wofscore >= 0
                            cgtext('  £',400,330);
                            cgalign('l','c')
                            cgtext(num2str(data.wofscore),428,330)
                            cgalign('c','c');
                        else
                            cgtext('- £',400,330);        
                            cgalign('l','c')
                            cgtext(num2str(abs(data.wofscore)),428,330);
                            cgalign('c','c');
                        end    
                    end
                    cgflip(data.background);        
                 end
                              
                cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); %cgtext('Outcomes added to total earnings',0,300);
                cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); 
                cgtext('No preference',0,originy-80);
                cgtext('Strongly confident left',-350, originy-80); %cgtext('|',-200-stimDimensions(1)/2, originy);
                cgtext('Strongly confident right',350, originy-80); %cgtext('|',200+stimDimensions(1)/2,originy);

                if blocknum == 2
                        cgpencol(0,0,0);
                        cgrect(420,350,200,100);
                        cgpencol(1,1,1)
                        cgfont('Arial',data.FontSizeText)
                        cgtext('Bonus',420,365);
                        if data.wofscore >= 0
                            cgtext('  £',400,330);
                            cgalign('l','c')
                            cgtext(num2str(data.wofscore),428,330)
                            cgalign('c','c');
                        else
                            cgtext('- £',400,330);        
                            cgalign('l','c')
                            cgtext(num2str(abs(data.wofscore)),428,330);
                            cgalign('c','c');
                        end    
                end
                    
                % Keep drawing stim
                cgpencol(0,0,0); cgrect(-200,0,stimDimensions(1),stimDimensions(2))
                cgdrawsprite(1,-200,0) % Left stim

                cgrect(200,0,stimDimensions(1),stimDimensions(2))
                cgdrawsprite(2,200,0) % Right stim

                cgrect(-175,originy-25,350,50)
                cgdrawsprite(3,-175,originy-25)
                cgrect(175,originy-25,.350,50)
                cgdrawsprite(4,175,originy-25) 

                % Draw rating line
                cgpencol(0,0,0); cgdraw(sx(1),sy(1),sx(2),sy(2)); 
                cgpencol(1,1,1); cgellipse(sx(1),sy(1),10,10,'f');
                cgellipse(sx(2),sy(2),10,10,'f');

                % Update mouse position
                [x, ~, ~, mp2] = cgmouse; y = originy;
                if x < sx(1), x = sx(1); end
                if x > sx(2), x = sx(2); end;
                rating = (x - sx(1)) / (sx(2) - sx(1)); %0 to 1
                rating_rt = (time - tstart)/1000; %in s since end of waittime
                
                cgpencol(1,1,0); cgellipse(x,y,10,10,'f');
                cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); %cgtext('+',0,0); 
               
                cgflip(data.background);    
                
                if rating > 0.49 && rating < 0.51
                    cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+5);
                    cgtext('No choice given.',0,250); 
                    cgpencol(1,1,1)
                    cgrect(0,-100,100,100);
                    cgpencol(0,0,0)
                    cgtext(num2str(0), 0, -100);
                    cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); %cgtext('Outcomes added to total earnings',0,300); 
                    cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); 
                    cgtext('No preference',0,originy-80);
                    cgtext('Strongly confident left',-350, originy-80); %cgtext('|',-200-stimDimensions(1)/2, originy);
                    cgtext('Strongly confident right',350, originy-80); %cgtext('|',200+stimDimensions(1)/2,originy);
 
                    if blocknum == 2
                        cgpencol(0,0,0);
                        cgrect(420,350,200,100);
                        cgpencol(1,1,1)
                        cgfont('Arial',data.FontSizeText)
                        cgtext('Bonus',420,365);
                        if data.wofscore >= 0
                            cgtext('  £',400,330);
                            cgalign('l','c')
                            cgtext(num2str(data.wofscore),428,330)
                            cgalign('c','c');
                        else
                            cgtext('- £',400,330);        
                            cgalign('l','c')
                            cgtext(num2str(abs(data.wofscore)),428,330);
                            cgalign('c','c');
                        end    
                    end       
                    
                    % Keep drawing stim
                    cgpencol(0,0,0); cgrect(-200,0,stimDimensions(1),stimDimensions(2))
                    cgdrawsprite(1,-200,0) % Left stim

                    cgrect(200,0,stimDimensions(1),stimDimensions(2))
                    cgdrawsprite(2,200,0) % Right stim
                    
                    cgrect(-175,originy-25,350,50)
                    cgdrawsprite(3,-175,originy-25)
                    cgrect(175,originy-25,.350,50)
                    cgdrawsprite(4,175,originy-25) 
                    
                    % Draw rating line
                    cgpencol(0,0,0); cgdraw(sx(1),sy(1),sx(2),sy(2)); 
                    cgpencol(1,1,1); cgellipse(sx(1),sy(1),10,10,'f');
                    cgellipse(sx(2),sy(2),10,10,'f');

                    % Update mouse position
                    [x, ~, ~, mp2] = cgmouse; y = originy;
                    if x < sx(1), x = sx(1); end
                    if x > sx(2), x = sx(2); end;
                    rating = (x - sx(1)) / (sx(2) - sx(1)); %0 to 1
                    rating_rt = (time - tstart)/1000; %in s since end of waittime

                    cgpencol(1,1,0); cgellipse(x,y,10,10,'f');
                    cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); %cgtext('+',0,0); 

                    cgflip(data.background);   
                    wait(data.time.outcometime);
                end   
          end     
      
                cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); %cgtext('Outcomes added to total earnings',0,300); 
                cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); 
                cgtext('No preference',0,originy-80);
                cgtext('Strongly confident left',-350, originy-80); %cgtext('|',-200-stimDimensions(1)/2, originy);
                cgtext('Strongly confident right',350, originy-80); %cgtext('|',200+stimDimensions(1)/2,originy);
                    if blocknum == 2
                        cgpencol(0,0,0);
                        cgrect(420,350,200,100);
                        cgpencol(1,1,1)
                        cgfont('Arial',data.FontSizeText)
                        cgtext('Bonus',420,365);
                        if data.wofscore >= 0
                            cgtext('  £',400,330);
                            cgalign('l','c')
                            cgtext(num2str(data.wofscore),428,330)
                            cgalign('c','c');
                        else
                            cgtext('- £',400,330);        
                            cgalign('l','c')
                            cgtext(num2str(abs(data.wofscore)),428,330);
                            cgalign('c','c');
                        end    
                    end
                % Keep drawing stim
                cgpencol(0,0,0); cgrect(-200,0,stimDimensions(1),stimDimensions(2))
                cgdrawsprite(1,-200,0) % Left stim

                cgrect(200,0,stimDimensions(1),stimDimensions(2))
                cgdrawsprite(2,200,0) % Right stim
                
                cgrect(-175,originy-25,350,50)
                cgdrawsprite(3,-175,originy-25)
                cgrect(175,originy-25,.350,50)
                cgdrawsprite(4,175,originy-25)                 

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


        rating = (x - sx(1)) / (sx(2) - sx(1)); %0 to 1
        rating_rt = (time - tstart)/1000; %in s since end of waittime
        
        if rating < .5, choicemade = 1;
        elseif rating > .5, choicemade = 2;
        else choicemade = 0;
        end
     
    rt = (time - tstart2)/1000; 
    experiment_data(itrial,24) = rating;
    experiment_data(itrial,15) = rating_rt; %rt in s %%%%NEW%%%%edited
    experiment_data(itrial,6) = choicemade; %which button pressed (1 left, 2 right) or 0 for none
   
     experiment_data(itrial,7) = experiment_data(itrial,1+choicemade); %%%%NEW%%%%
     
        if experiment_data(itrial,3+choicemade) == 1 %if chose the badmL option %%%%NEW%%%%
            experiment_data(itrial,8) = rewardmL(chosemL); chosemL=chosemL+1; %%%%NEW%%%%
        elseif experiment_data(itrial,3+choicemade) == 2 %if chose the goodmL option
            experiment_data(itrial,8) = rewardmHvL(chosemHvL); chosemHvL = chosemHvL + 1;
        %%%%NEW%%%%
        else %condition of the picture chosen == 3 (goodmH)
            experiment_data(itrial,8) = rewardmHvH(chosemHvH); chosemHvH=chosemHvH+1; %%%%NEW%%%%
        end; %%%%NEW%%%%

     
    %    experiment_data(itrial,8)=data.rewardsize*(experiment_data(itrial,3+choicemade)>rand); %determine outcome %%%%NEW%%%%
    %Add variable in col 23 cumulatively summing rewards
     experiment_data(itrial,23) = data.currentscore + experiment_data(itrial,8);
   
    
    experiment_data(itrial,9) = rt;   %rt in s
    if choicemade==1 
            if blocknum == 2
                cgpencol(0,0,0);
                cgrect(420,350,200,100);
                cgpencol(1,1,1)
                cgfont('Arial',data.FontSizeText)
                cgtext('Bonus',420,365);
             
                if data.wofscore >= 0
                    cgtext('  £',400,330);
                    cgalign('l','c')
                    cgtext(num2str(data.wofscore),428,330)
                    cgalign('c','c');
                else
                    cgtext('- £',400,330);        
                    cgalign('l','c')
                    cgtext(num2str(abs(data.wofscore)),428,330);
                    cgalign('c','c');
                end    
            end
       
        cgpencol(1,.8,0)
        cgrect(-200,0,stimDimensions(1)+20,stimDimensions(2)+20);
        if experiment_data(itrial,3)>0; cgdrawsprite(2,200,0); end %Other stimulus does not disappear
        cgdrawsprite(1,-200,0);
    elseif choicemade==2
           
            if blocknum == 2
                cgpencol(0,0,0);
                cgrect(420,350,200,100);
                cgpencol(1,1,1)
                cgfont('Arial',data.FontSizeText)
                cgtext('Bonus',420,365);
               
                if data.wofscore >= 0
                    cgtext('  £',400,330);
                    cgalign('l','c')
                    cgtext(num2str(data.wofscore),428,330)
                    cgalign('c','c');
                else
                    cgtext('- £',400,330);        
                    cgalign('l','c')
                    cgtext(num2str(abs(data.wofscore)),428,330);
                    cgalign('c','c');
                end    
            end
        cgpencol(1,.8,0)
        cgrect(200,0,stimDimensions(1)+20,stimDimensions(2)+20);
        if experiment_data(itrial,2)>0; cgdrawsprite(1,-200,0); end %Other stimulus does not disappear
        cgdrawsprite(2,200,0)     
    end
   
     cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10);
     cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); 
     cgtext('No preference',0,originy-80);
     cgtext('Strongly confident left',-350, originy-80); %cgtext('|',-200-stimDimensions(1)/2, originy);
     cgtext('Strongly confident right',350, originy-80); %cgtext('|',200+stimDimensions(1)/2,originy);
    cgrect(-175,originy-25,350,50)
    cgdrawsprite(3,-175,originy-25)
    cgrect(175,originy-25,.350,50)
    cgdrawsprite(4,175,originy-25) 
    
     % Draw rating line
     cgpencol(0,0,0); cgdraw(sx(1),sy(1),sx(2),sy(2)); 
     cgpencol(1,1,1); cgellipse(sx(1),sy(1),10,10,'f');
     cgellipse(sx(2),sy(2),10,10,'f');
     cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); %cgtext('+',0,0);
     % Update mouse position
        [x, ~, ~, mp2] = cgmouse; y = originy;
            if x < sx(1), x = sx(1); end
            if x > sx(2), x = sx(2); end;
            cgpencol(1,1,0); cgellipse(x,y,10,10,'f');
            cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); %cgtext('+',0,0);
    
   cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); %cgtext('+',0,0);        
    cgflip(data.background);
    wait(data.time.delay);
    
    data.currentscore = experiment_data(itrial,23);
    
    %%% display the outcome
    %cgpencol(1,1,1); %cgtext('+',0,0);
    %img = data.rewardstim.picture; % Load reward stimuli (dimensions of which will be used for background of loss stimuli too)
    %rewardstimDimensions = data.rewardstim.stimDimensions;
    %cgloadbmp(3,img,rewardstimDimensions(1),rewardstimDimensions(2))
    cgfont('Arial',data.FontSizeTask);
    cgpencol(1,1,1); % White text
        if choicemade==1
             cgpencol(1,.8,0)
             cgrect(-200,0,stimDimensions(1)+20,stimDimensions(2)+20);
            cgdrawsprite(1,-200,0)
            if experiment_data(itrial,3)>0; cgdrawsprite(2,200,0); end %Other stimulus does not disappear
                %cgdrawsprite(3,-200,-100)
                cgpencol(1,1,1)
                cgrect(-200,-100,100,100);
                cgpencol(0,0,0)
                cgtext(num2str(experiment_data(itrial,8)), -200, -100);
                if blocknum == 2
                    cgpencol(0,0,0);
                    cgrect(420,350,200,100);
                    cgpencol(1,1,1)
                    cgfont('Arial',data.FontSizeText)
                    cgtext('Bonus',420,365);
                  
                    if data.wofscore >= 0
                        cgtext('  £',400,330);
                        cgalign('l','c')
                        cgtext(num2str(data.wofscore),428,330)
                        cgalign('c','c');
                    else
                        cgtext('- £',400,330);        
                        cgalign('l','c')
                        cgtext(num2str(abs(data.wofscore)),428,330);
                        cgalign('c','c');
                    end    
                end
        elseif choicemade==2
             cgpencol(1,.8,0)
             cgrect(200,0,stimDimensions(1)+20,stimDimensions(2)+20);
             cgdrawsprite(2,200,0)
            if experiment_data(itrial,2)>0; cgdrawsprite(1,-200,0); end %Other stimulus does not disappear
                %cgdrawsprite(3,200,-100)
                cgpencol(1,1,1)
                cgrect(200,-100,100,100);
                cgpencol(0,0,0)
                cgtext(num2str(experiment_data(itrial,8)), 200, -100);
                
                if blocknum == 2
                    cgpencol(0,0,0);
                    cgrect(420,350,200,100);
                    cgpencol(1,1,1)
                    cgfont('Arial',data.FontSizeText)
                    cgtext('Bonus',420,365);
                    
                   if data.wofscore >= 0
                        cgtext('  £',400,330);
                        cgalign('l','c')
                        cgtext(num2str(data.wofscore),428,330)
                        cgalign('c','c');
                    else
                        cgtext('- £',400,330);        
                        cgalign('l','c')
                        cgtext(num2str(abs(data.wofscore)),428,330);
                        cgalign('c','c');
                   end    
                end
        end
        
    cgrect(-175,originy-25,350,50)
    cgdrawsprite(3,-175,originy-25)
    cgrect(175,originy-25,.350,50)
    cgdrawsprite(4,175,originy-25)  
    
     cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10);
     cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); 
     cgtext('No preference',0,originy-80);
     cgtext('Strongly confident left',-350, originy-80); %cgtext('|',-200-stimDimensions(1)/2, originy);
     cgtext('Strongly confident right',350, originy-80); %cgtext('|',200+stimDimensions(1)/2,originy);

     % Draw rating line
     cgpencol(0,0,0); cgdraw(sx(1),sy(1),sx(2),sy(2)); 
     cgpencol(1,1,1); cgellipse(sx(1),sy(1),10,10,'f');
     cgellipse(sx(2),sy(2),10,10,'f');
     cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask);% cgtext('+',0,0);
     % Update mouse position
        [x, ~, ~, mp2] = cgmouse; y = originy;
            if x < sx(1), x = sx(1); end
            if x > sx(2), x = sx(2); end;
            cgpencol(1,1,0); cgellipse(x,y,10,10,'f');
            cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); %cgtext('+',0,0);
            
    cgflip(data.background); experiment_data(itrial,16) = (time - tstart)/1000; %outcome revealed %%%%NEW%%%%edited
    wait(data.time.outcometime);

        
        if blocknum == 2
            cgpencol(0,0,0);
            cgrect(420,350,200,100);
            cgpencol(1,1,1)
            cgfont('Arial',data.FontSizeText)
            cgtext('Bonus',420,365);    
            if data.wofscore >= 0
                cgtext('  £',400,330);
                cgalign('l','c')
                cgtext(num2str(data.wofscore),428,330)
                cgalign('c','c');
            else
                cgtext('- £',400,330);        
                cgalign('l','c')
                cgtext(num2str(abs(data.wofscore)),428,330);
                cgalign('c','c');
            end    
        end
    cgrect(-175,originy-25,350,50)
    cgdrawsprite(3,-175,originy-25)
    cgrect(175,originy-25,.350,50)
    cgdrawsprite(4,175,originy-25)    
    
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); %cgtext('Outcomes added to total earnings',0,300);
    cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); 
    cgtext('No preference',0,originy-80);
    cgtext('Strongly confident left',-350, originy-80); %cgtext('|',-200-stimDimensions(1)/2, originy);
    cgtext('Strongly confident right',350, originy-80); %cgtext('|',200+stimDimensions(1)/2,originy);

    % Draw rating line
    cgpencol(0,0,0); cgdraw(sx(1),sy(1),sx(2),sy(2)); 
    cgpencol(1,1,1); cgellipse(sx(1),sy(1),10,10,'f');
    cgellipse(sx(2),sy(2),10,10,'f');
    
    cgflip(data.background); experiment_data(itrial,17) = (time - tstart)/1000; %iti starts %%%%NEW%%%%edited
    wait(data.time.iti+round((data.time.maxiti-data.time.iti)*rand)); %jitter between min and max iti
    experiment_data(itrial,18) = (time - tstart)/1000; %trials end - could jitter %%%%NEW%%%%edited
    
     if itrial == data.ntrials
        cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10);
        cgtext([' Good job! '],0,100);
        cgtext(['You have reached the exit of the this garden'],0,0);
        cgtext('Click to continue',0,-200);
        cgflip(data.background);
        wait(data.time.instructiontime);
        [~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end
        [experiment_data(itrial,10), experiment_data(itrial,11), experiment_data(itrial,12), alltime] = happy_rating('nowhappy',data.time.ratewait,data); %respkeys, waittime, resptime, startlocation, pixperkey);
        cgpencol(1,1,1); %cgtext('+', 0, 0); 
        cgflip(data.background);
        experiment_data(itrial,19) = (alltime(1) - tstart)/1000; %happy question onset time
        experiment_data(itrial,20) = (alltime(2) - tstart)/1000; %happy line onset
        experiment_data(itrial,21) = experiment_data(itrial,12) + experiment_data(itrial,19); %happy rating made
        experiment_data(itrial,22) = (time - tstart)/1000;       %ITI (0.7s) and 1s of looking at rating already built in
        data.maintask{blocknum} = experiment_data;
        save(data.savefile, 'data');    
     end
    

end
wait(data.time.iti);


