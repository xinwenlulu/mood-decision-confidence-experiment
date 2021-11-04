
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%practice block
function [experiment_data]   = Xin_PracticeBlock(data)

picprob = [1 2 1 2];
trialmtx = ones(data.npracticetrials,4);

for n = 1:data.npracticetrials
    trialmtx(n,:) = picprob;
end

for n = 1:data.npracticetrials
    if rand>0.5
        trialmtx(n,:) = trialmtx(n,[2 1 4 3]);
    end
end

tstart = time;  
 
wait(data.time.iti);

experiment_data = zeros(data.npracticetrials,24);
experiment_data(:,1) = 1:data.npracticetrials;
experiment_data(:,13) = 0; %wheel of fortune outcome %%%%NEW%%%%
experiment_data(:,2:5) = trialmtx; %zero anything that is empty - 0 stimulus and 0 prob
%%%%%%%%%%%%%%%%%%%%%%%%
bad = [10 20 25 30 40];
good = [45 60 70 80 95];
% shuffle the rewards 
bad = bad(randperm(length(bad)));
good = good(randperm(length(good)));

%predetermine rewards for each option
rewardmL = repmat(bad,1,ceil(data.npracticetrials/size(bad,2)));
rewardmH = repmat(good,1,ceil(data.npracticetrials/size(good,2)));
chosemL = 1; chosemH = 1; 

for itrial = 1:data.npracticetrials
    disp(experiment_data(:,1:14));
    allowleftclick = 1; %allowrightclick = -1; 
    stimDimensions = data.stimDimensions;
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); %cgtext('+',0,0);
    
    pic = data.practice_pic{experiment_data(itrial,2)};
    cgloadbmp(1,pic,stimDimensions(1),stimDimensions(2));
 
    cgrect(-200,0,stimDimensions(1),stimDimensions(2));
    cgdrawsprite(1,-200,0);
    
    pic = data.practice_pic{experiment_data(itrial,3)};
    cgloadbmp(2,pic,stimDimensions(1),stimDimensions(2));
  
    cgrect(200,0,stimDimensions(1),stimDimensions(2));
    cgdrawsprite(2,200,0)
    
    originy = -180; 
    pic = data.confidence{1};
    cgloadbmp(3,pic,350,50);
    cgrect(-175,originy-25,350,50);
    cgdrawsprite(3,-175,originy-25);
    pic = data.confidence{2};
    cgloadbmp(4,pic,350,50);
    cgrect(175,originy-25,350,50);
    cgdrawsprite(4,175,originy-25);
    
    cgflip(data.background);
    tstart2 = time; experiment_data(itrial,14) = (tstart2 - tstart)/1000;
    originy = -180; sx = [-350 350]; sy = [originy, originy]; % for lefthandside
    timelimit=Inf;
    %randstart = 0.5;               %skip the randstart and just start from the middle
    mp2 = 0; 
    x = -280+560*rand;
    y = 0; %initialize mouse to center of screen at moment line appears
    cgmouse(x, 0); %mouse jumps to random location on line
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
  
                cgflip(data.background);        
             end

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
    elseif experiment_data(itrial,3+choicemade) == 2 %if chose the goodmH option
        experiment_data(itrial,8) = rewardmH(chosemH); chosemH = chosemH + 1;
    %%%%NEW%%%%
    else %condition of the picture chosen == 3 (goodmH)
        experiment_data(itrial,8) = 0; %%%%NEW%%%%
    end %%%%NEW%%%%
     
    %Add variable in col 23 cumulatively summing rewards
    experiment_data(itrial,23) = data.currentscore + experiment_data(itrial,8);
    experiment_data(itrial,9) = rt;   %rt in s
    
    if choicemade==1 
       
        cgpencol(1,.8,0)
        cgrect(-200,0,stimDimensions(1)+20,stimDimensions(2)+20);
        if experiment_data(itrial,3)>0; cgdrawsprite(2,200,0); end %Other stimulus does not disappear
        cgdrawsprite(1,-200,0);
    elseif choicemade==2
        
        cgpencol(1,.8,0)
        cgrect(200,0,stimDimensions(1)+20,stimDimensions(2)+20);
        if experiment_data(itrial,2)>0; cgdrawsprite(1,-200,0); end %Other stimulus does not disappear
        cgdrawsprite(2,200,0)     
    end
    
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); %cgtext('Outcomes added to total earnings',0,300); 
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
    cgfont('Arial',data.FontSizeTask);
    cgpencol(1,1,1); % White text
        if choicemade==1
           cgpencol(1,.8,0)
           cgrect(-200,0,stimDimensions(1)+20,stimDimensions(2)+20);
           cgdrawsprite(1,-200,0)
           if experiment_data(itrial,3)>0; cgdrawsprite(2,200,0); end %Other stimulus does not disappear
            cgpencol(1,1,1)
            cgrect(-200,-100,100,100);
            cgpencol(0,0,0)
            cgtext(num2str(experiment_data(itrial,8)), -200, -100);
            
        elseif choicemade==2
            cgpencol(1,.8,0)
            cgrect(200,0,stimDimensions(1)+20,stimDimensions(2)+20);
            cgdrawsprite(2,200,0)
            if experiment_data(itrial,2)>0; cgdrawsprite(1,-200,0); end %Other stimulus does not disappear
            cgpencol(1,1,1)
            cgrect(200,-100,100,100);
            cgpencol(0,0,0)
            cgtext(num2str(experiment_data(itrial,8)), 200, -100);
           
        end
        
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); %cgtext('Outcomes added to total earnings',0,300); 
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
     cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask);% cgtext('+',0,0);
     % Update mouse position
     [x, ~, ~, mp2] = cgmouse; y = originy;
    if x < sx(1), x = sx(1); end
    if x > sx(2), x = sx(2); end;
    cgpencol(1,1,0); cgellipse(x,y,10,10,'f');
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeTask); %cgtext('+',0,0);

    cgflip(data.background); 
    experiment_data(itrial,16) = (time - tstart)/1000; %outcome revealed %%%%NEW%%%%edited
    wait(data.time.outcometime);
    
    
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); %cgtext('Outcomes added to total earnings',0,300); 
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
    cgflip(data.background); experiment_data(itrial,17) = (time - tstart)/1000; %iti starts %%%%NEW%%%%edited
    wait(data.time.iti+round((data.time.maxiti-data.time.iti)*rand)); %jitter between min and max iti
    experiment_data(itrial,18) = (time - tstart)/1000; %trials end - could jitter %%%%NEW%%%%edited
end
wait(data.time.iti);

% Variance
% tstart = time;
% pics = picprob(1:2);
% conditions = picprob(3:4);
% index = randperm(2);
% pics = pics(index);
% conditions = conditions(index);
% variance = nan(2,5);

% for i = 1:2
%     allowleftclick = 1; %allowrightclick = -1; 
%     pic = data.practice_pic{pics(i)};
%     cgloadbmp(1,pic,stimDimensions(1),stimDimensions(2));
%     cgrect(0,0,stimDimensions(1),stimDimensions(2));
%     cgdrawsprite(1,0,0);
%     
%     originy = -180; 
%     sx = [-350 350]; sy = [originy, originy];
%     tstart3 = time;
%     variance(i,5) = (tstart3 - tstart)/1000;
%     timelimit=Inf;
%     mp2 = 0; x = -280+560*rand; y = 0; %initialize mouse to center of screen at moment line appears
%     cgmouse(x, 0);
%        while (mp2 ~= 1) && time < tstart + timelimit
%             cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10);
%             cgtext('How predictable do you think this flower is?',0,300)
%             cgpencol(1,1,1); cgfont('Arial',data.FontSizeText);
%             cgtext('(in terms of how consistent the number of gems you can get from this flower)',0,250)
%             %cgtext('Outcomes added to total earnings',0,300); 
%             cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); 
%             %cgtext('No preference',0,originy-30);
%             cgtext('Not at all predictable',-350, originy-30); %cgtext('|',-200-stimDimensions(1)/2, originy);
%             cgtext('Very predictable',350, originy-30); %cgtext('|',200+stimDimensions(1)/2,originy);
% 
%             % Keep drawing stim
%             cgpencol(0,0,0); cgrect(0,0,stimDimensions(1),stimDimensions(2))
%             cgdrawsprite(1,0,0) % Left stim
%      
%             % Draw rating line
%             cgpencol(0,0,0); cgdraw(sx(1),sy(1),sx(2),sy(2)); 
%             cgpencol(1,1,1); cgellipse(sx(1),sy(1),10,10,'f');
%             cgellipse(sx(2),sy(2),10,10,'f');
% 
%             % Update mouse position
%             [x, ~, ~, mp2] = cgmouse; y = originy;
%             if x < sx(1), x = sx(1); end
%             if x > sx(2), x = sx(2); end;
%             cgpencol(1,1,0); cgellipse(x,y,10,10,'f');
% 
%             cgflip(data.background);        
%          end
%             cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); %cgtext('Outcomes added to total earnings',0,300); 
%             cgtext('How predictable do you think this flower is?',0,300)
%             cgpencol(1,1,1); cgfont('Arial',data.FontSizeText);
%             cgtext('(in terms of how consistent the number of gems you can get from this flower)',0,250)
%             cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); 
%             %cgtext('No preference',0,originy-30);
%             cgtext('Not at all predictable',-350, originy-30); %cgtext('|',-200-stimDimensions(1)/2, originy);
%             cgtext('Very predictable',350, originy-30); %cgtext('|',200+stimDimensions(1)/2,originy);
% 
%             % Keep drawing stim
%             cgpencol(0,0,0); cgrect(0,0,stimDimensions(1),stimDimensions(2))
%             cgdrawsprite(1,0,0) % Left stim
%      
%             % Draw rating line
%             cgpencol(0,0,0); cgdraw(sx(1),sy(1),sx(2),sy(2)); 
%             cgpencol(1,1,1); cgellipse(sx(1),sy(1),10,10,'f');
%             cgellipse(sx(2),sy(2),10,10,'f');
% 
%             % Update mouse position
%             [x, ~, ~, mp2] = cgmouse; y = originy;
%             if x < sx(1), x = sx(1); end
%             if x > sx(2), x = sx(2); end;
%             cgpencol(1,1,0); cgellipse(x,y,10,10,'f'); 
%             cgflip(data.background);
%         rating = (x - sx(1)) / (sx(2) - sx(1)); %0 to 1
%         rt = (time - tstart3)/1000; %in s since end of waittime
%      
%     variance(i,1) = pics(i);    
%     variance(i,2) = conditions(i);
%     variance(i,3)= rating;
%     variance(i,4) = rt; %rt in s %%%%NEW%%%%edited%%%%NEW%%%%
%     wait(data.time.delay)
%     
% % keep drawing everything apart from stim
%     cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); %cgtext('Outcomes added to total earnings',0,300); 
%     cgtext('How predictable do you think this flower is?',0,300)
%     cgpencol(1,1,1); cgfont('Arial',data.FontSizeText);
%     cgtext('(in terms of how consistent the number of gems you can get from this flower)',0,250)
%     cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); 
%     %cgtext('No preference',0,originy-30);
%     cgtext('Not at all predictable',-350, originy-30); %cgtext('|',-200-stimDimensions(1)/2, originy);
%     cgtext('Very predictable',350, originy-30); %cgtext('|',200+stimDimensions(1)/2,originy);
% 
%     % Draw rating line
%     cgpencol(0,0,0); cgdraw(sx(1),sy(1),sx(2),sy(2)); 
%     cgpencol(1,1,1); cgellipse(sx(1),sy(1),10,10,'f');
%     cgellipse(sx(2),sy(2),10,10,'f');
% 
%     % Update mouse position
%     [x, ~, ~, mp2] = cgmouse; y = originy;
%     if x < sx(1), x = sx(1); end
%     if x > sx(2), x = sx(2); end;
%     cgpencol(1,1,0); cgellipse(x,y,10,10,'f'); 
%     cgflip(data.background);
%     wait(data.time.iti);
%     
%     
% end
wait(data.time.iti);