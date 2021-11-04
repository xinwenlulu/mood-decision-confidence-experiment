function [rating, randstart, rt, alltime] = confidence_rating_sascha(ratingtype, waittime,data)

%[rating, randstart, rt] = happy_rating(ratingtype)
%if ratingtype unused, asks how happy you are at the moment with no wait
%if ratingtype is 'lifehappy', asks how happy you are +wait
%if ratingtype is 'nowhappy', asks how happy you are +wait
%if ratingtype is 'happy', ask how happy you are +wait
%before your first trial for money +wait
%waits for WAITTIME milliseconds until asking for rating

timelimit = Inf;

%--------------------------------------------------
sx = [-384 384]; sy = [0 0]; %in pixels, since it's 1024 wide, this is 75% of the screen

%random start with the mouse means occasionally they will have to pick up
%the mouse if they consistently rate above or below the midpoint because on
%average the mouse will end up in the center of the line between trials.

%randstart = rand;             %start in a random location on the line (0,1)
randstart = 0.5;               %skip the randstart and just start from the middle
rating = nan; rt = nan;        %in case early exit

cgfont('Arial',data.FontSizeText);
cgpencol(1,1,1);
if exist('ratingtype'),
    if strcmp(ratingtype,'lifeconfidence'),
        cgtext('Taken all together, how confident do you feel these days?',0,300);
        cgtext('Mark your rating relative to the least and most confident time of your life.',0,250);
        cgflip(data.background); 
        wait(waittime);
        cgtext('Taken all together, how confident do you feel these days?',0,300);
        cgtext('Mark your rating relative to the least and most confident time of your life.',0,250);
        cgtext('Click to give your rating',0,50);
        cgflip(data.background); 

        [~,~,~,mp]=cgmouse; mp=0; %no mouse press yet
        while mp~=1, %1 from leftmouse, 4 for rightmouse - if a free choice trial
            [~,~,~, mp] = cgmouse;
        end;
       
    elseif strcmp(ratingtype,'nowconfidence'),
        %cgtext('Think about right now. How happy are you at this moment?',0,300);
        %cgflip(data.background); 
        %wait(waittime);
        cgtext('Think about right now. How confident do you feel at this moment?',0,300);
        cgtext('Click to give your rating',0,50);
        cgflip(data.background); 
        
        [~,~,~,mp]=cgmouse; mp=0; %no mouse press yet
        while mp~=1, %1 from leftmouse, 4 for rightmouse - if a free choice trial
            [~,~,~, mp] = cgmouse;
        end;
        
    elseif strcmp(ratingtype,'choiceconfidence'),
        cgtext('How confident are you that you know which animal has the most gems?',0,100);
        cgflip(data.background); 
        wait(waittime);
    end;
    
    tstart = time;                %get time in ms from first time call
    
    tstart2 = time;
    alltime = [tstart tstart2];
else
    cgtext('How confident are you that you know which animal has the most gems?',0,100);
    cgflip(data.background);
    tstart2 = time;               %get time in ms from first time call
    alltime = [tstart2 tstart2];
end;


mp = 0; x = 0; y = 0; %initialize mouse to center of screen at moment line appears
cgmouse((sx(2)-sx(1))*randstart + sx(1), 0); %mouse jumps to random location on line
while (mp ~= 1) && time < tstart + timelimit;
    cgpencol(1,1,1);
    if exist('ratingtype'),
        if strcmp(ratingtype,'lifeconfidence'),
            cgtext('Taken all together, how confident do you feel these days?',0,300);
            cgtext('Mark your rating relative to the least and most confident time of your life.',0,250);
        elseif strcmp(ratingtype,'nowconfidence');
            cgtext('Think about right now. How confident do you feel at this moment?',0,300);
        elseif strcmp(ratingtype,'choiceconfidence'),
            cgtext('How confident are you that you know which animal has the most gems?',0,100);
        end;
        if sum(strfind(ratingtype,'choiceconfidence')), %not money
            cgtext('not at all',sx(1)-65,sy(1)+20); cgtext('confident',sx(1)-65,sy(1)-20);
            cgtext('very',sx(2)+50,sy(2)+20); cgtext('confident',sx(2)+50,sy(2)-20);
        end;
    else
        cgtext('How confident are you that you know which animal has the most gems?',0,100);
    end;
    
    cgpencol(0,0,0); cgdraw(sx(1),sy(1),sx(2),sy(2)); %draw rating line
    cgpencol(1,1,1); cgellipse(sx(1),sy(1),10,10,'f');
    cgellipse(sx(2),sy(2),10,10,'f');
    
    % update mouse position
    [x, ~, ~, mp] = cgmouse; y = 0;
    if x < sx(1), x = sx(1); end
    if x > sx(2), x = sx(2); end;
    cgpencol(1,1,0); cgellipse(x,y,10,10,'f');
    cgflip(data.background);
end
rating = (x - sx(1)) / (sx(2) - sx(1)); %0 to 1
rt = (time - tstart2)/1000; %in s since end of waittime

%%%%%%%%%%%%%
wait(700); %show rating made for 1s
cgfont('Arial',data.FontSizeTask);
cgpencol(1,1,1); cgtext('+', 0, 0);
cgflip(data.background);
wait(700); %show fixation cross for 700ms
