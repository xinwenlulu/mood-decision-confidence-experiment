function [variance]   = Xin_VarianceEstimate(data,blocknum)
tstart = time;
stimDimensions = data.stimDimensions;
variance = nan(3,5);
pics = data.picprob(blocknum,1:3);
conditions = data.picprob(blocknum,4:6);
index = randperm(3);
pics = pics(index);
conditions = conditions(index);

for i = 1:3
    allowleftclick = 1; %allowrightclick = -1; 
    pic = data.pictures{pics(i)};
    cgloadbmp(1,pic,stimDimensions(1),stimDimensions(2));
    cgrect(0,0,stimDimensions(1),stimDimensions(2));
    cgdrawsprite(1,0,0);
    
    originy = -180; 
    sx = [-350 350]; sy = [originy, originy];
    tstart3 = time;
    variance(i,5) = (tstart3 - tstart)/1000;
    timelimit=Inf;
    mp2 = 0; x = -280+560*rand; y = 0; %initialize mouse to center of screen at moment line appears
    cgmouse(x, 0);
       while (mp2 ~= 1) && time < tstart + timelimit
            cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10);
            cgtext('how much did the number of gems vary for this flower?',0,250)
            %cgtext('Outcomes added to total earnings',0,300); 
            cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); 
            %cgtext('No preference',0,originy-30);
            cgtext('It did not vary at all',-350, originy-30); %cgtext('|',-200-stimDimensions(1)/2, originy);
            cgtext('It varied a lot',350, originy-30); %cgtext('|',200+stimDimensions(1)/2,originy);

            % Keep drawing stim
            cgpencol(0,0,0); cgrect(0,0,stimDimensions(1),stimDimensions(2))
            cgdrawsprite(1,0,0) % Left stim
     
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
            cgtext('how much did the number of gems vary for this flower?',0,250)
            cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); 
            %cgtext('No preference',0,originy-30);
            cgtext('It did not vary at all',-350, originy-30); %cgtext('|',-200-stimDimensions(1)/2, originy);
            cgtext('It varied a lot',350, originy-30); %cgtext('|',200+stimDimensions(1)/2,originy);

            % Keep drawing stim
            cgpencol(0,0,0); cgrect(0,0,stimDimensions(1),stimDimensions(2))
            cgdrawsprite(1,0,0) % Left stim
     
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
        rating = (x - sx(1)) / (sx(2) - sx(1)); %0 to 1
        rt = (time - tstart3)/1000; %in s since end of waittime
     
    variance(i,1) = pics(i);    
    variance(i,2) = conditions(i);
    variance(i,3)= rating;
    variance(i,4) = rt; %rt in s %%%%NEW%%%%edited%%%%NEW%%%%
    wait(data.time.delay)
    
% keep drawing everything apart from stim
    cgpencol(1,1,1); cgfont('Arial',data.FontSizeText+10); %cgtext('Outcomes added to total earnings',0,300); 
    cgtext('how much did the number of gems vary for this flower?',0,250)
    cgfont('Arial',data.FontSizeText); cgtext('|',0, originy); 
    %cgtext('No preference',0,originy-30);
    cgtext('It did not vary at all',-350, originy-30); %cgtext('|',-200-stimDimensions(1)/2, originy);
    cgtext('It varied a lot',350, originy-30); %cgtext('|',200+stimDimensions(1)/2,originy);

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
    wait(data.time.iti);
    
    
end
wait(data.time.iti);