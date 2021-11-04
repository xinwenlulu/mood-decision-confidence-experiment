function [allcomparisons, rt, alltime] = testblock(data, nBlocks)

%data.picprob(blocknum,:)
timelimit = Inf;


% allcomparisons = zeros(8:data.ntrials);
% allcomparisons(:,1) = (1:8);
% allcomparisons(:,13) = wof; 
% 


%if after block 1 then, 6 comparisons - needs to be general enough to
%handle 3rd block, 4th blocks looking at past 3 blocks
% 
% trialmtx = [data.picprob(1:2,:); data.picprob
% 
% 11 vs 12
% 11 vs 21
% 21 vs 22
% 12 vs 22

%hard code trialmtx for 4 comparisons
% testmtx = picprob(randperm(size(picprob,1)),:); %shuffle order
% for n=1:size(testmtx,1), if rand>0.5, testmtx(n,:) = testmtx(n,[2 1 4 3]); end; end; %left-right shuffle
% 
% for i=1:nBlocks
%     for j=1:size(picprob,1);
%         jtest(i,j) = data.picprob(i,j);
%     end
% end

% % % Liam
% % % This is not correct because the second testblock just refers to the
% % % 2nd set of stimuli being tested -- NOT necessarily 75/25 split
% % disp('overwriting probabilities in trialmtx to 0.25 0.75')
% % data.picprob(:,3) = [0.25 0.25 0.25 0.25];
% % data.picprob(:,4) = [0.75 0.75 0.75 0.75];


% Good vs bad stim, B1 and B2
allcomparisons{1,1} = char( data.pictures( data.picprob(3,1) ) );
allcomparisons{1,2} = char( data.pictures( data.picprob(3,2) ) );
allcomparisons{1,3} = [data.picprob(3,3), data.picprob(3,4)] ;
allcomparisons{2,1} = char( data.pictures( data.picprob(4,1) ) );
allcomparisons{2,2} = char( data.pictures( data.picprob(4,2) ) );
allcomparisons{2,3} = [data.picprob(4,3), data.picprob(4,4)] ; 

% Stim matched for probability, B1 and B2
allcomparisons{3,1} = char( data.pictures( data.picprob(3,1) ));
allcomparisons{3,2} = char( data.pictures( data.picprob(4,1) ));
allcomparisons{3,3} = [data.picprob(3,3), data.picprob(4,3)];
allcomparisons{4,1} = char( data.pictures( data.picprob(3,2) ));
allcomparisons{4,2} = char( data.pictures( data.picprob(4,2) ));
allcomparisons{4,3} = [data.picprob(3,4), data.picprob(4,4)]; % Previously 1,3 and 2,3 (duplicated)
%allcomparisons{8,3} = [data.picprob(1,3), data.picprob(2,3)];


% allcomparisons{3,1} = data.b1.stim{1};
% allcomparisons{3,2} = data.b2.stim{1};
% allcomparisons{3,3} = [data.b1.stim{1,2}, data.b2.stim{1,2}];
% allcomparisons{4,1} = data.b1.stim{2};
% allcomparisons{4,2} = data.b2.stim{2};
% allcomparisons{4,3} = [data.b1.stim{2,2}, data.b2.stim{2,2}];
data.colorleft = [0,0,0];
data.colorright = [0,0,0];

if nBlocks == 3
%     allcomparisons = NaN(24,15)
%     allcomparisons(:,1) = 1:24
 
    allcomparisons{5,1} = char( data.pictures( data.picprob(3,1) ));
    allcomparisons{5,2} = char( data.pictures( data.picprob(3,2) ));
    allcomparisons{5,3} = [data.picprob(3,3), data.picprob(3,4)];
    
    allcomparisons{6,1} = char( data.pictures( data.picprob(3,1) ));
    allcomparisons{6,2} = char( data.pictures( data.picprob(1,1) ));
    allcomparisons{6,3} = [data.picprob(3,3), data.picprob(1,3)];

    allcomparisons{7,1} = char( data.pictures( data.picprob(3,1) ));
    allcomparisons{7,2} = char( data.pictures( data.picprob(1,2) ));
    allcomparisons{7,3} = [data.picprob(3,3), data.picprob(1,4)];
    
    allcomparisons{8,1} = char( data.pictures( data.picprob(3,1) ));
    allcomparisons{8,2} = char( data.pictures( data.picprob(2,1) ));
    allcomparisons{8,3} = [data.picprob(3,3), data.picprob(2,3)];

    allcomparisons{9,1} = char( data.pictures( data.picprob(3,2) ));
    allcomparisons{9,2} = char( data.pictures( data.picprob(1,1) ));
    allcomparisons{9,3} = [data.picprob(3,4), data.picprob(1,3)];
    
    allcomparisons{10,1} = char( data.pictures( data.picprob(3,2) ));
    allcomparisons{10,2} = char( data.pictures( data.picprob(1,2) ));
    allcomparisons{10,3} = [data.picprob(3,4), data.picprob(1,4)];

    allcomparisons{11,1} = char( data.pictures( data.picprob(3,2) ));
    allcomparisons{11,2} = char( data.pictures( data.picprob(2,1) ));
    allcomparisons{11,3} = [data.picprob(3,4), data.picprob(2,3)];
    
    allcomparisons{12,1} = char( data.pictures( data.picprob(3,2) ));
    allcomparisons{12,2} = char( data.pictures( data.picprob(2,2) ));
    allcomparisons{12,3} = [data.picprob(3,4), data.picprob(2,4)];
    
    
%     allcomparisons{5,1} = data.b1.stim{1};
%     allcomparisons{5,2} = data.b3.stim{1};
%     allcomparisons{5,3} = [data.b1.stim{1,2}, data.b3.stim{1,2}];
%     allcomparisons{6,1} = data.b1.stim{2};
%     allcomparisons{6,2} = data.b3.stim{2};
%     allcomparisons{6,3} = [data.b1.stim{2,2}, data.b3.stim{2,2}];
% 
%     allcomparisons{7,1} = data.b2.stim{1};
%     allcomparisons{7,2} = data.b3.stim{1};
%     allcomparisons{7,3} = [data.b2.stim{1,2}, data.b3.stim{1,2}];
%     allcomparisons{8,1} = data.b2.stim{2};
%     allcomparisons{8,2} = data.b3.stim{2};
%     allcomparisons{8,3} = [data.b2.stim{2,2}, data.b3.stim{2,2}];

end

compareIndex = randperm(size(allcomparisons,1));

if ~strcmp( data.initials, 'debugging' )
    cgfont('Arial',data.FontSizeText); cgpencol(1,1,1); 
    cgflip(data.background)
    cgtext('Choose the animal you think will give you the most gems',0,100);
    cgtext('You will find out your total earnings at the end',0,0);
    cgtext('Click to continue',0,-200);
    cgflip(data.background)

    [~,~,~,mp]=cgmouse; mp=0; %no mouse press yet
    while (mp ~=1) && (mp ~=4), %1 from leftmouse, 4 for rightmouse, 2 for middle mouse
        [~,~,~, mp] = cgmouse;
    end;

%     cgtext('In this block choose between the animals you have seen so far',0,100);
%     cgtext('',0,0);
%     
%     cgflip(data.background);
%     
%      [~,~,~,mp]=cgmouse; mp=0; %no mouse press yet
%     while (mp ~=1) && (mp ~=4), %1 from leftmouse, 4 for rightmouse, 2 for middle mouse
%         [~,~,~, mp] = cgmouse;
%     end;
    
    cgfont('Arial',data.FontSizeTask);
    cgtext('+',0,0);
    cgflip(data.background);
   wait(data.time.delay);
    
%     cgfont('Arial',data.FontSizeText);
%     cgtext('Now choose between the animals',0,0);
%     cgtext('You only get feedback on if you have won at the end',0,0);
%     cgtext('You only get feedback on if you have won at the end',0,0);
   % cgflip(data.background);
    
end

for i=1:size(allcomparisons,1)
    i
    disp( allcomparisons{compareIndex(i),1} ) % Print the comparison to console
    disp( allcomparisons{compareIndex(i),2} ) % Print the comparison to console
    disp( allcomparisons{compareIndex(i),3} ) % Print the comparison to console
    stimDimensions = data.stimDimensions;
    
    cgfont('Arial',data.FontSizeText+10);
    cgtext('Outcomes added to total earnings',0,300);
        
    % Draw left stim
    leftstim = char( allcomparisons{compareIndex(i),1} );
    cgloadbmp(1,leftstim,stimDimensions(1),stimDimensions(2))
    cgpencol(data.colorleft)
    cgrect(-200,0,stimDimensions(1)+10,stimDimensions(2)+10)
    cgdrawsprite(1,-200,0)
    
    
    % Draw right stim
    rightstim = char( allcomparisons{compareIndex(i),2} );
    cgpencol(data.colorright)
    cgrect(200,0,stimDimensions(1)+10,stimDimensions(2)+10)
    cgloadbmp(2,rightstim,stimDimensions(1),stimDimensions(2))
    cgdrawsprite(2,200,0)
    cgpencol(1,1,1)

    cgfont('Arial',data.FontSizeTask);
    cgpencol(1,1,1);
    cgtext('+',0,0);
    
  
    tstart = time;                %get time in ms from first time call
    %wait(waittime);
    tstart2 = time;
    alltime = [tstart tstart2];
    cgflip(data.background);
    cgfont('Arial',data.FontSizeTask);
    [~,~,~,mp]=cgmouse; mp=0; %no mouse press yet
    
    while (mp ~=1) && (mp ~=4), %1 from leftmouse, 4 for rightmouse, 2 for middle mouse
        [~,~,~, mp] = cgmouse;
    end;
    rt = (time - tstart2)/1000; %rt in s
    %choicemade = sqrt(mp); %1 for left, 2 for right
    
%     %% Get line measurement too
%     originy = 220;
%     
%     %sx = [-384 384]; sy = [originy, originy]; %in pixels, since it's 1024 wide, this is 75% of the screen
% if mp==1
%     sx = [-200 0]; sy = [originy, originy]; % for lefthandside
% elseif mp==4
%     sx = [200 400]; sy = [originy, originy]; % for righthandside
% end
% 
%     
%     %randstart = rand;             %start in a random location on the line (0,1)
%     randstart = 0.5;               %skip the randstart and just start from the middle
%     rating = nan; rt = nan;        %in case early exit
% 
%     %cgfont('Arial',data.FontSizeText);
%     cgpencol(1,1,1);
%     tstart = time;                %get time in ms from first time call
%     tstart2 = time;
%     alltime = [tstart tstart2];
%     
%     %mp2 = 0; x = 220; y = 0; %initialize mouse to center of screen at moment line appears
%     mp2 = 0; x = sx(1); y = 0; %initialize mouse to center of screen at moment line appears
%     %cgmouse((sx(2)-sx(1))*randstart + sx(1), 0); %mouse jumps to random location on line
%     cgmouse(sx(1), 0); %mouse jumps to random location on line
%     
%      while (mp2 ~= 1) && time < tstart + timelimit;
% %         
% %         % Draw left stim
% %     leftstim = allcomparisons{compareIndex(i),1};
% %     cgloadbmp(1,leftstim,stimDimensions(1),stimDimensions(2))
% %     cgpencol(data.colorleft)
% %     cgrect(-200,0,stimDimensions(1)+10,stimDimensions(2)+10)
% %     cgdrawsprite(1,-200,0)
% %     
% %     
% %     % Draw right stim
% %     rightstim = allcomparisons{compareIndex(i),2};
% %     cgpencol(data.colorright)
% %     cgrect(200,0,stimDimensions(1)+10,stimDimensions(2)+10)
% %     cgloadbmp(1,rightstim,stimDimensions(1),stimDimensions(2))
% %     cgdrawsprite(1,200,0)
% %     cgpencol(1,1,1)
%     if mp==1
%         cgdrawsprite(1,-200,0)
%     elseif mp==4
%         cgdrawsprite(2,200,0)
%     end  
%         
%         cgpencol(1,1,1);
%         cgfont('Arial',data.FontSizeTask);
%         cgtext('How much do you prefer it?',0,300);
%         cgpencol(0,0,0); cgdraw(sx(1),sy(1),sx(2),sy(2)); %draw rating line
%         cgpencol(1,1,1); cgellipse(sx(1),sy(1),10,10,'f');
%         cgellipse(sx(2),sy(2),10,10,'f');
% 
%         % update mouse position
%         [x, ~, ~, mp2] = cgmouse; y = originy;
%         if x < sx(1), x = sx(1); end
%         if x > sx(2), x = sx(2); end;
%         cgpencol(1,1,0); cgellipse(x,y,10,10,'f');
%         cgflip(data.background);
%     end
%     rating = (x - sx(1)) / (sx(2) - sx(1)); %0 to 1
%     rt = (time - tstart2)/1000; %in s since end of waittime
% 
%     allcomparisons{compareIndex(i),4} = mp; % log preference choice and strength
%     allcomparisons{compareIndex(i),5} = rating; % log preference choice and strength
%     
    %wait(700); %show rating made 
        
    %%
    
    
    %preferredstim=mp;
        
    if mp==1
        cgloadbmp(1,leftstim,stimDimensions(1),stimDimensions(2))
        cgpencol(1,.8,0)
        cgrect(-200,0,stimDimensions(1)+20,stimDimensions(2)+20);
        cgdrawsprite(1,-200,0)
    elseif mp==4
        cgpencol(1,.8,0)
        cgrect(200,0,stimDimensions(1)+20,stimDimensions(2)+20)
        cgloadbmp(1,rightstim,stimDimensions(1),stimDimensions(2))
        cgdrawsprite(1,200,0)
        cgpencol(1,1,1)
    end
    
    cgfont('Arial',data.FontSizeTask);
    cgpencol(1,1,1);
    cgtext('+',0,0);
    %cgtext('Which is best?',0,300);
cgfont('Arial',data.FontSizeText+10);
    cgtext('Outcomes added to total earnings',0,300);
    cgflip(data.background); 
    
%     for iLever=1:length(data.fLevers)
% 
%             if iLever<10
%                imagenumber=['0', num2str(iLever)];
%             else
%                 imagenumber=num2str(iLever);
%             end
% 
%         cgpencol(1,1,1);
%         cgtext('+',0,0);
%         %cgloadbmp(1,data.pictureleft,294,154)
%         if mp==1
%             img = strrep( fullfile(data.imgfolder, leftstim), '.bmp', ['-', imagenumber, '.bmp'] );
%             cgloadbmp(1,img,stimDimensions(1),stimDimensions(2))
%             cgdrawsprite(1,-200,0)%cgdrawsprite(1,-200,0)
%             cgpencol(data.colorleft)
%             cgrect(-200,0,stimDimensions(1)+10,stimDimensions(2)+10)
%             cgdrawsprite(1,-200,0)
% 
%         elseif mp==4
%             img = strrep( fullfile(data.imgfolder, rightstim), '.bmp', ['-', imagenumber, '.bmp'] );
%             cgloadbmp(1,img,stimDimensions(1),stimDimensions(2))
%             cgdrawsprite(1,200,0)
%             cgpencol(data.colorright)
%             cgrect(200,0,stimDimensions(1)+10,stimDimensions(2)+10)
%             cgdrawsprite(1,200,0)
%         end
% 
%         cgflip(data.background); 
%         wait(5)
%    
%     end
allcomparisons{compareIndex(i),4} = i; % Log the order in which each permutation was displayed   
allcomparisons{compareIndex(i),5} = mp; % No mouse press for slider ratings
allcomparisons{compareIndex(i),6} = rt; % No mouse press for slider ratings
    %allcomparisons{compareIndex(i),5} = rating; % log preference choice and strength
    mp
   wait(data.time.outcometime);
    
    cgfont('Arial',data.FontSizeTask);
    cgpencol(1,1,1);
    cgtext('+',0,0);
    cgfont('Arial',data.FontSizeText+10);
    cgtext('Outcomes added to total earnings',0,300);
    %cgtext('Which is best?',0,300);
    cgflip(data.background);
   wait(data.time.delay);

end

compareIndex = randperm(size(allcomparisons,1)); % Reshuffled and run again, switching order that it samples (column 2 first, 1 second)
compareIndex = [compareIndex, compareIndex]; % Concatenate and go from end of first block
allcomparisons = [allcomparisons;allcomparisons] % Concatenate and only go from end of first block

starti = i+1;

%%

for i=starti: size(allcomparisons,1)
    
    disp( allcomparisons{compareIndex(i),3} ) % Print the comparison to console
    stimDimensions = data.stimDimensions;
    
    cgfont('Arial',data.FontSizeText+10);
    cgtext('Outcomes added to total earnings',0,300);
    cgfont('Arial',data.FontSizeTask);
    cgpencol(1,1,1);
    cgtext('+',0,0);
   
    % Draw left stim
    leftstim = char( allcomparisons{compareIndex(i),1} );
    cgloadbmp(1,leftstim,stimDimensions(1),stimDimensions(2))
%     cgpencol(data.colorleft)
%     cgrect(-200,0,stimDimensions(1)+10,stimDimensions(2)+10)
    cgdrawsprite(1,-200,0)
    
    
    % Draw right stim
    rightstim = char( allcomparisons{compareIndex(i),2} );
%     cgpencol(data.colorright)
%     cgrect(200,0,stimDimensions(1)+10,stimDimensions(2)+10)
    cgloadbmp(2,rightstim,stimDimensions(1),stimDimensions(2))
    cgdrawsprite(2,200,0)
    cgpencol(1,1,1)
    
    
%     tstart = time;                %get time in ms from first time call
%     %wait(waittime);
%     tstart2 = time;
%     alltime = [tstart tstart2];
    cgflip(data.background);
    cgfont('Arial',data.FontSizeTask);
%     [~,~,~,mp]=cgmouse; mp=0; %no mouse press yet
%     
%     while (mp ~=1) && (mp ~=4), %1 from leftmouse, 4 for rightmouse, 2 for middle mouse
%         [~,~,~, mp] = cgmouse;
%     end;
%     rt = (time - tstart2)/1000; %rt in s
    %choicemade = sqrt(mp); %1 for left, 2 for right
    
    % Get line measurement too
    originy = -200;
    
    %sx = [-384 384]; sy = [originy, originy]; %in pixels, since it's 1024 wide, this is 75% of the screen
% if mp==1
%     sx = [-200-stimDimensions(1)/2 0]; sy = [originy, originy]; % for lefthandside
%     cgfont('Arial',data.FontSizeText);
%     cgtext('|',-200-stimDimensions(1)/2, 0);
%     cgtext('Strongly prefer',20,-300);
%     cgtext('left',40,-300);
% elseif mp==4
%     sx = [0 200+stimDimensions(1)/2]; sy = [originy, originy]; % for righthandside
%     cgtext('|',0,300);
%     cgtext('Strongly prefer',20,300);
%     cgtext('right',40,300);
% end

sx = [-200-stimDimensions(1)/2 200+stimDimensions(1)/2]; sy = [originy, originy]; % for lefthandside


    
    %randstart = rand;             %start in a random location on the line (0,1)
    randstart = 0.5;               %skip the randstart and just start from the middle
    rating = nan; rt = nan;        %in case early exit

    %cgfont('Arial',data.FontSizeText);
    cgpencol(1,1,1);
    tstart = time;                %get time in ms from first time call
    tstart2 = time;
    alltime = [tstart tstart2];
    
    %mp2 = 0; x = 220; y = 0; %initialize mouse to center of screen at moment line appears
    mp2 = 0; x = 0; y = 0; %initialize mouse to center of screen at moment line appears
    %cgmouse((sx(2)-sx(1))*randstart + sx(1), 0); %mouse jumps to random location on line
    cgmouse(0, 0); %mouse jumps to random location on line
    
     while (mp2 ~= 1) && time < tstart + timelimit;
         
    cgfont('Arial',data.FontSizeText+10);
    cgtext('Outcomes added to total earnings',0,300);
         cgfont('Arial',data.FontSizeText);
    cgtext('|',0, originy);
    %cgtext('No preference',0,originy-30);

%cgtext('|',-200-stimDimensions(1)/2, originy);
cgtext('Strongly confident left',-350, originy-30);
%cgtext('|',200+stimDimensions(1)/2,originy);
cgtext('Strongly confident right',350, originy-30);
%cgflip(data.background);
%         
        % Draw left stim
        
    leftstim = char( allcomparisons{compareIndex(i),1} );
    cgloadbmp(1,leftstim,stimDimensions(1),stimDimensions(2))
    cgpencol(data.colorleft)
    cgrect(-200,0,stimDimensions(1)+10,stimDimensions(2)+10)
    cgdrawsprite(1,-200,0)
    
    
    % Draw right stim
    rightstim = char( allcomparisons{compareIndex(i),2} );
    cgpencol(data.colorright)
    cgrect(200,0,stimDimensions(1)+10,stimDimensions(2)+10)
    cgloadbmp(2,rightstim,stimDimensions(1),stimDimensions(2))
    cgdrawsprite(2,200,0)
    cgpencol(1,1,1)
    cgpencol(0,0,0);
%    if mp==1
       % cgdrawsprite(1,-200,0)
       % cgrect(200,0,stimDimensions(1)+10,stimDimensions(2)+10)
 %   elseif mp==4
       % cgdrawsprite(2,200,0)
        %cgrect(-200,0,stimDimensions(1)+10,stimDimensions(2)+10)
  %  end  
%     
% cgfont('Arial',data.FontSizeText);
% cgpencol(1,1,1);
% if mp==1
%     sx = [-200-stimDimensions(1)/2 0]; sy = [originy, originy]; % for lefthandside
%     %cgtext('|',-200-stimDimensions(1)/2, originy);
%     cgtext('Strongly prefer',-350,originy-30);
%     %cgtext('LEFT',70,originy-90);
% elseif mp==4
%     sx = [0 200+stimDimensions(1)/2]; sy = [originy, originy]; % for righthandside
%     %cgtext('|',300,originy);
%     cgtext('Strongly prefer',350, originy-30);
%     %cgtext('RIGHT',300,originy-60);
% end
% %cgtext('|', 0, originy);
% cgtext('No preference',0, originy-30);
% %cgflip(data.background);

    
        cgpencol(1,1,1);
        cgfont('Arial',data.FontSizeTask);
        %cgtext('Which do you prefer?',0,300);
        cgpencol(0,0,0); cgdraw(sx(1),sy(1),sx(2),sy(2)); %draw rating line
        cgpencol(1,1,1); cgellipse(sx(1),sy(1),10,10,'f');
        cgellipse(sx(2),sy(2),10,10,'f');

        % update mouse position
        [x, ~, ~, mp2] = cgmouse; y = originy;
        if x < sx(1), x = sx(1); end
        if x > sx(2), x = sx(2); end;
        cgpencol(1,1,0); cgellipse(x,y,10,10,'f');
        cgfont('Arial',data.FontSizeTask);
        cgpencol(1,1,1);
        cgtext('+',0,0);

        cgflip(data.background);
    end
    rating = (x - sx(1)) / (sx(2) - sx(1)); %0 to 1
    rt = (time - tstart2)/1000; %in s since end of waittime

    allcomparisons{compareIndex(i),4} = i;
    allcomparisons{compareIndex(i),5} = rating; % log preference choice and strength
    allcomparisons{compareIndex(i),6} = rt;% log preference choice and strength
    
    mp
    rating
    wait(data.time.outcometime);
        
   
    
    
    %preferredstim=mp;
    
    
    
%     if rating==1
%         cgloadbmp(1,leftstim,stimDimensions(1),stimDimensions(2))
%         cgpencol(data.colorleft)
%         cgrect(-200,0,stimDimensions(1)+10,stimDimensions(2)+10)
%         cgdrawsprite(1,-200,0)
%     elseif mp==4
%         cgpencol(data.colorright)
%         cgrect(200,0,stimDimensions(1)+10,stimDimensions(2)+10)
%         cgloadbmp(2,rightstim,stimDimensions(1),stimDimensions(2))
%         cgdrawsprite(2,200,0)
%         cgpencol(1,1,1)
%     end
    
    cgfont('Arial',data.FontSizeTask);
    cgpencol(1,1,1);
    cgtext('+',0,0);
cgfont('Arial',data.FontSizeText+10);
    cgtext('Outcomes added to total earnings',0,300);
    cgflip(data.background); 
    
%     for iLever=1:length(data.fLevers)
% 
%             if iLever<10
%                imagenumber=['0', num2str(iLever)];
%             else
%                 imagenumber=num2str(iLever);
%             end
% 
%         cgpencol(1,1,1);
%         cgtext('+',0,0);
%         %cgloadbmp(1,data.pictureleft,294,154)
%         if mp==1
%             img = strrep( fullfile(data.imgfolder, leftstim), '.bmp', ['-', imagenumber, '.bmp'] );
%             cgloadbmp(1,img,stimDimensions(1),stimDimensions(2))
%             cgdrawsprite(1,-200,0)%cgdrawsprite(1,-200,0)
%             cgpencol(data.colorleft)
%             cgrect(-200,0,stimDimensions(1)+10,stimDimensions(2)+10)
%             cgdrawsprite(1,-200,0)
% 
%         elseif mp==4
%             img = strrep( fullfile(data.imgfolder, rightstim), '.bmp', ['-', imagenumber, '.bmp'] );
%             cgloadbmp(1,img,stimDimensions(1),stimDimensions(2))
%             cgdrawsprite(1,200,0)
%             cgpencol(data.colorright)
%             cgrect(200,0,stimDimensions(1)+10,stimDimensions(2)+10)
%             cgdrawsprite(1,200,0)
%         end
% 
%         cgflip(data.background); 
%         wait(5)
%    
%     end
    
 wait(data.time.outcometime);
    
    cgfont('Arial',data.FontSizeTask);
    cgpencol(1,1,1);
    cgpencol(1,1,1);
    cgtext('+',0,0);
    cgfont('Arial',data.FontSizeText+10);
    cgtext('Outcomes added to total earnings',0,300);
    cgflip(data.background);

    wait(data.time.delay);

end




