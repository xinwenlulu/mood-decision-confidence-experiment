files = dir('Mood&Confidence*.mat');
ques = dir('Questionnaires*.mat');
gender = [];
age = [];
for i=1:length(files),
    %subjectid = i;
    %subjectid = str2num(strrep(files(i).name(21:22),'_',''))
    load(files(i).name);
    load(ques(i).name);
    HPS = questionnaires.HPS.HPS.HPS_score;
    Altman = questionnaires.Altman.Altman.Altman_score;
    alldata(i).maintask = data.maintask;
    alldata(i).test = data.test;
    alldata(i).subjectid = str2num(strrep(files(i).name(21:22),'_',''));
    alldata(i).HPS = HPS;
    alldata(i).Altman = Altman;
    if data.gender(1) == 'f' || data.gender(1) == 'F'
        gender(end+1) = 0;
    else
        gender(end+1) = 1;
    end
    
    if length(data.dob) == 2
        age(end+1) = str2num(data.dob);
    else
        age(end+1) = 2019 - str2num(data.dob(end-3:end));
    end
end;
%%
%alldata = alldata(1);
for n=1:length(alldata), % Loop through all participants
    alldata(n).mLconfidence = nan(40,2);
    alldata(n).mHvLconfidence = nan(20,2);
    alldata(n).mHvHconfidence = nan(20,2);
    
    for m=1:length(alldata(n).maintask), %Loop over blocks
        t = alldata(n).maintask{m}(:,1:13); % get trial matrix for current subj and block 
        t(:,end+1) = alldata(n).maintask{m}(:,24); % append confidence 
        happy = t(:,10); % column containing happy rating
        alldata(n).wof_order(m) = t(1,13); %wof identity
        alldata(n).happy(:,m) = happy(~isnan(happy)); %happiness ratings for each block, nans removed
%         
        alldata(n).blockhappy(m) = mean(alldata(n).happy(:,m));
        alldata(n).hapchangeNeg = nan;
        alldata(n).hapchangePos = nan;
        if alldata(n).wof_order(m) == 0
            alldata(n).hapNeu1 = alldata(n).happy(1,m);
            alldata(n).hapNeuend = alldata(n).happy(end,m);
            alldata(n).hapNeu = alldata(n).blockhappy(m);
            alldata(n).hapNeg1 = nan;
            alldata(n).hapNeg = nan;
            alldata(n).hapPos1 = nan;
            alldata(n).hapPos = nan;
        elseif alldata(n).wof_order(m) == -1
            alldata(n).hapNeg1 = alldata(n).happy(1,m);
            alldata(n).hapNeg = alldata(n).blockhappy(m);
            alldata(n).hapchangeNeg = alldata(n).hapNeuend-alldata(n).hapNeg1;
        else
            alldata(n).hapPos1 = alldata(n).happy(1,m);
            alldata(n).hapPos = alldata(n).blockhappy(m);
            alldata(n).hapchangePos = alldata(n).hapPos1-alldata(n).hapNeuend;
        end
        %alldata(n).winhappy(m) = mean(t(t(:,8)>0,10)); %mean happy at next rating after win
        %alldata(n).losshappy(m) = mean(t(t(:,8)==0,10)); %after no win
        
        
        
        % Work out stim number for each condition
        for trial = 1:size(t,1)
            if t(trial,14) < 0.5
                t(trial,14) = 1 - t(trial,14);
            end
            if t(trial,4) == 1
                mLstim = t(trial,2);
            elseif t(trial,4) == 2
                mHvLstim = t(trial,2);
            else %condition == 3
                mHvHstim = t(trial,2);
            end
        end
        
        t(:,14) = 2.*(t(:,14)-0.5);
        alldata(n).confidence(:,m) = t(:,14);
        
        if m == 1
            alldata(n).confidenceNeu = nanmean(t(:,14));
            alldata(n).confidenceNeg = nan;
            alldata(n).confidencePos = nan;
        else
            alldata(n).confidenceNeg = nanmean(t(t(:,13)==-1,14));
            alldata(n).confidencePos = nanmean(t(t(:,13)==1,14));
        end
        
        
        mL_mHvL = [];
        mL_mHvH = [];
        mHvL_mHvH = [];
        allmL = [];
        allmHvL = [];
        allmHvH = [];
        for trial = 1:size(t,1)
            if t(trial,2) == mLstim || t(trial,3) == mLstim
                allmL(end+1,:) = t(trial,:);
                if t(trial,2) == mHvHstim || t(trial,3) == mHvHstim
                    mL_mHvH(end+1,:) = t(trial,:);         
                    allmHvH(end+1,:) = t(trial,:);
                else % mHvL
                    mL_mHvL(end+1,:) = t(trial,:);
                    allmHvL(end+1,:) = t(trial,:);
                end
            else %no mL means mHvL_mHvH
                mHvL_mHvH(end+1,:) = t(trial,:);
                allmHvH(end+1,:) = t(trial,:);
                allmHvL(end+1,:) = t(trial,:);   
            end
        end
        
%         for trial = 1:size(mL_mHvL,1)
%             if mL_mHvL(trial,5) == 1
%                 mL_mHvL(trial,14) = 1 - mL_mHvL(trial,14);
%             end
%              alldata(n).con_mL_mHvL(trial,m) = mL_mHvL(trial,14);
%         end
%         
%         for trial = 1:size(mL_mHvH,1)
%             if mL_mHvH(trial,5) == 1
%                 mL_mHvH(trial,14) = 1 - mL_mHvH(trial,14);
%             end
%              alldata(n).con_mL_mHvH(trial,m) = mL_mHvH(trial,14);
%         end
        % Store in all data which is hi and low stim for this subj and this
        % block (and average earnings for each to check)
        % And whether chose correct stim on each trial or not
        alldata(n).mLstim(m) = mLstim; 
        alldata(n).mHvLstim(m) = mHvLstim;
        alldata(n).mHvHstim(m) = mHvHstim;
        alldata(n).mLreward(m) = mean(t(t(:,7)==mLstim,8));
        
        alldata(n).mLconfidence(allmL(:,7)== mLstim,m) = allmL(allmL(:,7)==mLstim,14);
        alldata(n).mHvLconfidence(mL_mHvL(:,7)== mHvLstim,m) = mL_mHvL(mL_mHvL(:,7)==mHvLstim,14);
        alldata(n).mHvHconfidence(mL_mHvH(:,7)== mHvHstim,m) = mL_mHvH(mL_mHvH(:,7)==mHvHstim,14);
        alldata(n).mHvLreward(m) = mean(t(t(:,7)==mHvLstim,8));
        alldata(n).mHvHreward(m) = mean(t(t(:,7)==mHvHstim,8));
         % Mean earnings (per trial) 

        alldata(n).maintask_mLchoice(:,m) = allmL(:,7)==mLstim; % whether chose the bad option
        alldata(n).maintask_mHvLchoice(:,m) = mL_mHvL(:,7) == mHvLstim;
        alldata(n).maintask_mHvHchoice(:,m) = mL_mHvH(:,7) == mHvHstim;
        
        %%% Count up %switch after a loss
        %cWinstay=0; cLoseswitch=0;cWinswitch=0; cLosestay=0;
        %trials=t; % all trials
        %trials=t2; % free choice only
        
        %for itrial=2:size(trials, 1)
        %    if trials(itrial-1,8)>0 % Previous trial was a win
        %        if ( trials(itrial,7) == trials(itrial-1,7) ) && ( and(trials(itrial,2)>0,trials(itrial,3)>0) ) % Win stay, unforced
        %            cWinstay = cWinstay+1;
        %        elseif ( trials(itrial,7) ~= trials(itrial-1,7) ) && ( and(trials(itrial,2)>0,trials(itrial,3)>0) ) % Win switch, unforced
        %             cWinswitch=cWinswitch+1;
        %        end
        %    else % Previous trial was a no-win
        %        if ( trials(itrial,7) ~= trials(itrial-1,7) ) && ( and(trials(itrial,2)>0,trials(itrial,3)>0) ) % Lose switch, unforced
        %            cLoseswitch=cLoseswitch+1;
        %        elseif ( trials(itrial,7) == trials(itrial-1,7) ) && ( and(trials(itrial,2)>0,trials(itrial,3)>0) ) % Lose stay, unforced
        %             cLosestay = cLosestay+1;
        %        end
        %    end
        %end
    %alldata(n).winstay(m) = cWinstay/(cWinstay+cWinswitch); %size(trials,1);
    %alldata(n).loseswitch(m) = cLoseswitch/(cLoseswitch+cLosestay); %size(trials,1);
    %alldata(n).winswitch(m) = cWinswitch/(cWinstay+cWinswitch); %size(trials,1);
    %alldata(n).losestay(m) = cLosestay/(cLoseswitch+cLosestay); %size(trials,1);
    
%     trials(2:end,9:10)       = trials(1:end-1,7:8); trials(1,:) = []; %make columns 9 and 10 previous trial choice/outcome, toss first trial
%     trials                   = trials(trials(:,2)&trials(:,3),1:10); %only need to know what happened on choice trials relative to previous trial
%     alldata(n).winstay2(m)   = mean(trials(trials(:,10)>0,9)==trials(trials(:,10)>0,7));
%     alldata(n).loseswitch2(m) = mean(trials(trials(:,10)==0,9)~=trials(trials(:,10)==0,7));

%flips=[]; cFlipcorrect=0; cFlipincorrect=0;
%cNoflipcorrect=0;cNoflipincorrect=0;
%for itrial=2:size(trials, 1)
    
     %if trials(itrial-1,2)~=trials(itrial,2) && trials(itrial-1,3)~=trials(itrial,3) && sum(sum(trials(itrial-1:itrial,3:4)==0))<3  % Both LHS and RHS have changed therefore flip (not just forced choice)
%     if ( trials(itrial,2) ~= trials(itrial-1, 2) && ~any([trials(itrial, 2), trials(itrial-1, 2)]==0) )  || ( trials(itrial,3) ~= trials(itrial-1, 3) && ~any([trials(itrial, 3), trials(itrial-1, 3)]==0) )
         %if sum(sum(trials(itrial-1:itrial,2:3)==0))<2     
%         flips = [flips, itrial];
         
%         if trials(itrial,2) > 0 && trials(itrial,3) > 0 % free choice so count it
%             if trials(itrial,7)==histim % They chose correct
%                 cFlipcorrect = cFlipcorrect +1;
%             else
%                 cFlipincorrect =  cFlipincorrect+1;
%             end

%         end
    
%         if trials(itrial-1,2) > 0 && trials(itrial-1,3) > 0 
%             if trials(itrial-1,7)==histim % They chose correct on trial before flip
%                 cNoflipcorrect = cNoflipcorrect +1;
%             else
%                 cNoflipincorrect =  cNoflipincorrect+1;
%             end
%         end
             
%     end
     
%end
    %alldata(n).flips(:,m) = flips;
    %disp(['s', num2str(n), ' - flips b', num2str(m), ': ', num2str(flips)])
    %alldata(n).cFlipcorrect(m) = cFlipcorrect/(cFlipcorrect + cFlipincorrect ); 
    %alldata(n).cFlipincorrect(m) = cFlipincorrect /(cFlipcorrect + cFlipincorrect ); 
    %alldata(n).cNoflipcorrect(m) = cNoflipcorrect/(cNoflipcorrect + cNoflipincorrect ); 
    %alldata(n).cNoflipincorrect(m) = cNoflipincorrect /(cNoflipcorrect + cNoflipincorrect ); 
    
    end; %end of block loop
    
    % store mean performance and happiness for each block in alldata
    alldata(n).choosing_mL = mean(alldata(n).maintask_mLchoice);
    alldata(n).choosing_mHvL = mean(alldata(n).maintask_mHvLchoice);
    alldata(n).choosing_mHvH = mean(alldata(n).maintask_mHvHchoice); 
    %alldata(n).learning_correct = mean(alldata(n).maintask_choice); % calculate mean performance (accuracy %) for each block

    
%%% This wof part assumes all subjs have a positive and negative (i.e 3 blocks)
%%% However doesn't apply for smartphone where participant either had a positive or a negative mood induction

%      if alldata(n).wof_order(2)<0 %if wof loss first
%          alldata(n).happy_reorder = alldata(n).happy(:,[1 2]);
%          alldata(n).meanhappy_reorder = alldata(n).meanhappy([1 3 2]);
%     else
%         alldata(n).happy_reorder = alldata(n).happy;
%         alldata(n).meanhappy_reorder = alldata(n).meanhappy;

end;

%disp('--Mean correct/incorrect upon flip:')
%[nanmean(nanmean(vertcat(alldata.cFlipcorrect))), nanmean(nanmean(vertcat(alldata.cFlipincorrect))) ]
%disp('--Mean correct/incorrect on flip-1:')
%[nanmean(nanmean(vertcat(alldata.cNoflipcorrect))), nanmean(nanmean(vertcat(alldata.cNoflipincorrect))) ]
%%
% Plot choice curve averaged across subjects and blocks
disp('---fraction of choices correct in each block for all subjects:');
%disp(vertcat(alldata.choosing_mL));
%disp(vertcat(alldata.choosing_mHvL));
%disp(vertcat(alldata.choosing_mHvH));
%disp(vertcat(alldata.learning_correct));
figure('color',[1 1 1]); %subplot(1,1,1);     

%t=horzcat(alldata.maintask_choice);% maintask_choice is trialwise correct or not. This concatenates all blocks one after the other
for i = 1:24
    subplot(4,6,i)
    %subplot(1,1,1)
%t=horzcat(alldata(i).maintask_mLchoice);
%t=horzcat(alldata(i).maintask_mHvLchoice);
t=horzcat(alldata(i).maintask_mHvHchoice);
%     confidence = nan(size(alldata(1).confidence,1),2*length(alldata));
%     for k = 1:length(alldata)
%         confidence(:,2*k-1) = alldata(k).confidence(:,1);
%         confidence(:,2*k) = alldata(k).confidence(:,2);
%     end
%     t = horzcat(confidence);
%     confidence = nan(size(alldata(1).con_mL_mHvL,1),2);%*length(alldata));
%     for k = 1:length(alldata)
%         confidence(:,2*k-1) = alldata(i).con_mL_mHvL(:,1);
%         confidence(:,2*k) = alldata(i).con_mL_mHvL(:,2);
%     end
%     t = horzcat(confidence);

indB1=1; indB2 = 2;
  %indB1 = 1:2:length(alldata)*2; % B1 = Block 1
  %indB2 = 2:2:length(alldata)*2; % B2 = Block 2
%indB3 = 3:3:length(alldata)*3; % B3 = Block 3

% block1plot = mean(confidence(:,indB1),2);
% block2plot = mean(confidence(:,indB2),2);

% % block = zeros(size(t,1),2*length(alldata));
 block = zeros(size(t,1),2);
for b = 1:2%*length(alldata)
    choice = 0;
    cumchoice = [];
    for n = 1:size(t,1)
        choice = choice + t(n,b);
        cumchoice(end + 1) = choice/n;  
    end
   block(:,b) = cumchoice;
end
block1plot = mean(block(:,indB1),2);
block2plot = mean(block(:,indB2),2);
%  
%   block1plot = mean(t(:,indB1),2);
%   block2plot = mean(t(:,indB2),2);
%block3plot = mean(t(:,indB3),2);
Avgplot = mean(block,2);
%Avgplot = mean(t,2); 
%Avgplot = mean(confidence,2);
    hold on
    plot(Avgplot, 'k', 'LineWidth',2); 
    plot(block1plot,'color','b'); 
    plot(block2plot,'color','g'); 
%plot(block3plot,'color','r'); 


%%%Plot with error bars instead

    %legend('Average', 'Block 1', 'Block 2');
%  semB1 = std(t(:,indB1),1)./sqrt(size(t(:,indB1),2));]
%  semB2 = std(t(:,indB2),1)./sqrt(size(t(:,indB2),2));
%  errorbar( mean(t(:,indB1),2), semB1 )
%  errorbar( mean(t(:,indB2),2), semB2 )
%axis([0 18 0 1]); % This assumes 18 free choice
    axis([0 size(t,1) 0 1]); % Rescale to correct number of free choice trials
    hold on; plot(xlim,[0.5 0.5],'k--'); 
    title(sprintf('id = %d', alldata(i).subjectid))
    %title(sprintf('n = %d',length(alldata))); 
    %title(sprintf('n = 1'));
    xlabel('choice trial number'); 
    ylabel('mHvH');
    %ylabel('confidence');

% %%% Throw out sessions where subjects completed multiple runs - 
% %%% N.B. 'firstrun' NEEDS SPECIFYING MANUALLY!
% t=vertcat(alldata.learning_correct);
% subjsfirstrunonly = [1;0;0;1;1;0;0;0;1;0;0;0;1;0;1;1;1;1;1]; % This corresponds to the .mats in the 'files' variable
% disp('Block mean performance throwing out any repeated runs')
% disp(mean(t(find(subjsfirstrunonly==1),:)))
% disp(mean(mean(t(find(subjsfirstrunonly==1),:))))
end
%% Plot choice curves for first and second block
%t=horzcat(alldata.maintask_choice); % PCENT correct, concatenated for each block and each subj
%t=horzcat(alldata.choosing_mL);
%t=horzcat(alldata.choosing_mHvL);
%t=horzcat(alldata.choosing_mHvH);

% %%% I think this has a bug, see comments - LM
% halfway = round(size(t,1)/2); % halfway through number of subjects...
% t2=[mean(t(1:halfway,:),1); mean(t(halfway+1:end,:),1)]; % This actually averages blocks separately for first and second half of subjects
% %%
t2=[mean(t(:,indB1),1); mean(t(:,indB2),1)]; % Average of B1 and Average of B2 and Average of B3 together

subplot(1,1,1); plot(t2); set(gca,'xtick',[1 2],'xticklabel',{'1st block','2nd block'}); 
 axis([0 3 0 1]); 
ylabel('percent mHvH'); hold on; plot(xlim,[0.5 0.5],'k--');

%% Scatterplot choice correct across blocks
%t=vertcat(alldata.learning_correct); % Average pcent correct per block for each subj
t=vertcat(alldata.choosing_mL);
%t=vertcat(alldata.choosing_mHvL);
%t=vertcat(alldata.choosing_mHvH);

%subplot(2,3,3); hold on; plot(t(:,2),t(:,3),'.'); % Assumes 3 blocks
subplot(1,1,1); hold on; plot(t(:,1),t(:,2),'.'); % Assumes 2 blocks
axis([0 1 0 1]); plot([0.5 0.5],ylim,'k'); plot(xlim,[0.5 0.5],'k');
xlabel('block 1'); ylabel('block 2'); axis square;

%% Plot mean happiness for each of the blocks by WoF

% %%% Much of this assumed that all subjs have both pos and neg and relies on blocks
% %%% having been reordered. Not true for smartphone, commented out and
% replaced.
% disp('happiness after neutral, positive, and negative blocks:');
% disp(vertcat(alldata.meanhappy_reorder));


% run1=vertcat(alldata.taskver)==1
% run2=vertcat(alldata.taskver)==2
% t(run1)
%t=vertcat(alldata.meanhappy); % All mean happiness ratings, per block and subject

for i = 1:length(alldata)
    if length(alldata(i).happy) == 4
        alldata(i).happy(5,:) = alldata(i).happy(4,:);
    end
end
    
    
t = vertcat(alldata.happy);
% happyRating = nan(size(t,1)/2, size(t,2)*2);
% for i = 1:size(t,1)/2
%     happyRating(i,:) = [t(2*i-1,1) t(2*i,1) t(2*i-1,2) t(2*i,2)]
% end
happyRating = nan(length(alldata),10);
for i = 1:length(alldata)
        %happyRating(i,:) = [t(4*i-3,1) t(4*i-2,1) t(4*i-1,1) t(4*i,1) t(4*i,1) t(4*i-3,2) t(4*i-2,2) t(4*i-1,2) t(4*i,2) t(4*i,1)];
     happyRating(i,:) = [t(5*i-4,1) t(5*i-3,1) t(5*i-2,1) t(5*i-1,1) t(5*i,1) t(5*i-4,2) t(5*i-3,2) t(5*i-2,2) t(5*i-1,2) t(5*i,1)];
end



wofoutcomes=vertcat(alldata.wof_order); % Mood state each block was completed under, per subject
Indposblock=find(wofoutcomes(:,2)== 1);
Indnegblock=find(wofoutcomes(:,2)== -1);

poshappyavg = nan(length(Indposblock),10);
neghappyavg = nan(length(Indnegblock),10);
%allhappyavg=nan(length(alldata),2); % Make a matrix of all block averages for happiness ratings
poshappyavg(:,:) = happyRating(Indposblock,:);
neghappyavg(:,:) = happyRating(Indnegblock,:);
%allhappyavg(:,1) = t(:,1); % Column 1 is neutral/B1
%posblock1=t(:,2); % temporarily fill both mood block averages
%negblock1=t(:,2);
%posblock2=t(:,3);
%negblock2=t(:,3);
%posblock1(Indnegblock)=nan; % Set the other mood block averages to nan
%negblock1(Indposblock)=nan;
%posblock2(Indposblock)=nan;
%negblock2(Indnegblock)=nan;
%posblock1(isnan(posblock1)) = posblock2(isnan(posblock1))
%negblock1(isnan(negblock1)) = negblock2(isnan(negblock1))

%allhappyavg(:,2) = posblock1; % Column 2 is positive
%allhappyavg(:,3) = negblock1; % Column 3 is negative

% subplot(2,3,4); plot(allhappyavg);
% legend('Baseline', 'Positive', 'Negative')

disp('---happiness for neutral and positive blocks:')
%disp(poshappyavg)
disp('---means for neutral and positive blocks:')
%disp(nanmean(poshappyavg))

disp('happiness for neutral and negative blocks:')
disp(neghappyavg)
disp('means for neutral and negative blocks')
disp(nanmean(neghappyavg))

% t = horzcat(alldata.happy_reorder);
% happy_reorder = {mean(t(:,1:3:end)),mean(t(:,2:3:end)),mean(t(:,3:3:end))};

%t=horzcat(alldata.happy); % Concatenate all

happyB1=neghappyavg(:,1:5); % Work out which is B1 and B2
happyB2=neghappyavg(:,6:10); %Block 2
%happyB3=t(:,indB3); %Block 3
everyhappy = [happyB1 happyB2];

happypos1=poshappyavg;%(Indposblock); % temporarily fill both mood block averages
happyneg1=neghappyavg;%(Indnegblock);
%happypos2=happyB3;
%happyneg2=happyB3;

%happypos1(:,Indnegblock)=nan; % Set the other mood block averages to nan
%happyneg1(:,Indposblock)=nan;
%happypos2(:,Indposblock)=nan;
%happyneg2(:,Indnegblock)=nan;

%happypos1(isnan(happypos1)) = happypos2(isnan(happypos1))
%happyneg1(isnan(happyneg1)) = happyneg2(isnan(happyneg1))
subplot(1,1,1); hold on; mike={'b','y','c','m'}; %color order
%for n=1:length(happy_reorder), plot(happy_reorder{n},mike{n}); end; ylim([0 1]);
for i = 1:length(Indposblock)
    plot(poshappyavg(i,:))
end
% plot(poshappyavg(2,:))
% plot(poshappyavg(3,:))
% plot(poshappyavg(4,:))
 plot(mean(poshappyavg),'k')
% for i = 1:length(Indnegblock)
% plot(neghappyavg(i,:))%,'b')
% end
% plot(neghappyavg(6,:))
% plot(neghappyavg(7,:))
% plot(neghappyavg(8,:))
% plot(neghappyavg(9,:))
%plot(mean(neghappyavg),'k')
%plot(neghappyavg,'r')
%plot(mean(happyB1,2),'k')
%plot(nanmean(happypos1,2),'g')
%plot(nanmean(happyneg1,2),'r')

ylim([0 1]);
%ylim([0 .6]);
set(gca,'xtick',[1 2 3 4 5 6 7 8 9 10],'xticklabel',{'1','2','3', '4','5','6','7','8','9','10'}); 
axis([1 10 0 1]); 
xlabel('happiness rating'); ylabel('happiness');
%legend(['participant',num2str(1)])%,['participant',num2str(2)],['participant',num2str(3)],['participant',num2str(4)]);

x=(happypos1-repmat(mean(everyhappy),size(happypos1,1),1));
y=repmat(std(everyhappy),size(happypos1,1),1);
zscorehappypos = x-y;
subplot(1,1,1); hold on
plot(zscorehappypos(1,:))
plot(zscorehappypos(2,:))
plot(mean(zscorehappypos))

x=(happyneg1-repmat(mean(everyhappy),size(happyneg1,1),1));
y=repmat(std(everyhappy),size(happyneg1,1),1);
zscorehappyneg = x-y;
subplot(1,1,1); hold on
plot(zscorehappyneg(1,:))
plot(zscorehappyneg(2,:))
plot(zscorehappyneg(3,:))
plot(zscorehappyneg(4,:))
plot(zscorehappyneg(5,:))
plot(zscorehappyneg(6,:))
plot(zscorehappyneg(7,:))
plot(zscorehappyneg(8,:))
plot(zscorehappyneg(9,:))
plot(mean(zscorehappyneg),'k')

x=(happyB1-repmat(mean(everyhappy),size(happyB1,1),1));
y=repmat(std(everyhappy),size(happyB1,1),1);
zscorehappyB1 = x-y;

mooddiff_PoverN = mean(happypos1)-mean(happyneg1);
mooddiff_PoverNeu = mean(happypos1)-mean(happyB1);
mooddiff_NoverNeu = mean(happyB1)-mean(happyneg1);
[mean(mooddiff_PoverN), mean(mooddiff_PoverNeu), mean(mooddiff_NoverNeu)]

mooddiff_PoverN_zscore = mean(zscore(happypos1))-mean(zscore(happyneg1));
mooddiff_PoverNeu_zscore = mean(zscore(happypos1))-mean(zscore(happyB1));
mooddiff_NoverNeu_zscore = mean(zscore(happyB1))-mean(zscore(happyneg1));
[mean(mooddiff_PoverN_zscore), mean(mooddiff_PoverNeu_zscore), mean(mooddiff_NoverNeu_zscore)]

mooddiff_PoverN_zscore = mean( zscore(happypos1) - zscore(happyneg1) );
mooddiff_PoverNeu_zscore = mean( zscore(happypos1)- zscore(happyB1) );
mooddiff_NoverNeu_zscore = mean( zscore(happyB1)- zscore(happyneg1) );
[mean(mooddiff_PoverN_zscore), mean(mooddiff_PoverNeu_zscore), mean(mooddiff_NoverNeu_zscore)]

%% Scatterplot happiness block after positive vs negative relative to neutral

%t=vertcat(alldata.meanhappy_reorder);
subplot(2,3,5); hold on; 
%plot(t(:,2)-t(:,1),t(:,3)-t(:,1),'.');

%happydiff_posneu = nanmean(happyPos,1)-mean(happyB1,1);
%happydiff_negneu = nanmean(happyNeg,1)-mean(happyB1,1);
%plot(happydiff_posneu,happydiff_negneu,'.'); 
%xlabel('happiness Positive-Neutral block'); ylabel('happiness Negative - Neutral block'); axis square;
%axis([-1 1 -1 1]*0.5); plot([0 0],ylim,'k'); plot(xlim,[0 0],'k');


plot( nanmean(happyB1,1), nanmean(happypos1,1), '.', 'Color', 'g', 'MarkerSize',15)
plot( nanmean(happyB1,1), nanmean(happyneg1,1), '.', 'Color', 'r', 'MarkerSize',15 )
axis([0 1 0 1]); plot(xlim, ylim, 'k'); axis square;
xlabel('happiness Neutral block'); ylabel('happiness Moood block'); legend('Positive', 'Negative')



%% Backfill happiness, check if higher after positive than negative outcome

winhappy = vertcat(alldata.winhappy); % backfilled happy for B1 and B2
losshappy = vertcat(alldata.losshappy); % backfilled happy for B1 and B2
subplot(2,3,6); hold on; 
plot(winhappy(:,2)-losshappy(:,2),winhappy(:,3)-losshappy(:,3),'.'); %
%Assumes 3 blocks
%plot(winhappy(:,1)-losshappy(:,1),winhappy(:,2)-losshappy(:,2),'.');
axis([-0.05 0.5 -0.05 0.5]); plot([0 0],ylim,'k'); plot(xlim,[0 0],'k'); axis square;
xlabel('win-loss happiness (block 1)'); ylabel('win-loss happiness (block 2)');

verify the high and low stimuli were paired with right number of rewards
disp('---- rate of reward for hi prob stim')
disp(vertcat(alldata.hireward)); %looks good - only off because of variance in choices
disp('---- rate of reward for lo prob stim')
disp(vertcat(alldata.loreward)); 
%%
for n=1:length(alldata) % Loop through all participants        
    t = alldata(n).test{1}; % get trial matrix for current subj and block 
        
    mL = nan(0,14);
    mHvL = nan(0,14);
    mHvH = nan(0,14);

    mL_mHvL1 = [];
    mL_mHvL2 = [];
    mL_mHvH1 = [];
    mL_mHvH2 = [];
    mHvL_mHvH1 = [];
    mHvL_mHvH2 = [];
    
        for trial = 1:size(t,1)
            if t(trial,4) == t(trial,5) %different blocks
                if t(trial,4) == 1
                    mL(end+1,:) = t(trial,:);
                elseif t(trial,4) == 2
                    mHvL(end+1,:) = t(trial,:);
                else % condition = 3
                    mHvH(end+1,:) = t(trial,:);
                end
            elseif t(trial,6) == t(trial,7) == 0 %block1 same block
                if t(trial,2) == alldata(n).mLstim(1) || t(trial,3) == alldata(n).mLstim(1)
                    if t(trial,2) == alldata(n).mHvLstim(1) || t(trial,3) == alldata(n).mHvLstim(1)
                        mL_mHvL1(end+1,:) = t(trial,:);
                    else % mHvHstim
                        mL_mHvH1(end+1,:) = t(trial,:);
                    end
                else %mHvL_mHvH
                    mHvL_mHvH1(end+1,:) = t(trial,:);
                end
            else %block2 same block
                if t(trial,2) == alldata(n).mLstim(2) || t(trial,3) == alldata(n).mLstim(2)
                    if t(trial,2) == alldata(n).mHvLstim(2) || t(trial,3) == alldata(n).mHvLstim(2)
                        mL_mHvL2(end+1,:) = t(trial,:);
                    else % mHvHstim
                        mL_mHvH2(end+1,:) = t(trial,:);
                    end
                else %mHvL_mHvH
                    mHvL_mHvH2(end+1,:) = t(trial,:);
                end
                
            end
        end
        
        
        mLchoice = mL(:,12);
        mHvLchoice = mHvL(:,12);
        mHvHchoice = mHvH(:,12);
        
        preferpos = [];
        preferneg = [];
        
        for b = 1:2
            if mLchoice(b) == alldata(n).mLstim(1) %chose block1
                if alldata(n).wof_order(2) == -1
                    preferpos(end+1) = 1;
                    preferneg(end+1) = 0;
                else %block2 wof == 1
                    preferpos(end+1) = 0;
                    preferneg(end+1) = 1;
                end
            else % chose block2
                if alldata(n).wof_order(2) == 1
                    preferpos(end+1) = 1;
                    preferneg(end+1) = 0;
                else %wof =-1
                    preferpos(end+1) = 0;
                    preferneg(end+1) = 1;
                end
            end
            
             if mHvLchoice(b) == alldata(n).mLstim(1) %chose block1
                if alldata(n).wof_order(2) == -1
                    preferpos(end+1) = 1;
                    preferneg(end+1) = 0;
                else %block2 wof == 1
                    preferpos(end+1) = 0;
                    preferneg(end+1) = 1;
                end
            else % chose block2
                if alldata(n).wof_order(2) == 1
                    preferpos(end+1) = 1;
                    preferneg(end+1) = 0;
                else %wof =-1
                    preferpos(end+1) = 0;
                    preferneg(end+1) = 1;
                end
             end
            
             if mHvHchoice(b) == alldata(n).mLstim(1) %chose block1
                if alldata(n).wof_order(2) == -1
                    preferpos(end+1) = 1;
                    preferneg(end+1) = 0;
                else %block2 wof == 1
                    preferpos(end+1) = 0;
                    preferneg(end+1) = 1;
                end
            else % chose block2
                if alldata(n).wof_order(2) == 1
                    preferpos(end+1) = 1;
                    preferneg(end+1) = 0;
                else %wof =-1
                    preferpos(end+1) = 0;
                    preferneg(end+1) = 1;
                end
            end
        end

        t = [mL_mHvL1;mL_mHvL2;mL_mHvH1;mL_mHvH2];              
        correct = [];
        for x = 1:size(t,1)
            if t(x,12) ~= mL
                correct(end+1) = 1;
            else
                correct(end+1) = 0;
            end
        end
        
         mHvL = [];
         mHvH = [];
         
    for y = 1:2
        if mHvL_mHvH1(y,12) == alldata(n).mHvLstim(1)
            mHvL(end+1) = 1;
            mHvH(end+1) = 0;
        else
            mHvL(end+1) = 0;
            mHvH(end+1) = 1;
        end
        
        if mHvL_mHvH2(y,12) == alldata(n).mHvLstim(2)
            mHvL(end+1) = 1;
            mHvH(end+1) = 0;
        else
            mHvL(end+1) = 0;
            mHvH(end+1) = 1;
        end
    end
    
        
        alldata(n).preferpos = sum(preferpos)/length(preferpos);
        alldata(n).preferneg = sum(preferneg)/length(preferneg);
        alldata(n).correct = sum(correct)/length(correct);
        alldata(n).mHvL = sum(mHvL)/length(mHvL);
        alldata(n).mHvH = sum(mHvH)/length(mHvH);
        
    
end;


