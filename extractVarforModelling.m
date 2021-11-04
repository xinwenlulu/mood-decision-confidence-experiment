%extract Choices,LeftStim,RightStim,Outcome
files = dir('Mood&Confidence*.mat');
for i=1:length(files)
    load(files(i).name);
    alldata(i).maintask = data.maintask;
    alldata(i).subjectid = data.subjectid;
    alldata(i).wof = data.wofdraws;
end;
%%

for n=1:length(alldata) % Loop through all participants
    mHvLcon = nan(0,1);
    mHvHcon = nan(0,1);

    for m=1:length(alldata(n).maintask), %Loop over blocks
        
        t = alldata(n).maintask{m}(:,1:13); % get trial matrix for current subj and block
        % raw confidence at column 24 (<.5 chose left >0.5 chose right )
        confidence = alldata(n).maintask{m}(:,24); 
        
        % make all the confidence data the confidence for the chosen choice
        for trial = 1:length(confidence)
            if confidence(trial) < 0.5
                confidence(trial) = 1 - confidence(trial); %now all confidence between 0.5 and 1
            end
        end
        confidence = 2.*(confidence - 0.5); % rescale the confidence data so it's between 0 adn 1
        
        alldata(n).modellingVar{m} = nan(size(t,1),4);
        
        for itrial = 1:size(t,1)
            leftStim = t(itrial,4);
            rightStim = t(itrial,5);
            if t(itrial,6) == 1 % chose left
                choice = t(itrial,4); % condition of the leftStim
            elseif t(itrial,6) == 2 % chose rigt
                choice = t(itrial,5); % condition of the rightStim
            end
            outcome = t(itrial,8);
            
            % writing these variables in the ModellngVar mtx
                alldata(n).modellingVar{m}(itrial,1) = choice;
                alldata(n).modellingVar{m}(itrial,2) = leftStim;
                alldata(n).modellingVar{m}(itrial,3) = rightStim;
                alldata(n).modellingVar{m}(itrial,4) = outcome;
                
                if itrial >=34 && itrial <=43 % get confidence for the last 10 bad vs good trials
                    if choice == 2 
                        mHvLcon(end+1,1) = confidence(itrial); % confidence for mHvL in the last 10 trials
                    elseif choice == 3
                        mHvHcon(end+1,1) = confidence(itrial); % cofidence for mHvH in the last 10 trials
                    end
                end
                
        end
        
      %calculate mean confidence for mHvL and mHvH for each participant
                alldata(n).mHvLcon(m) = nanmean(mHvLcon);
                alldata(n).mHvHcon(m) = nanmean(mHvHcon);
            
    end; %end of block loop

end;


