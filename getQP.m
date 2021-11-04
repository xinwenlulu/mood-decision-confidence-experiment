load resultsQalltrials.mat
load resultsSSATalltrials.mat
%these results filese were obtained using ScriptForXin.m

for i = 1:length(alldata)
    modelconfidence(i).wof = alldata(i).wof;
    for m = 1:2 % 2 blocks (1 before wof and 2 after wof)

Choices = alldata(i).modellingVar{m}(:,1);
LeftStim = alldata(i).modellingVar{m}(:,2);
RightStim = alldata(i).modellingVar{m}(:,3);
Outcome = alldata(i).modellingVar{m}(:,4);
%using the parameters fitted by the QLearning models
Param = [resultsQ(i).beta(m) resultsQ(i).LR(m) resultsQ(m).drift(m)];

% get Q values for all 44 trials
Q = QLearn(Choices,LeftStim,RightStim,Outcome,Param);

% using the parameters fitted by the SSAT model
Param = [resultsSSAT(i).beta(m) resultsSSAT(i).LR(m) resultsSSAT(i).drift(m) resultsSSAT(i).LRVar(m) resultsSSAT(i).threshold(m)];
% get P values for all 44 trials
PChoice = PSSAT(Choices,LeftStim,RightStim,Outcome,Param);

% get the Q and P values for the chosen option in the last 10 bad VS good
% trials
Qchosen = nan(10,3);
Pchosen = nan(10,3);
for trial = 1:10
    Qchosen(trial,1) = Q(trial+33,Q(trial+33,4));
    Qchosen(trial,2) = Q(trial+33,4); % chosen option
    Qchosen(trial,3) = Q(trial+33,5); % unchosen option
    Pchosen(trial,1) = PChoice(trial+33,PChoice(trial+33,4));
    Pchosen(trial,2) = PChoice(trial+33,4);
    Pchosen(trial,3) = PChoice(trial+33,5);
end

% separate the Q and P values for mL_mHvL and mL_mHvH comparisons
PmHvL = nan(0,1);
PmHvH = nan(0,1);
QmHvL = nan(0,1);
QmHvH = nan(0,1);
for trial = 1:size(Qchosen,1)
    if Qchosen(trial,2) == 2 % if chose mHvL
        QmHvL(end+1,1) = Qchosen(trial,1);
    elseif Qchosen(trial,2) == 3 % if chose mHvH
        QmHvH(end+1,1) = Qchosen(trial,1);
    end
    
    if Pchosen(trial,2) == 2
        PmHvL(end+1,1) = Pchosen(trial,1);
    elseif Qchosen(trial,2) == 3
        PmHvH(end+1,1) = Pchosen(trial,1);
    end
end

% calculate a mean Q value for mHvL and mHvH from the last 10 trials for
% each participant
QmHvL = mean(QmHvL);
QmHvH = mean(QmHvH);

modelconfidence(i).QmHvL(m) = QmHvL;
modelconfidence(i).QmHvH(m) = QmHvH;

% calculate a mean P value for mHvL and mHvH from the last 10 trials for each participant      
PmHvL = mean(PmHvL);
PmHvH = mean(PmHvH);

 modelconfidence(i).PmHvL(m) = PmHvL;
 modelconfidence(i).PmHvH(m) = PmHvH;
    end
end

