% induced_PoverN = mean(happypos1)>mean(happyneg1);
% induced_PoverNeu = mean(happypos1)>mean(happyB1);
% induced_NoverNeu = mean(happyneg1)<mean(happyB1);
% [sum(induced_PoverN), sum(induced_PoverNeu), sum(induced_NoverNeu)]
% sum(and(and(induced_PoverN, induced_PoverNeu),induced_NoverNeu))

mooddiff_PoverN = mean(happypos1)-mean(happyneg1);
mooddiff_PoverNeu = mean(happypos1)-mean(happyB1);
mooddiff_NoverNeu = mean(happyB1)-mean(happyneg1);
[mean(mooddiff_PoverN), mean(mooddiff_PoverNeu), mean(mooddiff_NoverNeu)]


for i=1:length(mooddiff_PoverN)
%     alldata(i).induced_PoverN = induced_PoverN(i);
%     alldata(i).induced_PoverNeu = induced_PoverNeu(i);
%     alldata(i).induced_NoverNeu = induced_NoverNeu(i);
    alldata(i).mooddiff_PoverN = mooddiff_PoverN(i);
    alldata(i).mooddiff_PoverNeu = mooddiff_PoverNeu(i);
    alldata(i).mooddiff_NoverNeu = mooddiff_NoverNeu(i);
end


save('results_learning.mat', 'alldata')

%% Put into results var
for i=1:length(vertcat(alldata.mooddiff_PoverN))
    
%     results(i).induced_PoverN = alldata(i).induced_PoverN;
%     results(i).induced_PoverNeu = alldata(i).induced_PoverNeu;
%     results(i).induced_NoverNeu = alldata(i).induced_NoverNeu;
    results(i).mooddiff_PoverN = alldata(i).mooddiff_PoverN;
    results(i).mooddiff_PoverNeu = alldata(i).mooddiff_PoverNeu;
    results(i).mooddiff_NoverNeu = alldata(i).mooddiff_NoverNeu;
    
end

 induced_PoverN = mean(happypos1)>mean(happyneg1);
 induced_PoverNeu = mean(happypos1)>mean(happyB1);
 induced_NoverNeu = mean(happyneg1)<mean(happyB1);

mooddiff_PoverN = vertcat(alldata.mooddiff_PoverN);
mooddiff_PoverNeu = vertcat(alldata.mooddiff_PoverNeu);
mooddiff_NoverNeu = vertcat(alldata.mooddiff_NoverNeu);
clearvars -EXCEPT mooddiff_* happy*


%%

disp('Average no of totes unsure (.5) per subject (T1, T2)')
t=vertcat(results.n_unsure);
[sum(t(mooddiff_PoverN>0,:)); sum(t(mooddiff_PoverN<=0,:))]
disp('Average no of basically unsure (.48-.52) per subject (T1, T2)')
disp(mean( vertcat(results.n_unsureish)))
t=vertcat(results.n_unsureish);
[sum(t(mooddiff_PoverN>0,:)); sum(t(mooddiff_PoverN<=0,:))]
% ? Participants susceptible to mood induction more unsure

%%
disp('--- % chose Neg over Neutral (hiprob, loprob)')
[sum(horzcat(results.moodpref_NeuNeg_hiprob)==-1), sum(horzcat(results.moodpref_NeuNeg_hiprob)==0)];
disp( [ sum(horzcat(results.moodpref_NeuNeg_hiprob)==-1) / sum(~isnan(horzcat(results.moodpref_NeuNeg_hiprob))), sum(horzcat(results.moodpref_NeuNeg_loprob)==-1) / sum(~isnan(horzcat(results.moodpref_NeuNeg_loprob))) ] )

disp('--- % chose Pos over Neutral (hiprob, loprob)')
[sum(horzcat(results.moodpref_NeuPos_hiprob)==1), sum(horzcat(results.moodpref_NeuPos_hiprob)==0)];
disp( [sum(horzcat(results.moodpref_NeuPos_hiprob)==1) / sum(~isnan(horzcat(results.moodpref_NeuPos_hiprob))), sum(horzcat(results.moodpref_NeuPos_loprob)==1) / sum(~isnan(horzcat(results.moodpref_NeuPos_loprob))) ] )

disp('--- % chose Pos over Neg (hiprob, loprob)')
[sum(horzcat(results.moodpref_NegPos_hiprob)==1), sum(horzcat(results.moodpref_NegPos_hiprob)==-1)]
disp( [sum(horzcat(results.moodpref_NegPos_hiprob)==1) / sum(~isnan(horzcat(results.moodpref_NegPos_hiprob))), sum(horzcat(results.moodpref_NegPos_loprob)==1) / sum(~isnan(horzcat(results.moodpref_NegPos_loprob))) ] )

disp('--- score (out of +/-4) chose Neg>Neu, Pos>Neu, Pos>Neg:')
disp([nanmean(horzcat(results.allprobs_neg_over_neu)) ,nanmean(horzcat(results.allprobs_pos_over_neu)), nanmean(horzcat(results.allprobs_pos_over_neg))])

%%

disp('Average confidence rating for P>N: ')
mean( vertcat(results.allprefs_PN_conf) )
max( vertcat(results.allprefs_PN_conf) )
min( vertcat(results.allprefs_PN_conf) )



%%
t=vertcat(results.allprobs_neg_over_neu);
[nansum(t(mooddiff_NoverNeu>0,:)); nansum(t(mooddiff_NoverNeu<=0,:))]
  
t= vertcat(results.allprobs_pos_over_neu);
[nansum(t(mooddiff_PoverNeu>0,:)); nansum(t(mooddiff_PoverNeu<=0,:))]

t= vertcat(results.allprobs_pos_over_neg);
[nansum(t(mooddiff_PoverN>0,:)); nansum(t(mooddiff_PoverN<=0,:))]

disp('--- % chose Pos over Neg (hiprob):')
t = vertcat(results.moodpref_NegPos_hiprob);
disp( [nansum(t(mooddiff_PoverN>0,2)); nansum(t(mooddiff_PoverN<=0,2))] )
t = vertcat(results.moodpref_NegPos_loprob);
disp('--- % chose Pos over Neg (loprob):')
disp( [nansum(t(mooddiff_PoverN>0,2)); nansum(t(mooddiff_PoverN<=0,2))] )

%disp( [sum(horzcat(results.moodpref_NegPos_hiprob)==1) / sum(~isnan(horzcat(results.moodpref_NegPos_hiprob))), sum(horzcat(results.moodpref_NegPos_loprob)==1) / sum(~isnan(horzcat(results.moodpref_NegPos_loprob))) ] )


%Confidence
t=vertcat(results.allprefs_PN_conf);
[mean(t(mooddiff_PoverN>0,:)); mean(t(mooddiff_PoverN<=0,:))]

%%
y= vertcat(results.allprefs_PN_conf);
y=y(:,2);
% Contrasts vs Neutral can be in either test block
t2= vertcat(results.allprobs_pos_over_neg);
t2(isnan(t2(:,1))) = t2(~isnan(t2(:,2)),2); % First merge two test blocks 
x=t2(:,1);
conf=vertcat(results.allprefs_PN_conf);
hold off
figure; plot(t2, conf, 'k.', 'MarkerSize',12); xlabel('Total choices (Positive - Negative)'); ylabel('Confidence rating'); hold on
coeffs = polyfit(x, y, 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y)

%%
% P vs N contrast only possible in 2nd test block
x=mooddiff_PoverN'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
%x=mooddiff_PoverN_zscore'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
% Contrasts vs Neutral can be in either test block
%t2=t; t2(isnan(t2(:,1))) = t2(~isnan(t2(:,2)),2); % First merge two test blocks 
%x=mooddiff_NoverNeu; y=t2(:,1); % N>Neu
%x=mooddiff_PoverNeu; y=t2(:,1); % P>Neu
figure;plot(x,y, '.','MarkerSize',12); xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4])

hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y)

y= vertcat(results.allprefs_PN_conf)*100;y=y(:,2);
x= mooddiff_PoverN_zscore';
%x= mooddiff_PoverN';

figure;plot(x,y, '.','MarkerSize',12); xlabel('Mood induced (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y)



%%
%hps=vertcat(results.HPS)

x=hps(~isnan(hps))';
%y=mooddiff_PoverN(~isnan(hps)); % P>N

y=vertcat(results.allprobs_pos_over_neg); y=y(:,2);
y=y(~isnan(hps))';
figure;plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
xlabel('HPS (mean item score)'); ylabel('Preference (-4 to +4)');set(gca, 'YLim', [-4 4]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y)

y= vertcat(results.allprefs_PN_conf)*100;y=y(:,2);
y=y(~isnan(hps))';
figure;plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS (mean item score)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y)
%% Mood bias in confidence ratings


