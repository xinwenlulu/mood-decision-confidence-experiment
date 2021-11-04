load resultsQalltrials.mat
beta = vertcat(resultsQ.beta);
LR = vertcat(resultsQ.LR);
Drift = vertcat(resultsQ.drift);
fval = vertcat(resultsQ.loglikelihood);

wof = [1;-1;1;-1;1;-1;1;-1;1;-1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;1;-1;1;-1;1;-1;1;1;-1;-1;1;-1;1;1;1;-1;1;1;-1;-1;1;-1];

databeta = mean(mean(beta,2));
sdbeta = std(mean(beta,2));
sebeta = sdbeta/sqrt(44);

dataLR = mean(mean(LR,2));
sdLR = std(mean(LR,2));
seLR = sdLR/sqrt(44);

dataDrift = mean(mean(Drift,2));
sdDrift = std(mean(Drift,2));
seDrift = sdDrift/sqrt(44);

alldata = [databeta dataLR dataDrift];
allsd = [sdbeta sdLR sdDrift];
allse = [sebeta seLR seDrift];

disp(alldata)
disp(allsd)

load resultsSSATalltrials.mat
beta = vertcat(resultsSSAT.beta);
LR = vertcat(resultsSSAT.LR);
Drift = vertcat(resultsSSAT.drift);
LRVar = vertcat(resultsSSAT.LRVar);
threshold = vertcat(resultsSSAT.threshold);
fval = vertcat(resultsSSAT.loglikelihood);

databeta = mean(mean(beta,2));
sdbeta = std(mean(beta,2));
sebeta = sdbeta/sqrt(44);

dataLR = mean(mean(LR,2));
sdLR = std(mean(LR,2));
seLR = sdLR/sqrt(44);

dataDrift = mean(mean(Drift,2));
sdDrift = std(mean(Drift,2));
seDrift = sdDrift/sqrt(44);

dataLRVar = mean(mean(LRVar,2));
sdLRVar = std(mean(LRVar,2));
seLRVar = sdLRVar/sqrt(44);

dataThreshold = mean(mean(threshold,2));
sdThreshold = std(mean(threshold,2));
seThreshold = sdThreshold/sqrt(44);

alldata = [databeta dataLR dataDrift dataLRVar dataThreshold];
allsd = [sdbeta sdLR sdDrift sdLRVar sdThreshold];
allse = [sebeta seLR seDrift seLRVar seThreshold];

disp(alldata)
disp(allsd)
% subplot(1,3,1)
% %bar(1:2,databeta)
% xlabel('parameters'); ylabel('parameter values)'); 
% %axis([0 7 0 1]); 
% xlim([0 3]);
% set(gca,'xtick',[1 2],'xticklabel',{'beta1'; 'beta2';}); %set(gca, 'YLim', [-4 4])
% hold on
% errorbar([databeta(1) databeta(2)],sdbeta,'.')
%%
load resultsQmLvsmH.mat
beta = vertcat(resultsQ.beta);
LR = vertcat(resultsQ.LR);
Drift = vertcat(resultsQ.drift);
fval = vertcat(resultsQ.loglikelihood);

wof = [1;-1;1;-1;1;-1;1;-1;1;-1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;1;-1;1;-1;1;-1;1;1;-1;-1;1;-1;1;1;1;-1;1;1;-1;-1;1;-1];

databeta = mean(mean(beta,2));
sdbeta = std(mean(beta,2));
sebeta = sdbeta/sqrt(44);

dataLR = mean(mean(LR,2));
sdLR = std(mean(LR,2));
seLR = sdLR/sqrt(44);

dataDrift = mean(mean(Drift,2));
sdDrift = std(mean(Drift,2));
seDrift = sdDrift/sqrt(44);

alldata = [databeta dataLR dataDrift];
allsd = [sdbeta sdLR sdDrift];
allse = [sebeta seLR seDrift];

disp(alldata)
disp(allsd)

load resultsSSATmLvsmH.mat
beta = vertcat(resultsSSAT.beta);
LR = vertcat(resultsSSAT.LR);
Drift = vertcat(resultsSSAT.drift);
LRVar = vertcat(resultsSSAT.LRVar);
threshold = vertcat(resultsSSAT.threshold);
fval = vertcat(resultsSSAT.loglikelihood);

databeta = mean(mean(beta,2));
sdbeta = std(mean(beta,2));
sebeta = sdbeta/sqrt(44);

dataLR = mean(mean(LR,2));
sdLR = std(mean(LR,2));
seLR = sdLR/sqrt(44);

dataDrift = mean(mean(Drift,2));
sdDrift = std(mean(Drift,2));
seDrift = sdDrift/sqrt(44);

dataLRVar = mean(mean(LRVar,2));
sdLRVar = std(mean(LRVar,2));
seLRVar = sdLRVar/sqrt(44);

dataThreshold = mean(mean(threshold,2));
sdThreshold = std(mean(threshold,2));
seThreshold = sdThreshold/sqrt(44);

alldata = [databeta dataLR dataDrift dataLRVar dataThreshold];
allsd = [sdbeta sdLR sdDrift sdLRVar sdThreshold];
allse = [sebeta seLR seDrift seLRVar seThreshold];

disp(alldata)
disp(allsd)
% subplot(1,3,1)
% %bar(1:2,databeta)
% xlabel('parameters'); ylabel('parameter values)'); 
% %axis([0 7 0 1]); 
% xlim([0 3]);
% set(gca,'xtick',[1 2],'xticklabel',{'beta1'; 'beta2';}); %set(gca, 'YLim', [-4 4])
% hold on
% errorbar([databeta(1) databeta(2)],sdbeta,'.')
%%
load resultsQ11-44.mat
beta = vertcat(resultsQ.beta);
LR = vertcat(resultsQ.LR);
Drift = vertcat(resultsQ.drift);
fval = vertcat(resultsQ.loglikelihood);

wof = [1;-1;1;-1;1;-1;1;-1;1;-1;-1;1;-1;1;-1;1;-1;1;-1;1;-1;1;1;-1;1;-1;1;-1;1;1;-1;-1;1;-1;1;1;1;-1;1;1;-1;-1;1;-1];

databeta = mean(mean(beta,2));
sdbeta = std(mean(beta,2));
sebeta = sdbeta/sqrt(44);

dataLR = mean(mean(LR,2));
sdLR = std(mean(LR,2));
seLR = sdLR/sqrt(44);

dataDrift = mean(mean(Drift,2));
sdDrift = std(mean(Drift,2));
seDrift = sdDrift/sqrt(44);

alldata = [databeta dataLR dataDrift];
allsd = [sdbeta sdLR sdDrift];
allse = [sebeta seLR seDrift];

disp(alldata)
disp(allsd)

load resultsSSAT11-44.mat
beta = vertcat(resultsSSAT.beta);
LR = vertcat(resultsSSAT.LR);
Drift = vertcat(resultsSSAT.drift);
LRVar = vertcat(resultsSSAT.LRVar);
threshold = vertcat(resultsSSAT.threshold);
fval = vertcat(resultsSSAT.loglikelihood);

databeta = mean(mean(beta,2));
sdbeta = std(mean(beta,2));
sebeta = sdbeta/sqrt(44);

dataLR = mean(mean(LR,2));
sdLR = std(mean(LR,2));
seLR = sdLR/sqrt(44);

dataDrift = mean(mean(Drift,2));
sdDrift = std(mean(Drift,2));
seDrift = sdDrift/sqrt(44);

dataLRVar = mean(mean(LRVar,2));
sdLRVar = std(mean(LRVar,2));
seLRVar = sdLRVar/sqrt(44);

dataThreshold = mean(mean(threshold,2));
sdThreshold = std(mean(threshold,2));
seThreshold = sdThreshold/sqrt(44);

alldata = [databeta dataLR dataDrift dataLRVar dataThreshold];
allsd = [sdbeta sdLR sdDrift sdLRVar sdThreshold];
allse = [sebeta seLR seDrift seLRVar seThreshold];

disp(alldata)
disp(allsd)
% subplot(1,3,1)
% %bar(1:2,databeta)
% xlabel('parameters'); ylabel('parameter values)'); 
% %axis([0 7 0 1]); 
% xlim([0 3]);
% set(gca,'xtick',[1 2],'xticklabel',{'beta1'; 'beta2';}); %set(gca, 'YLim', [-4 4])
% hold on
% errorbar([databeta(1) databeta(2)],sdbeta,'.')