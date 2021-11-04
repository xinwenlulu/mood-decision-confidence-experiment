%% Quality control
subjectid = vertcat(alldata.subjectid);
t=vertcat(alldata.choosing_mHvL);
disp('Subjs chose mHvL stims <50% of the time at the end of 1st and 2nd block:')
disp( subjectid(find(t(:,1)<= .5))')
disp ( subjectid(find(t(:,2)<= .5))' )

disp('Subjs chose mHvL stims <50% of the time on average of 2 learning blocks:')
disp( subjectid(find(mean(t,2)<= .5))')

t=vertcat(alldata.choosing_mHvH);
disp('Subjs chose mHvH stim <50% of the time at the end of 1st and 2nd block:')
disp( subjectid(find(t(:,1)<= .5))' )
disp ( subjectid(find(t(:,2)<= .5))' )
disp('Subjs chose mHvH stims <50% of the time on average of 2 learning blocks:')
disp( subjectid(find(mean(t,2)<= .5))')

t=vertcat(alldata.choosing_mL);
disp('Subjs chose mL stim >50% of the time at the end of 1st and 2nd block:')
disp( subjectid(find(t(:,1)>= .5))' )
disp ( subjectid(find(t(:,2)>= .5))' )
disp('Subjs chose mL stims >50% of the time on average of 2 learning blocks:')
disp( subjectid(find(mean(t,2)>= .5))')


disp('Performance: choosing mHvL')
t=mean(vertcat(alldata.choosing_mHvL),2);
disp(mean(t))

disp('Performance: choosing mHvH')
t=mean(vertcat(alldata.choosing_mHvH),2);
disp(mean(t))

disp('Performance: choosing mL')
t=mean(vertcat(alldata.choosing_mL),2);
disp(mean(t))
%excludesubs = [7 10 13] % Repeat offenders manually from above checks
%%
allconfidence = nan(44,0);
subjectid = [];
for i = 1:length(alldata)
    allconfidence(:,end+1) = vertcat(alldata(i).confidence(:,1));
    allconfidence(:,end+1) = vertcat(alldata(i).confidence(:,2));
    subjectid(end + 1) = alldata(i).subjectid;
end

id_std1_std2 = [vertcat(subjectid); vertcat(std(allconfidence(:,1:2:end))); vertcat(std(allconfidence(:,2:2:end)))];
id_std1_std2 = id_std1_std2';

x = id_std1_std2(id_std1_std2(:,2)<=0.05,:);
y = id_std1_std2(id_std1_std2(:,3)<=0.05,:);
z = id_std1_std2((id_std1_std2(:,2)/2+id_std1_std2(:,3)/2)<=0.05);
disp([id_std1_std2(id_std1_std2(:,2)<=0.05,:)]) 
disp([id_std1_std2(id_std1_std2(:,3)<=0.05,:)])
disp(z)
%%
pos = [];
neg = [];
allmLconfidence = horzcat(alldata.mLconfidence);
mLconfidence = nan(size(allmLconfidence,1)/5,size(allmLconfidence,2));
for i = 1:size(allmLconfidence,1)/5
    mLconfidence(i,:) = nanmean(allmLconfidence(5*i-4:5*i,:));
end

allmHvLconfidence = horzcat(alldata.mHvLconfidence);
mHvLconfidence = nan(size(allmHvLconfidence,1)/5,size(allmHvLconfidence,2));
for i = 1:size(allmHvLconfidence,1)/5
    mHvLconfidence(i,:) = nanmean(allmHvLconfidence(5*i-4:5*i,:));
end

allmHvHconfidence = horzcat(alldata.mHvHconfidence);
mHvHconfidence = nan(size(allmHvHconfidence,1)/5,size(allmHvHconfidence,2));
for i = 1:size(allmHvHconfidence,1)/5
    mHvHconfidence(i,:) = nanmean(allmHvHconfidence(5*i-4:5*i,:));
end

for i = 1:length(alldata)
    
    if alldata(i).wof_order(2) == 1
        pos(end+1) = i;
    else
        neg(end+1) = i;
    end
    
end

neutral = 1:2:size(allmLconfidence,2);
positive = 2.*pos;
negative = 2.*neg;

confidenceNeumLall = allmLconfidence(:,neutral);
confidenceNeumL = nan(size(confidenceNeumLall,1)/5,size(confidenceNeumLall,2));
for i = 1:size(confidenceNeumLall,1)/5
    confidenceNeumL(i,:) = nanmean(confidenceNeumLall(5*i-4:5*i,:));
end
    
confidencePosmLall = allmLconfidence(:,positive);
confidencePosmL = nan(size(confidencePosmLall,1)/5,size(confidencePosmLall,2));
for i = 1:size(confidencePosmLall,1)/5
    confidencePosmL(i,:) = nanmean(confidencePosmLall(5*i-4:5*i,:));
end

confidenceNegmLall = allmLconfidence(:,negative);
confidenceNegmL = nan(size(confidenceNegmLall,1)/5,size(confidenceNegmLall,2));
for i = 1:size(confidenceNegmLall,1)/5
    confidenceNegmL(i,:) = nanmean(confidenceNegmLall(5*i-4:5*i,:));
end


subplot(2,3,1); hold on;  %color order
mLneutralplot = nanmean(confidenceNeumL,2);
mLposplot = nanmean(confidencePosmL,2);
mLnegplot = nanmean(confidenceNegmL,2);

Avgplot = nanmean(mLconfidence,2);

    hold on
    plot(Avgplot, 'k', 'LineWidth',2); 
    plot(mLneutralplot,'color','b'); 
    plot(mLnegplot,'color','r');
    plot(mLposplot,'color','g'); 
ylim([0 1]);
set(gca,'xtick',[1 2 3 4 5 6 7 8],'xticklabel',{'1-5','6-10','11-15','16-20','21-25','26-30','31-35','36-40'}); 
xlabel('trial number'); ylabel('confidence mL');
legend('Average', 'mL neutral', 'mL negative', 'mL positive');

confidenceNeumHvLall = allmHvLconfidence(:,neutral);
confidenceNeumHvL = nan(size(confidenceNeumHvLall,1)/5,size(confidenceNeumHvLall,2));
for i = 1:size(confidenceNeumHvLall,1)/5
    confidenceNeumHvL(i,:) = nanmean(confidenceNeumHvLall(5*i-4:5*i,:));
end
confidencePosmHvLall = allmHvLconfidence(:,positive);
confidencePosmHvL = nan(size(confidencePosmHvLall,1)/5,size(confidencePosmHvLall,2));
for i = 1:size(confidencePosmHvLall,1)/5
    confidencePosmHvL(i,:) = nanmean(confidencePosmHvLall(5*i-4:5*i,:));
end

confidenceNegmHvLall = allmHvLconfidence(:,negative);
confidenceNegmHvL = nan(size(confidenceNegmHvLall,1)/5,size(confidenceNegmHvLall,2));
for i = 1:size(confidenceNegmHvLall,1)/5
    confidenceNegmHvL(i,:) = nanmean(confidenceNegmHvLall(5*i-4:5*i,:));
end

subplot(2,3,2); hold on; mike={'b','y','c','m'}; %color order
mHvLneutralplot = nanmean(confidenceNeumHvL,2);
mHvLposplot = nanmean(confidencePosmHvL,2);
mHvLnegplot = nanmean(confidenceNegmHvL,2);

Avgplot = nanmean(mHvLconfidence,2);

    hold on
    plot(Avgplot, 'k', 'LineWidth',2); 
    plot(mHvLneutralplot,'color','b'); 
    plot(mHvLnegplot,'color','r');
    plot(mHvLposplot,'color','g'); 
ylim([0 1]);

xlabel('trial number'); ylabel('confidence mHvL');
set(gca,'xtick',[1 2 3 4],'xticklabel',{'1-5','6-10','11-15','16-20'}); 
legend('Average', 'mHvL neutral', 'mHvL negative', 'mHvL positive');


confidenceNeumHvHall = allmHvHconfidence(:,neutral);
confidenceNeumHvH = nan(size(confidenceNeumHvHall,1)/5,size(confidenceNeumHvHall,2));
for i = 1:size(confidenceNeumHvHall,1)/5
    confidenceNeumHvH(i,:) = nanmean(confidenceNeumHvHall(5*i-4:5*i,:));
end
confidencePosmHvHall = allmHvHconfidence(:,positive);
confidencePosmHvH = nan(size(confidencePosmHvHall,1)/5,size(confidencePosmHvHall,2));
for i = 1:size(confidencePosmHvHall,1)/5
    confidencePosmHvH(i,:) = nanmean(confidencePosmHvHall(5*i-4:5*i,:));
end

confidenceNegmHvHall = allmHvHconfidence(:,negative);
confidenceNegmHvH = nan(size(confidenceNegmHvHall,1)/5,size(confidenceNegmHvHall,2));
for i = 1:size(confidenceNegmHvHall,1)/5
    confidenceNegmHvH(i,:) = nanmean(confidenceNegmHvHall(5*i-4:5*i,:));
end

subplot(2,3,3); hold on; mike={'b','y','c','m'}; %color order
mHvHneutralplot = nanmean(confidenceNeumHvH,2);
mHvHposplot = nanmean(confidencePosmHvH,2);
mHvHnegplot = nanmean(confidenceNegmHvH,2);

Avgplot = nanmean(mHvHconfidence,2);

    hold on
    plot(Avgplot, 'k', 'LineWidth',2); 
    plot(mHvHneutralplot,'color','b'); 
    plot(mHvHnegplot,'color','r');
    plot(mHvHposplot,'color','g'); 
ylim([0 1]);

xlabel('trial number'); ylabel('confidence mHvH');
set(gca,'xtick',[1 2 3 4],'xticklabel',{'1-5','6-10','11-15','16-20'});
legend('Average', 'mHvH neutral', 'mHvH negative', 'mHvH positive');


%%
mLconfidenceNeu = horzcat(mLneutralplot(1:20));
mLconfidenceNeg = horzcat(mLnegplot(1:20));
mLconfidencePos = horzcat(mLposplot(1:20));

disp('confidence mL neutral negative positive:')
disp([nanmean(mLconfidenceNeu),nanmean(mLconfidenceNeg),nanmean(mLconfidencePos)]);

disp('std confidence mL neutral, negative,')
disp([nanstd(mLconfidenceNeu) ,nanstd(mLconfidenceNeg),nanstd(mLconfidencePos)]);


data=[nanmean(mLconfidenceNeu) ,nanmean(mLconfidenceNeg),nanmean(mLconfidencePos)];
sd = [nanstd(mLconfidenceNeu) ,nanstd(mLconfidenceNeg),nanstd(mLconfidencePos)];
se = sd/length(mLconfidenceNeu);

figure
%bar(1:3,data)
subplot(2,3,1)
xlabel('happiness blocks'); ylabel('confidence (0 - 1)'); 
title('mL confidence')
axis([0 4 0 1]); 
set(gca,'xtick',[1 2 3],'xticklabel',{'neutral'; 'negative'; 'positive'}); %set(gca, 'YLim', [-4 4])
hold on
errorbar([data(1) data(2) data(3)],se,'.')


mHvLconfidenceNeu = horzcat(mHvLneutralplot(1:10));
mHvLconfidenceNeg = horzcat(mHvLnegplot(1:10));
mHvLconfidencePos = horzcat(mHvLposplot(1:10));

disp('confidence mHvL neutral negative positive:')
disp([nanmean(mHvLconfidenceNeu),nanmean(mHvLconfidenceNeg),nanmean(mHvLconfidencePos)]);

disp('std confidence mHvL neutral, negative,')
disp([nanstd(mHvLconfidenceNeu) ,nanstd(mHvLconfidenceNeg),nanstd(mHvLconfidencePos)]);


data1=[nanmean(mHvLconfidenceNeu) ,nanmean(mHvLconfidenceNeg),nanmean(mHvLconfidencePos)];
sd1 = [nanstd(mHvLconfidenceNeu) ,nanstd(mHvLconfidenceNeg),nanstd(mHvLconfidencePos)];
se1 = sd/length(mHvLconfidenceNeu);

%bar(1:3,data)
subplot(2,3,2)
xlabel('happiness blocks'); ylabel('confidence (0 - 1)'); 
title('mHvL confidence')
axis([0 4 0 1]); 
set(gca,'xtick',[1 2 3],'xticklabel',{'neutral'; 'negative'; 'positive'}); %set(gca, 'YLim', [-4 4])
hold on
errorbar([data1(1) data1(2) data1(3)],se1,'.')

mHvHconfidenceNeu = horzcat(mHvHneutralplot(1:10));
mHvHconfidenceNeg = horzcat(mHvHnegplot(1:10));
mHvHconfidencePos = horzcat(mHvHposplot(1:10));

disp('confidence mHvH neutral negative positive:')
disp([nanmean(mHvHconfidenceNeu),nanmean(mHvHconfidenceNeg),nanmean(mHvHconfidencePos)]);

disp('std confidence mHvH neutral, negative,')
disp([nanstd(mHvHconfidenceNeu) ,nanstd(mHvHconfidenceNeg),nanstd(mHvHconfidencePos)]);


data1=[nanmean(mHvHconfidenceNeu) ,nanmean(mHvHconfidenceNeg),nanmean(mHvHconfidencePos)];
sd1 = [nanstd(mHvHconfidenceNeu) ,nanstd(mHvHconfidenceNeg),nanstd(mHvHconfidencePos)];
se1 = sd/length(mHvHconfidenceNeu);

%bar(1:3,data)
subplot(2,3,3)
xlabel('happiness blocks'); ylabel('confidence (0 - 1)'); 
title('mHvH confidence')
axis([0 4 0 1]); 
set(gca,'xtick',[1 2 3],'xticklabel',{'neutral'; 'negative'; 'positive'}); %set(gca, 'YLim', [-4 4])
hold on
errorbar([data1(1) data1(2) data1(3)],se1,'.')

%%
hapchangeNeg = horzcat(alldata.hapchangeNeg);
x=hapchangeNeg(~isnan(hapchangeNeg));


neg_over_neu = vertcat(results.allprobs_neg_over_neu);
y=neg_over_neu(~isnan(hapchangeNeg)); % P>N

%oldhps=vertcat(Olddata.HPS)
%x=[x, oldhps(~isnan(oldhps))'];
%y=[y, mooddiff_PoverNeuOLD(~isnan(oldhps))];
figure('color',[1 1 1]); subplot(2,2,1); 
plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y', 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('Mood change'); ylabel('Mood bias binary(Neg>Neu)');%set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Neg > Neu: r=%.02f, p=%.02f\n',r(2), p(2));

hapchangeNeg = horzcat(alldata.hapchangeNeg);
x=hapchangeNeg(~isnan(hapchangeNeg));


neg_over_neu = vertcat(results.conf_Neg_over_Neu);
y=neg_over_neu(~isnan(hapchangeNeg))'; 


%oldhps=vertcat(Olddata.HPS)
%x=[x, oldhps(~isnan(oldhps))'];
%y=[y, mooddiff_PoverNeuOLD(~isnan(oldhps))];
subplot(2,2,2); 
plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('Mood change'); ylabel('Mood bias confidence(Neg>Neu)');%set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Neg > Neu: r=%.02f, p=%.02f\n',r(2), p(2));


hapchangePos = horzcat(alldata.hapchangePos);
x=hapchangePos(~isnan(hapchangePos));


pos_over_neu = vertcat(results.allprobs_pos_over_neu);
y=pos_over_neu(~isnan(hapchangePos)); % P>N


subplot(2,2,3); 
plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y', 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('Mood change'); ylabel('Mood bias binary(Pos>Neu)');%set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Pos > Neu: r=%.02f, p=%.02f\n',r(2), p(2));

hapchangePos = horzcat(alldata.hapchangePos);
x=hapchangePos(~isnan(hapchangePos));


pos_over_neu = vertcat(results.conf_Pos_over_Neu);
y=pos_over_neu(~isnan(hapchangePos))'; 


%oldhps=vertcat(Olddata.HPS)
%x=[x, oldhps(~isnan(oldhps))'];
%y=[y, mooddiff_PoverNeuOLD(~isnan(oldhps))];
subplot(2,2,4); 
plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('Mood change'); ylabel('Mood bias confidence(Pos>Neu)');%set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Pos > Neu: r=%.02f, p=%.02f\n',r(2), p(2));

%%
hapchangeNeg = horzcat(alldata.hapchangeNeg);
y=hapchangeNeg(~isnan(hapchangeNeg));


hps_neg = vertcat(results.hps_neg);
x=hps_neg(~isnan(hapchangeNeg)); % P>N

%oldhps=vertcat(Olddata.HPS)
%x=[x, oldhps(~isnan(oldhps))'];
%y=[y, mooddiff_PoverNeuOLD(~isnan(oldhps))];
figure('color',[1 1 1]); subplot(2,2,1); 
plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x', y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Mood change');%set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
title('Negative Wof')
[r, p]=corrcoef(x,y); fprintf('HPS and Mood change: r=%.02f, p=%.02f\n',r(2), p(2));


hapchangePos = horzcat(alldata.hapchangePos);
y=hapchangePos(~isnan(hapchangePos));


hps_pos = vertcat(results.hps_pos);
x=hps_pos(~isnan(hapchangePos)); % P>N


subplot(2,2,2); 
plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x', y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Mood Change');%set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
title('Positive Wof')
[r, p]=corrcoef(x,y); fprintf('Mood change: r=%.02f, p=%.02f\n',r(2), p(2));

%%


hapchangePos = vertcat(alldata.hapchangePos);
hapchangeNeg = vertcat(alldata.hapchangeNeg);
hapchange = hapchangePos;
hapchange(isnan(hapchange))= hapchangeNeg(~isnan(hapchangeNeg));
y=hapchange;


hps_pos = vertcat(results.hps_pos);
hps_neg = vertcat(results.hps_neg);
hps = hps_pos;
hps(isnan(hps)) = hps_neg(~isnan(hps_neg));
x=hps; % P>N


subplot(1,1,1); 
plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Mood Change');%set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
title('all Wof')
[r, p]=corrcoef(x,y); fprintf('Mood change: r=%.02f, p=%.02f\n',r(2), p(2));

%%

hps_neg = vertcat(results.hps_neg);
x=hps_neg(~isnan(hps_neg)); % P>N


neg_over_neu = vertcat(results.allprobs_neg_over_neu);
y=neg_over_neu(~isnan(hps_neg)); % P>N

%oldhps=vertcat(Olddata.HPS)
%x=[x, oldhps(~isnan(oldhps))'];
%y=[y, mooddiff_PoverNeuOLD(~isnan(oldhps))];
figure('color',[1 1 1]); subplot(2,2,1); 
plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Mood bias binary(Neg>Neu)');%set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Neg > Neu: r=%.02f, p=%.02f\n',r(2), p(2));


hps_neg = vertcat(results.hps_neg);
x=hps_neg(~isnan(hps_neg)); % P>N

neg_over_neu = vertcat(results.conf_Neg_over_Neu);
y=neg_over_neu(~isnan(hps_neg)); 


%oldhps=vertcat(Olddata.HPS)
%x=[x, oldhps(~isnan(oldhps))'];
%y=[y, mooddiff_PoverNeuOLD(~isnan(oldhps))];
subplot(2,2,2); 
plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Mood bias confidence(Neg>Neu)');%set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Neg > Neu: r=%.02f, p=%.02f\n',r(2), p(2));


hps_pos = vertcat(results.hps_pos);
x=hps_pos(~isnan(hps_pos)); % P>N

pos_over_neu = vertcat(results.allprobs_pos_over_neu);
y=pos_over_neu(~isnan(hps_pos)); % P>N


subplot(2,2,3); 
plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Mood bias binary(Pos>Neu)');%set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Pos > Neu: r=%.02f, p=%.02f\n',r(2), p(2));

hps_pos = vertcat(results.hps_pos);
x=hps_pos(~isnan(hps_pos)); % P>N


pos_over_neu = vertcat(results.conf_Pos_over_Neu);
y=pos_over_neu(~isnan(hps_pos)); 


%oldhps=vertcat(Olddata.HPS)
%x=[x, oldhps(~isnan(oldhps))'];
%y=[y, mooddiff_PoverNeuOLD(~isnan(oldhps))];
subplot(2,2,4); 
plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Mood bias confidence(Pos>Neu)');%set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Pos > Neu: r=%.02f, p=%.02f\n',r(2), p(2));
