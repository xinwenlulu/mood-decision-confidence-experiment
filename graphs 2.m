mLconfidence = horzcat(alldata.mLconfidence);
for i = 1:length(mLconfidence)
    if isnan(mLconfidence(i))
        mLconfidence(i) = 0;
    end
end

mHvLconfidence = horzcat(alldata.mHvLconfidence);
mHvHconfidence = horzcat(alldata.mHvHconfidence);

disp('confidence mL, mHvL, mHvH:')
disp([mean(mLconfidence) ,mean(mHvLconfidence),mean(mHvHconfidence)]);

disp('std confidence mL, mHvL, mHvH')
disp([std(mLconfidence) ,std(mHvLconfidence),std(mHvHconfidence)]);


data=[nanmean(mLconfidence) ,nanmean(mHvLconfidence),nanmean(mHvHconfidence)];
sd = [nanstd(mLconfidence) ,nanstd(mHvLconfidence),nanstd(mHvHconfidence)];
se = sd / sqrt( length(mLconfidence));

figure
%bar(1:3,data)
xlabel('choice'); ylabel('confidence (0 - 1)'); 
axis([0 4 0 1]); 
set(gca,'xtick',[1 2 3],'xticklabel',{'mL'; 'mHvL'; 'mHvH'}); %set(gca, 'YLim', [-4 4])
hold on
errorbar([data(1) data(2) data(3)],se,'.')

%%
confidenceNeu = horzcat(alldata.confidenceNeu);
confidenceNeg = horzcat(alldata.confidenceNeg);
confidenceNeg = confidenceNeg(~isnan(confidenceNeg));
confidencePos = horzcat(alldata.confidencePos);
confidencePos = confidencePos(~isnan(confidencePos));


disp('confidence mL, mHvL, mHvH:')
disp([mean(confidenceNeu) ,nanmean(confidenceNeg),nanmean(confidencePos)]);

disp('std confidence mL, mHvL, mHvH')
disp([std(confidenceNeu) ,nanstd(confidenceNeg),nanstd(confidencePos)]);


data=[nanmean(confidenceNeu) ,nanmean(confidenceNeg),nanmean(confidencePos)];
sd = [nanstd(confidenceNeu) ,nanstd(confidenceNeg),nanstd(confidencePos)];
se = [sd(1) / sqrt( length(confidenceNeu)), sd(2)/sqrt(length(confidenceNeg)), sd(3)/sqrt(length(confidencePos))];

figure
%bar(1:3,data)
xlabel('happiness block'); ylabel('confidence (0 - 1)'); 
axis([0 4 0 1]); 
set(gca,'xtick',[1 2 3],'xticklabel',{'Neutral'; 'Negative'; 'Positive'}); %set(gca, 'YLim', [-4 4])
hold on
errorbar([data(1) data(2) data(3)],se,'.')