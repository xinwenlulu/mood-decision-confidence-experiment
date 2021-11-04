%%% Analyse test

    
% Consistency between binary and ratings

% Confirm learning of low vs high stim (within-blocks)

% Confirm learning of low vs high stim (between-blocks)

% Check for mood bias (between-blocks) in binary

% Check for mood bias (between-blocks) in ratings

%for isub=1:length(files)

%taskverstring = strrep(taskverstring,'_', ' ');

for isub=1:length(alldata)
    
    %load(files(isub).name);
    %if isempty(data.WoF{1})
    %eval(['tests = ', testsvar, ';']);
    tests=alldata(isub).test;
    results(isub).subjectid = alldata(isub).subjectid;
   
    
    for itest=1:size(tests,2)
        
        clearvars -EXCEPT alldata* mooddiff_* Old* moodpref* tests isub itest results taskverstring
        
        trialmtx = tests{itest}(:,2:end); halfway=size(trialmtx,1)/2; % Only take first half (binary choice)

        diffEV = trialmtx(:,3) ~= trialmtx(:,4);
        sameEV = trialmtx(:,3) == trialmtx(:,4);
        withinblock = trialmtx(:,5) == trialmtx(:,6);
        betweenblock = trialmtx(:,5) ~= trialmtx(:,6);

        % Confidence
        binaryresp = trialmtx(1:halfway,11);
        for trial = 1:size(binaryresp,1)
            if binaryresp(trial) == trialmtx(trial,1)
                binaryresp(trial) = 1;
            else
                binaryresp(trial) = 2;
            end
        end
        confidence = trialmtx(halfway+1:end,10);
        choseleft_ratedleft = confidence(binaryresp==1)<.5;
        choseleft_ratedleft = and(confidence<.5,binaryresp==1);
        choseright_ratedright = confidence(binaryresp==2)>.5;
        choseright_ratedright = and(confidence>.5, binaryresp==2);
        
        consistent_binaryandconfidence = or(and(confidence<.48,binaryresp==1),and(confidence>.52, binaryresp==2)) ;
        results(isub).confidenceconcordance_all(itest) = (sum(choseleft_ratedleft)+sum(choseright_ratedright)) / length(confidence);
        results(isub).confidenceconcordance_notindiff(itest) = sum(consistent_binaryandconfidence)/length(consistent_binaryandconfidence);
        %if sum(consistent_binaryandconfidence)/length(consistent_binaryandconfidence)==0
        %    disp( [isub, itest] ); pause; disp([binaryresp, confidence]); pause
        %end

        results(isub).n_unsure(itest) = sum(confidence==.5);
        results(isub).n_unsureish(itest) = sum( confidence>.48 & confidence<.52 );
         
        for itrial=1:halfway
            resp =trialmtx(itrial,11);
            if resp == trialmtx(itrial,1)
                resp = 1;
            else
                resp = 2;
            end
            stimprobchosen(itrial) = trialmtx(itrial,resp+2); %condition chosen
            stimprobnotchosen(itrial) = trialmtx(itrial,5-resp); %condition not chosen
            wofpref(itrial) = trialmtx(itrial,4+resp);
            wofnonpref(itrial) = trialmtx(itrial,7-resp);
            
            
            if confidence(itrial) ~= .5
                if confidence(itrial) < .5 
                    equiv_LR = 1;
                elseif confidence(itrial) > .5 
                    equiv_LR = 2;
                end
                
                wofpref_conf(itrial) = trialmtx(halfway+itrial,4+equiv_LR);
                stimprobchosen_conf(itrial) = trialmtx(halfway+itrial,2+equiv_LR); 
                stimprobnotchosen_conf(itrial) = trialmtx(halfway+itrial,5-equiv_LR);
            else
                wofpref_conf(itrial) = nan;
                stimprobchosen_conf(itrial) = nan;
                stimprobnotchosen_conf(itrial) = nan;
            end
        end
        
        
        confidence(confidence(:)<0.5) = 1 - confidence(confidence(:)<0.5);
        confidence(:) = 2.*(confidence(:) - 0.5);
        con_sameEV = confidence(sameEV(10:18));
        wofpref_con_sameEV = wofpref_conf(sameEV(10:18));
        
        hps = alldata(isub).HPS;
        if sum(wofpref_conf(wofpref_conf == 1))==0
            results(isub).conf_Neg_over_Neu(itest) = (sum(con_sameEV(wofpref_con_sameEV == -1))-sum(con_sameEV(wofpref_con_sameEV == 0)))/3;
            neg_over_neu_con2binary = sum(wofpref_con_sameEV == -1) - sum(wofpref_con_sameEV == 0);
            results(isub).hps_neg = hps;
            results(isub).conf_Pos_over_Neu(itest) = nan;
            results(isub).hps_pos = nan;
        else
            results(isub).conf_Pos_over_Neu(itest) = (sum(con_sameEV(wofpref_con_sameEV == 1))-sum(con_sameEV(wofpref_con_sameEV == 0)))/3;
            pos_over_neu_con2binary = sum(wofpref_con_sameEV == 1) - sum(wofpref_con_sameEV == 0);
            results(isub).conf_Neg_over_Neu(itest) = nan;
            results(isub).hps_pos = hps;
            results(isub).hps_neg = nan;
        end
        
        wofpref_merged = [wofpref,wofpref_conf];
        stimprobchosen_merged = [stimprobchosen,stimprobchosen_conf];
        stimprobnotchosen_merged = [stimprobnotchosen,stimprobnotchosen_conf];
            
        
        results(isub).testall_chosehiprob(itest) = sum(stimprobnotchosen == 1) / length(stimprobchosen);
        results(isub).testall_chosehiprob_conf(itest) = nansum(stimprobnotchosen_conf == 1) / length(stimprobchosen_conf);
        results(isub).testall_chosehiprob_merged(itest) = nansum(stimprobnotchosen_merged== 1) / length(stimprobnotchosen_merged(~isnan(stimprobnotchosen_merged)));
        
        testlearning_binary = diffEV;
        testlearning_binary = testlearning_binary(1:halfway);
        testlearning_binary_chose = stimprobchosen(testlearning_binary);
        testlearning_binary_notchose = stimprobnotchosen(testlearning_binary);
        results(isub).testlearning_chosehiprob(itest) = sum(testlearning_binary_notchose == 1) /(length(testlearning_binary_chose)-2);
        results(isub).testlearning_chosehiprob_conf(itest) = nansum( stimprobnotchosen_conf(testlearning_binary) == 1) /( length(testlearning_binary_chose)-2);
        results(isub).testlearning_chosehiprob_merged(itest) = nansum( stimprobnotchosen_merged(diffEV)==1) / (length(stimprobnotchosen_merged(~isnan(stimprobnotchosen_merged)))-10) ;
        
        testmood_binary = sameEV;
        testmood_binary = testmood_binary(1:halfway);

        % Identify test comparisons
        choice_neutral = or( trialmtx(1:halfway,5) == 0, trialmtx(1:halfway,6) == 0 ); 
        choice_pos = or( trialmtx(1:halfway,5) == 1,  trialmtx(1:halfway,6) == 1 );
        choice_neg = or( trialmtx(1:halfway,5) == -1,  trialmtx(1:halfway,6) == -1 ); 

        choice_neutralvspos = and(choice_neutral, choice_pos);
        choice_neutralvsneg = and(choice_neutral, choice_neg);
        %choice_negvspos = and(choice_neg, choice_pos);

        % Test mood bias for binary choices
        NeuNegprefs = wofpref( and( choice_neutralvsneg, sameEV(1:halfway)) );
        NeuNegprobs = stimprobchosen( and( choice_neutralvsneg, sameEV(1:halfway)) );
        NeuNegprefs_constpref = wofpref( and(and( choice_neutralvsneg, sameEV(1:halfway)),consistent_binaryandconfidence) );
        NeuNegprobs_constpref = stimprobchosen( and( and( choice_neutralvsneg, sameEV(1:halfway)), consistent_binaryandconfidence) );
                
        if ~isempty( NeuNegprefs( NeuNegprobs== .6) ), results(isub).moodpref_NeuNeg_hiprob(itest) = NeuNegprefs( NeuNegprobs== .6); 
        else results(isub).moodpref_NeuNeg_hiprob(itest) = nan; end
        
        if ~isempty( NeuNegprefs( NeuNegprobs== .2) ), results(isub).moodpref_NeuNeg_loprob(itest) = NeuNegprefs( NeuNegprobs== .2); 
        else results(isub).moodpref_NeuNeg_loprob(itest) = nan; end

        NeuPosprefs = wofpref( and( choice_neutralvspos, sameEV(1:halfway)) );
        NeuPosprobs = stimprobchosen( and( choice_neutralvspos, sameEV(1:halfway)) );
        NeuPosprefs_constpref = wofpref( and( and( choice_neutralvspos, sameEV(1:halfway)), consistent_binaryandconfidence) );
        NeuPosprobs_constpref = stimprobchosen( and( and( choice_neutralvspos, sameEV(1:halfway)), consistent_binaryandconfidence ) );
        
        if ~isempty( NeuPosprefs( NeuPosprobs== .6) ), results(isub).moodpref_NeuPos_hiprob(itest) = NeuPosprefs( NeuPosprobs== .6); 
        else results(isub).moodpref_NeuPos_hiprob(itest) = nan; end
        if ~isempty( NeuPosprefs( NeuPosprobs== .2) ), results(isub).moodpref_NeuPos_loprob(itest) = NeuPosprefs( NeuPosprobs== .2); 
        else results(isub).moodpref_NeuPos_loprob(itest) = nan; end

        %NegPosprefs = wofpref( and( choice_negvspos, sameEV(1:halfway)) );
        %NegPosprobs = stimprobchosen( and( choice_negvspos, sameEV(1:halfway)) );
        %NegPosprefs_constpref = wofpref( and( and( choice_negvspos, sameEV(1:halfway)), consistent_binaryandconfidence ) );
        %NegPosprobs_constpref = stimprobchosen( and( and( choice_negvspos, sameEV(1:halfway)), consistent_binaryandconfidence) );
        
%         if ~isempty( NegPosprefs( NegPosprobs== .6) ), results(isub).moodpref_NegPos_hiprob(itest) = NegPosprefs( NegPosprobs== .6); 
%         else results(isub).moodpref_NegPos_hiprob(itest) = nan; end
%         if ~isempty( NegPosprefs( NegPosprobs== .2) ), results(isub).moodpref_NegPos_loprob(itest) = NegPosprefs( NegPosprobs== .2); 
%         else results(isub).moodpref_NegPos_loprob(itest) = nan; end
        
    
        allprefs_NB = wofpref( choice_neutralvsneg );
        allprefs_NB_constpref = wofpref( and(choice_neutralvsneg,consistent_binaryandconfidence) );
        
        if ~isempty( allprefs_NB ), 
            results(isub).allprobs_neg_over_neu(itest) = sum( allprefs_NB == -1) - sum( allprefs_NB == 0) + neg_over_neu_con2binary ; 
            results(isub).allprobs_constpref_neg_over_neu(itest) = sum( allprefs_NB_constpref == -1) - sum( allprefs_NB_constpref == 0); 
            
        else results(isub).allprobs_neg_over_neu(itest) = nan; end
        
        allprefs_PB = wofpref( choice_neutralvspos );
        if ~isempty( allprefs_PB ), results(isub).allprobs_pos_over_neu(itest) = sum( allprefs_PB == 1) - sum( allprefs_PB == 0)+ pos_over_neu_con2binary; 
        else results(isub).allprobs_pos_over_neu(itest) = nan; end
        
%         allprefs_PN = wofpref( choice_negvspos );
%         if ~isempty( allprefs_PN ), results(isub).allprobs_pos_over_neg(itest) = sum( allprefs_PN == 1) - sum( allprefs_PN == -1); 
%         else results(isub).allprobs_pos_over_neg(itest) = nan; end
        
        % Mood effect in confidence
        % Want ratings to vary from 0 (strongly prefer neg) to 1 (strongly prefer pos)
        % So reverse rating if rating<.5 but wofpref is +1, or if rating>.5 but wofpref is -1
        
        
        recode= or ( and( wofpref_conf==-1,confidence'>.5 ), and( wofpref_conf==1,confidence'<.5 ) );
        allprefs_conf_recode = confidence;
        allprefs_conf_recode(recode) = 1-allprefs_conf_recode(recode);
        
        
        %P > N
%         recode1= and( wofpref_conf(choice_negvspos)==-1,confidence(choice_negvspos)'>.5 ); % Prefer -1 but rating >.5
%         recode2= and( wofpref_conf(choice_negvspos)==1,confidence(choice_negvspos)'<.5 ); % Prefer +1 but rating <.5
%         recode= or(recode1,recode2);
%         allprefs_PN_conf =confidence(choice_negvspos);
%         allprefs_PN_conf(recode) = 1-allprefs_PN_conf(recode);
%         if ~isempty( allprefs_PN_conf ), results(isub).allprefs_PN_conf(itest) = mean(allprefs_PN_conf);
%         else results(isub).allprefs_PN_conf(itest) = nan; end
        
        
%         %P > N conf as binary
%         for k=1:length(allprefs_PN_conf)
%             if allprefs_PN_conf(k) >= .52
%                 conf_PN_binary(k) = 1;
%             elseif allprefs_PN_conf(k) <= .48
%                 conf_PN_binary(k) = -1;
%             else
%                 conf_PN_binary(k) = 0;
%             end
%         end
%         if ~isempty( allprefs_PN_conf ), results(isub).conf_PN_binary(itest) = sum(conf_PN_binary);
%         else results(isub).conf_PN_binary(itest) = nan; end
                
        %need to make binary a decimal
%         %-4 = 0, while 4 = 100
%         tmp = results(isub).allprobs_pos_over_neg(itest); 
%         tmp=tmp + length(allprefs_PN); 
%         if ~isempty( allprefs_PN_conf ), results(isub).allprobs_PN_pcnt(itest) = tmp*.125;
%         else results(isub).allprobs_PN_pcnt(itest) = nan; end;
%         if ~isempty( allprefs_PN_conf ), results(isub).allprefs_PN_merged(itest) = mean([results(isub).allprobs_PN_pcnt(itest); results(isub).allprefs_PN_conf(itest)])
%         else results(isub).allprefs_PN_merged(itest) = nan; end;
        
%         %Same EV
%         recode1= and( wofpref_conf( and(sameEV(halfway+1:end),choice_negvspos) ) == -1, confidence( and(sameEV(halfway+1:end),choice_negvspos) )'>.5 ); % Prefer -1 but rating >.5
%         recode2= and( wofpref_conf( and(sameEV(halfway+1:end),choice_negvspos) ) == 1, confidence( and(sameEV(halfway+1:end),choice_negvspos) )'<.5 ); % Prefer +1 but rating <.5
%         recode= or(recode1,recode2);
%         allprefs_PN_conf_sameEV =confidence( and(sameEV(halfway+1:end),choice_negvspos) );
%         allprefs_PN_conf_sameEV(recode) = 1-allprefs_PN_conf_sameEV(recode);
%         %Diff EV
%         allprefs_PN_conf_loprob = allprefs_PN_conf_sameEV( trialmtx( and(sameEV(halfway+1:end),choice_negvspos),1)==.2 );
%         allprefs_PN_conf_hiprob = allprefs_PN_conf_sameEV( trialmtx( and(sameEV(halfway+1:end),choice_negvspos),1)==.6 );
%         if ~isempty( allprefs_PN_conf_loprob), results(isub).allprefs_PN_conf_loprob(itest) = allprefs_PN_conf_loprob; end
%         if ~isempty( allprefs_PN_conf_hiprob), results(isub).allprefs_PN_conf_hiprob(itest) = allprefs_PN_conf_hiprob; end
        
                
    end

end


%% Quality control
subjectid = vertcat(results.subjectid);
t=vertcat(results.testlearning_chosehiprob);
disp('Subjs unable to identify mH stims >50% of the time (binary):')
disp( subjectid(find(t(:,1)<= .5))')
%disp ( find(t(:,2)<= .5)' )

t=vertcat(results.testlearning_chosehiprob_merged);
disp('Subjs unable to identify mH stim >50% of the time (binary+confidence):')
disp( subjectid(find(t(:,1)<= .5))' )
%disp ( find(t(:,2)<= .5)' )

t=vertcat(results.testlearning_chosehiprob_conf);
%t=vertcat(results.testall_chosehiprob_conf);
disp('Subjs unable to identify mH stim >50% of the time (confidence):')
disp( subjectid(find(t(:,1)<= .5))' )
%disp ( find(t(:,2)<= .5)' )

disp('Concordance of binary choice with coonfidence rating: ')
disp(mean( vertcat(results.confidenceconcordance_all)))

disp('Concordance of binary choice with coonfidence rating when not close to .5: ')
disp(mean( vertcat(results.confidenceconcordance_notindiff)))

disp('Average no of totally unsure (confidence =.5) per subject')
disp(mean( vertcat(results.n_unsure)))
disp('Average no of basically unsure (.48-.52) per subject')
disp(mean( vertcat(results.n_unsureish)))

disp('Subjs rating unsure more than half the time:')
t=vertcat(results.n_unsureish);
disp( subjectid(find(t(:,1)> 4))' ) ; %disp( find(t(:,2)> 5) );

disp('Subjs with <=50% concordance in ratings:')
t=vertcat(results.confidenceconcordance_all);
disp( subjectid(find(t(:,1)<= .5))'); %disp( find(t(:,2)<= .5) )

disp('Performance: mean test learnig accuracy')
t=vertcat(results.testlearning_chosehiprob_merged);
disp(mean(t))

%excludesubs = [7 10 13] % Repeat offenders manually from above checks

%%
% disp('--- % chose Neg over Neutral (hiprob, loprob)')
% [sum(horzcat(results.moodpref_NeuNeg_hiprob)==-1), sum(horzcat(results.moodpref_NeuNeg_hiprob)==0)]
% disp( [ sum(horzcat(results.moodpref_NeuNeg_hiprob)==-1) / sum(~isnan(horzcat(results.moodpref_NeuNeg_hiprob))), sum(horzcat(results.moodpref_NeuNeg_loprob)==-1) / sum(~isnan(horzcat(results.moodpref_NeuNeg_loprob))) ] )
% 
% disp('--- % chose Pos over Neutral (hiprob, loprob)')
% [sum(horzcat(results.moodpref_NeuPos_hiprob)==1), sum(horzcat(results.moodpref_NeuPos_hiprob)==0)]
% disp( [sum(horzcat(results.moodpref_NeuPos_hiprob)==1) / sum(~isnan(horzcat(results.moodpref_NeuPos_hiprob))), sum(horzcat(results.moodpref_NeuPos_loprob)==1) / sum(~isnan(horzcat(results.moodpref_NeuPos_loprob))) ] )
% 
% disp('--- % chose Pos over Neg (hiprob, loprob)')
% [sum(horzcat(results.moodpref_NegPos_hiprob)==1), sum(horzcat(results.moodpref_NegPos_hiprob)==-1)]
% disp( [sum(horzcat(results.moodpref_NegPos_hiprob)==1) / sum(~isnan(horzcat(results.moodpref_NegPos_hiprob))), sum(horzcat(results.moodpref_NegPos_loprob)==1) / sum(~isnan(horzcat(results.moodpref_NegPos_loprob))) ] )

disp('--- Chose Neg>Neu, Pos>Neu (-4 to +4):')
disp([nanmean(horzcat(results.allprobs_neg_over_neu)) ,nanmean(horzcat(results.allprobs_pos_over_neu))])

disp('--- % chose Neg>Neu, Pos>Neu:')
disp([nanstd(horzcat(results.allprobs_neg_over_neu)) ,nanstd(horzcat(results.allprobs_pos_over_neu))])

data=[nanmean(horzcat(results.allprobs_neg_over_neu)),nanmean(horzcat(results.allprobs_pos_over_neu))]
sd = [nanstd(horzcat(results.allprobs_neg_over_neu)) ,nanstd(horzcat(results.allprobs_pos_over_neu))]
se = sd / sqrt( length( horzcat(results.allprobs_neg_over_neu) ))

figure
bar(1:2,data)
xlabel('Mood contrast'); ylabel('Preference (-4 to 4)'); set(gca,'xticklabel',{'Neg>Neu'; 'Pos>Neu';}); %set(gca, 'YLim', [-4 4])
hold on
errorbar([data(1) data(2)],sd,'.')

% figure
% bar(1,-data(3))
% xlabel('Mood contrast'); ylabel('Preference (-4 to 4)'); set(gca,'xticklabel',{'Pos>Neg'}); set(gca, 'YLim', [-4 4])
% hold on
% errorbar(-data(3),sd(1),'.')





%%

% load 'allHPS.mat'
% 
% for i=1:length(results)
%     subjectid = str2num(strrep(files(i).name(19:20), '_', ''));
%     results(i).subjectid = subjectid;
%     results(i).HPS = xHPS(subjectid);
% end
%     
% hps=horzcat(results.HPS);
hps=nan(size(alldata,2),1);
ind_validhps=find(arrayfun(@(alldata) ~isempty(alldata.HPS),alldata));
hps(ind_validhps) = vertcat(alldata.HPS);
hpssplit = nanmedian(hps);
hihps = hps>hpssplit;

%hps_stretch(1:2:length(hps)*2) = hps; hps_stretch(2:2:length(hps)*2) = hps; 
%hihps_stretch(1:2:length(hihps)*2) = hihps; hps_stretch(2:2:length(hihps)*2) = hihps


%%
negneu = horzcat(results.allprobs_neg_over_neu);
posneu = horzcat(results.allprobs_pos_over_neu);
%posneg = horzcat(results.allprobs_pos_over_neg);

negneu(isnan(negneu)) = negneu(~isnan(negneu)); % Merge two test blocks to get rid of NaN 
posneu(isnan(posneu)) = posneu(~isnan(posneu)); % Merge two test blocks to get rid of NaN 
%posneg(isnan(posneg)) = posneg(~isnan(posneg)); % Merge two test blocks to get rid of NaN 
negneu=negneu(1:2:end);
posneu=posneu(1:2:end);
%posneg=posneg(1:2:end);

t=vertcat(results.testlearning_chosehiprob);
disp('--- %chose hi prob stim (t1; t2):')
disp('        Lo HPS                Hi HPS')
disp( [mean(t(~hihps,:)), mean(t(hihps,:))] )
[std(t(~hihps,:)), std(t(hihps,:))]

disp('--- Diff choice Neg>Neu, Pos>Neu, split HPS group (lo, hi):')
disp([nanmean(negneu(~hihps)), nanmean(negneu(hihps))])
disp([nanmean(posneu(~hihps)), nanmean(posneu(hihps))])
%disp([nanmean(posneg(~hihps)), nanmean(posneg(hihps))])
%disp('--- %pref for P>N, split HPS group (lo, hi):')
%conf_PN=vertcat(results.allprefs_PN_conf); conf_PN=conf_PN(:,2);
%disp( [nanmean(conf_PN(~hihps)), nanmean(conf_PN(hihps))] )


%%

moodpref_NeuNeg_hiprob = horzcat(results.moodpref_NeuNeg_hiprob);
%moodpref_NegPos_hiprob = horzcat(results.moodpref_NegPos_hiprob);
moodpref_NeuPos_hiprob = horzcat(results.moodpref_NeuPos_hiprob);

moodpref_NeuNeg_hiprob(isnan(moodpref_NeuNeg_hiprob)) = moodpref_NeuNeg_hiprob(~isnan(moodpref_NeuNeg_hiprob)); % Merge two test blocks to get rid of NaN 
%moodpref_NegPos_hiprob(isnan(moodpref_NegPos_hiprob)) = moodpref_NegPos_hiprob(~isnan(moodpref_NegPos_hiprob)); % Merge two test blocks to get rid of NaN 
moodpref_NeuPos_hiprob(isnan(moodpref_NeuPos_hiprob)) = moodpref_NeuPos_hiprob(~isnan(moodpref_NeuPos_hiprob)); % Merge two test blocks to get rid of NaN 
moodpref_NeuNeg_hiprob=moodpref_NeuNeg_hiprob(1:2:end);
%moodpref_NegPos_hiprob=moodpref_NegPos_hiprob(1:2:end);
moodpref_NeuPos_hiprob=moodpref_NeuPos_hiprob(1:2:end);

disp('--- % choice Neg>Neu, Pos>Neu for matched EV, HIPROB:')
disp(sum(moodpref_NeuNeg_hiprob==-1)/sum(~isnan(moodpref_NeuNeg_hiprob)))
disp(sum(moodpref_NeuPos_hiprob==1)/sum(~isnan(moodpref_NeuPos_hiprob)))
%disp(sum(moodpref_NegPos_hiprob==1)/sum(~isnan(moodpref_NegPos_hiprob)))

disp('above by HPS group:')
disp([sum(moodpref_NeuNeg_hiprob(~hihps)==-1)/sum(~hihps), sum(moodpref_NeuNeg_hiprob(hihps)==-1)/sum(hihps)])
disp([sum(moodpref_NeuPos_hiprob(~hihps)==1)/sum(~hihps), sum(moodpref_NeuPos_hiprob(hihps)==1)/sum(hihps)])
disp([sum(moodpref_NegPos_hiprob(~hihps)==1)/sum(~hihps),sum(moodpref_NegPos_hiprob(hihps)==1)/sum(hihps)])

moodpref_NeuNeg_loprob = horzcat(results.moodpref_NeuNeg_loprob);
moodpref_NegPos_loprob = horzcat(results.moodpref_NegPos_loprob);
moodpref_NeuPos_loprob = horzcat(results.moodpref_NeuPos_loprob);

moodpref_NeuNeg_loprob(isnan(moodpref_NeuNeg_loprob)) = moodpref_NeuNeg_loprob(~isnan(moodpref_NeuNeg_loprob)); % Merge two test blocks to get rid of NaN 
moodpref_NegPos_loprob(isnan(moodpref_NegPos_loprob)) = moodpref_NegPos_loprob(~isnan(moodpref_NegPos_loprob)); % Merge two test blocks to get rid of NaN 
moodpref_NeuPos_loprob(isnan(moodpref_NeuPos_loprob)) = moodpref_NeuPos_loprob(~isnan(moodpref_NeuPos_loprob)); % Merge two test blocks to get rid of NaN 
moodpref_NeuNeg_loprob=moodpref_NeuNeg_loprob(1:2:end);
moodpref_NegPos_loprob=moodpref_NegPos_loprob(1:2:end);
moodpref_NeuPos_loprob=moodpref_NeuPos_loprob(1:2:end);


disp('--- % choice Neg>Neu, Pos>Neu, Pos>Neg for matched EV, LOPROB:')
disp(sum(moodpref_NeuNeg_loprob==-1)/sum(~isnan(moodpref_NeuNeg_loprob)))
disp(sum(moodpref_NeuPos_loprob==1)/sum(~isnan(moodpref_NeuPos_loprob)))
disp(sum(moodpref_NegPos_loprob==1)/sum(~isnan(moodpref_NegPos_loprob)))

disp('above by HPS group (lo, hi):')
disp([sum(moodpref_NeuNeg_loprob(~hihps)==-1)/sum(~hihps),sum(moodpref_NeuNeg_loprob(hihps)==-1)/sum(hihps)])
disp([sum(moodpref_NeuPos_loprob(~hihps)==1)/sum(~hihps),sum(moodpref_NeuPos_loprob(hihps)==1)/sum(hihps)])
disp([sum(moodpref_NegPos_loprob(~hihps)==1)/sum(~hihps),sum(moodpref_NegPos_loprob(hihps)==1)/sum(hihps)])


%% Plot scatter
% figure;plot(hps, moodpref_NeuNeg_hiprob==-1,'.', 'MarkerSize',30)
% figure;plot(hps, moodpref_NeuPos_hiprob==1,'.', 'MarkerSize',30)
% figure;plot(hps, moodpref_NegPos_hiprob==1,'.', 'MarkerSize',30)

%% Plot HPS against mood manipulation, mood bias (binary and continuous)
x=hps(~isnan(hps))';

y=mooddiff_PoverN(~isnan(hps)); % P>N

%oldhps=vertcat(Olddata.HPS)
%x=[x, oldhps(~isnan(oldhps))'];
%y=[y, mooddiff_PoverNeuOLD(~isnan(oldhps))];
figure('color',[1 1 1]); subplot(2,3,1); 
plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Mood bias (Pos>Neg)');%set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('HPS-mooddiffP>N: r=%.02f, p=%.02f\n',r(2), p(2));

subplot(2,3,2);
y=vertcat(results.allprobs_pos_over_neg); y=y(:,2);
y=y(~isnan(hps))';

%yy=vertcat(Oldresults.allprobs_pos_over_neg); yy=yy(:,2);
%yy=yy(~isnan(oldhps)');
%y=[y;yy]';

%y=moodpref_NegPos_loprob' % HPS is correlated with mood bias in lo prob!
%y=moodpref_NegPos_hiprob' % HPS not correlated with mood bias in hi prob 

plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
xlabel('HPS'); ylabel('Preference (-4 to +4)');set(gca, 'YLim', [-4 4]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('HPS-choseP>N: r=%.02f, p=%.02f\n',r(2), p(2));

subplot(2,3,3);
y= vertcat(results.allprefs_PN_conf)*100;y=y(:,2);
y=y(~isnan(hps))';
plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('HPS-prefstrengthP>N: r=%.02f, p=%.02f\n',r(2), p(2));

subplot(2,3,4); title('HPS against mood bias (lo prob only)')
y=moodpref_NegPos_loprob'; % HPS is correlated with mood bias in lo prob!
%y=moodpref_NegPos_hiprob' % HPS not correlated with mood bias in hi prob 
y=y(~isnan(hps))';
plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
xlabel('HPS'); ylabel('Preference (-4 to +4)');set(gca, 'YLim', [-4 4]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
%[r, p]=corrcoef(x,y); fprintf('HPS-choseP>N: r=%.02f, p=%.02f\n',r(2), p(2));

subplot(2,3,5);title('HPS against mood bias (hi prob only)')
%y=moodpref_NegPos_loprob'; % HPS is correlated with mood bias in lo prob!
y=moodpref_NegPos_hiprob'; % HPS not correlated with mood bias in hi prob 
y=y(~isnan(hps))';
plot(x,y, '.','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
xlabel('HPS'); ylabel('Preference (-4 to +4)');set(gca, 'YLim', [-4 4]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
%[r, p]=corrcoef(x,y); fprintf('HPS-choseP>N: r=%.02f, p=%.02f\n',r(2), p(2));




%% Plot mood induction against mood bias

% P vs N contrast only possible in 2nd test block
x=mooddiff_PoverN'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
%x=mooddiff_PoverN_zscore'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
% Contrasts vs Neutral can be in either test block
%t2=t; t2(isnan(t2(:,1))) = t2(~isnan(t2(:,2)),2); % First merge two test blocks 
%x=mooddiff_NoverNeu; y=t2(:,1); % N>Neu
%x=mooddiff_PoverNeu; y=t2(:,1); % P>Neu
figure('color',[1 1 1]); subplot(2,3,1); plot(x,y, '.','MarkerSize',12); title('Mood induction against mood bias in binary ALL PROB');
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));


x=mooddiff_PoverN; y=moodpref_NegPos_hiprob;
%x=mooddiff_PoverN_zscore'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
% Contrasts vs Neutral can be in either test block
%t2=t; t2(isnan(t2(:,1))) = t2(~isnan(t2(:,2)),2); % First merge two test blocks 
%x=mooddiff_NoverNeu; y=t2(:,1); % N>Neu
%x=mooddiff_PoverNeu; y=t2(:,1); % P>Neu
subplot(2,3,2);plot(x,y, '.','MarkerSize',12); title('Mood induction against mood bias in binary HI PROB');
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));

x=mooddiff_PoverN; y=moodpref_NegPos_loprob;
%x=mooddiff_PoverN_zscore'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
% Contrasts vs Neutral can be in either test block
%t2=t; t2(isnan(t2(:,1))) = t2(~isnan(t2(:,2)),2); % First merge two test blocks 
%x=mooddiff_NoverNeu; y=t2(:,1); % N>Neu
%x=mooddiff_PoverNeu; y=t2(:,1); % P>Neu
subplot(2,3,3);plot(x,y, '.','MarkerSize',12); title('Mood induction against mood bias in binary LO PROB');
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));



y= vertcat(results.allprefs_PN_conf)*100;y=y(:,2);
%x= mooddiff_PoverN_zscore';
x= mooddiff_PoverN';
subplot(2,3,4); plot(x,y, '.','MarkerSize',12); title('Mood induction against mood bias in confidence ALL PROB');
xlabel('Mood induced (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias pref ratingP>N: r=%.02f, p=%.02f\n',r(2), p(2));


y= vertcat(results.allprefs_PN_conf_hiprob)*100;y=y(:,2);
%x= mooddiff_PoverN_zscore';
x= mooddiff_PoverN';
subplot(2,3,5); plot(x,y, '.','MarkerSize',12); title('Mood induction against mood bias in confidence HI PROB');
xlabel('Mood induced (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias pref ratingP>N: r=%.02f, p=%.02f\n',r(2), p(2));

y= vertcat(results.allprefs_PN_conf_loprob)*100;y=y(:,2);
%x= mooddiff_PoverN_zscore';
x= mooddiff_PoverN';
subplot(2,3,6); plot(x,y, '.','MarkerSize',12); title('Mood induction against mood bias in confidence LO PROB');
xlabel('Mood induced (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias pref ratingP>N: r=%.02f, p=%.02f\n',r(2), p(2));

%% Plot mood induction against mood effect separately for low and high HPS

y= vertcat(results.allprefs_PN_conf)*100;y=y(:,2);
%x= mooddiff_PoverN_zscore';
x= mooddiff_PoverN';

figure('color',[1 1 1]);plot(x(hihps),y(hihps), 'r.','MarkerSize',12)
hold on; plot(x(~hihps),y(~hihps), 'b.','MarkerSize',12)
xlabel('Mood induced (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);legend('hi hps', 'lo hps')
% Get fitted values
coeffs = polyfit(x(hihps), y(hihps), 1); fittedX = linspace(min(x(hihps)), max(x(hihps)), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
coeffs = polyfit(x(~hihps), y(~hihps), 1); fittedX = linspace(min(x(~hihps)), max(x(~hihps)), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);

%[r, p]=corrcoef(x,y)


%% Plot HPS against mood manipulation, mood bias (binary and continuous) by RUN

% RUN 1
tmp=vertcat(alldata.taskrun);
x=hps; x=x(tmp==1); t=find(isnan(x)); x=x(~isnan(x))';
y=mooddiff_PoverN; y=y(tmp==1); y(t)=[]; % P>N
figure('color',[1 1 1]); suptitle('HPS against mood manipulation & bias by RUN'); subplot(2,3,1); 
plot(x,y, '.r','MarkerSize',12); hold on;
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Mood bias (Pos>Neg)');%set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('HPS-mooddiffP>N: r=%.02f, p=%.02f\n',r(2), p(2));
%RUN 2
x=hps; x=x(tmp==2); t=find(isnan(x)); x=x(~isnan(x))';
y=mooddiff_PoverN; y=y(tmp==2); y(t)=[]; % P>N
hold on; subplot(2,3,1); plot(x,y, '.b','MarkerSize',12); hold on;
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Mood bias (Pos>Neg)');%set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('HPS-mooddiffP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'Run 1', 'Run 2');


%RUN 1
subplot(2,3,2);
x=hps; x=x(tmp==1); t=find(isnan(x)); x=x(~isnan(x))';
y=vertcat(results.allprobs_pos_over_neg); y=y(:,2); y=y(tmp==1)'; y(t)=[];
%y=moodpref_NegPos_loprob' % HPS is correlated with mood bias in lo prob!
%y=moodpref_NegPos_hiprob' % HPS not correlated with mood bias in hi prob 
plot(x,y, '.r','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
xlabel('HPS'); ylabel('Preference (-4 to +4)');set(gca, 'YLim', [-4 4]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('HPS-choseP>N: r=%.02f, p=%.02f\n',r(2), p(2));

%RUN 2
hold on; subplot(2,3,2);
x=hps; x=x(tmp==2); t=find(isnan(x)); x=x(~isnan(x))';
y=vertcat(results.allprobs_pos_over_neg); y=y(:,2); y=y(tmp==2)'; y(t)=[];
%y=moodpref_NegPos_loprob' % HPS is correlated with mood bias in lo prob!
%y=moodpref_NegPos_hiprob' % HPS not correlated with mood bias in hi prob 
plot(x,y, '.b','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
xlabel('HPS'); ylabel('Preference (-4 to +4)');set(gca, 'YLim', [-4 4]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('HPS-choseP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'Run 1', 'Run 2');


%RUN 1
subplot(2,3,3);
x=hps; x=x(tmp==1); t=find(isnan(x)); x=x(~isnan(x))';
y= vertcat(results.allprefs_PN_conf)*100;y=y(:,2); y=y(tmp==1)'; y(t)=[];
plot(x,y, '.r','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('HPS-prefstrengthP>N: r=%.02f, p=%.02f\n',r(2), p(2));

%RUN 2
subplot(2,3,3);
x=hps; x=x(tmp==2); t=find(isnan(x)); x=x(~isnan(x))';
y= vertcat(results.allprefs_PN_conf)*100;y=y(:,2); y=y(tmp==2)'; y(t)=[];
plot(x,y, '.b','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('HPS-prefstrengthP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'Run 1', 'Run 2');

%RUN 1
subplot(2,3,4); title('HPS against mood bias (lo prob only)')
x=hps; x=x(tmp==1); t=find(isnan(x)); x=x(~isnan(x))';
y=moodpref_NegPos_loprob'; % HPS is correlated with mood bias in lo prob!
%y=moodpref_NegPos_hiprob' % HPS not correlated with mood bias in hi prob 
y=y(tmp==1)'; y(t)=[]; 
plot(x,y, '.r','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
xlabel('HPS'); ylabel('Preference (-4 to +4)');set(gca, 'YLim', [-4 4]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
%[r, p]=corrcoef(x,y); fprintf('HPS-choseP>N: r=%.02f, p=%.02f\n',r(2), p(2));

%RUN 2
subplot(2,3,4); title('HPS against mood bias (lo prob only)')
x=hps; x=x(tmp==2); t=find(isnan(x)); x=x(~isnan(x))';
y=moodpref_NegPos_loprob'; % HPS is correlated with mood bias in lo prob!
%y=moodpref_NegPos_hiprob' % HPS not correlated with mood bias in hi prob 
y=y(tmp==2)'; y(t)=[]; 
plot(x,y, '.b','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
xlabel('HPS'); ylabel('Preference (-4 to +4)');set(gca, 'YLim', [-4 4]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
%[r, p]=corrcoef(x,y); fprintf('HPS-choseP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'Run 1', 'Run 2');


%RUN 1
subplot(2,3,5);title('HPS against mood bias (hi prob only)')
x=hps; x=x(tmp==1); t=find(isnan(x)); x=x(~isnan(x))';
%y=moodpref_NegPos_loprob'; % HPS is correlated with mood bias in lo prob!
y=moodpref_NegPos_hiprob'; % HPS not correlated with mood bias in hi prob 
y=y(tmp==1)'; y(t)=[]; 
plot(x,y, '.r','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
xlabel('HPS'); ylabel('Preference (-4 to +4)');set(gca, 'YLim', [-4 4]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
%[r, p]=corrcoef(x,y); fprintf('HPS-choseP>N: r=%.02f, p=%.02f\n',r(2), p(2));

%RUN 2
subplot(2,3,5);title('HPS against mood bias (hi prob only)')
x=hps; x=x(tmp==2); t=find(isnan(x)); x=x(~isnan(x))';
%y=moodpref_NegPos_loprob'; % HPS is correlated with mood bias in lo prob!
y=moodpref_NegPos_hiprob'; % HPS not correlated with mood bias in hi prob 
y=y(tmp==2)'; y(t)=[]; 
plot(x,y, '.b','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
xlabel('HPS'); ylabel('Preference (-4 to +4)');set(gca, 'YLim', [-4 4]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
%[r, p]=corrcoef(x,y); fprintf('HPS-choseP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'Run 1', 'Run 2');
%% Plot mood induction against mood bias Split by RUN

% P vs N contrast only possible in 2nd test block
% Separating plot of mood induction against bias by run
% Run 1
tmp=vertcat(alldata.taskrun);
x=mooddiff_PoverN(tmp==1)'; y=vertcat(results.allprobs_pos_over_neg); y=y(:,2); y=y(tmp==1); % P>N
figure('color',[1 1 1]); suptitle('Mood induction against mood bias by RUN'); subplot(2,3,1); plot(x,y, '.r','MarkerSize',12); 
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4]);
hold on
% Get fitted values: Run 1
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));


%Run 2
tmp=vertcat(alldata.taskrun);
x=mooddiff_PoverN(tmp==2)'; y=vertcat(results.allprobs_pos_over_neg); y=y(:,2); y=y(tmp==2); % P>N
hold on; plot(x,y, '.b','MarkerSize',12);
% Get fitted values: Run 2
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'Run 1', 'Run 2');


%RUN 1
x=mooddiff_PoverN; x=x(tmp==1); y=moodpref_NegPos_hiprob; y=y(tmp==1);
%x=mooddiff_PoverN_zscore'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
% Contrasts vs Neutral can be in either test block
%t2=t; t2(isnan(t2(:,1))) = t2(~isnan(t2(:,2)),2); % First merge two test blocks 
%x=mooddiff_NoverNeu; y=t2(:,1); % N>Neu
%x=mooddiff_PoverNeu; y=t2(:,1); % P>Neu
subplot(2,3,2);plot(x,y, '.r','MarkerSize',12); title('Mood induction against mood bias in binary HI PROB');
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));

%RUN 2
x=mooddiff_PoverN; x=x(tmp==2); y=moodpref_NegPos_hiprob; y=y(tmp==2);
%x=mooddiff_PoverN_zscore'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
% Contrasts vs Neutral can be in either test block
%t2=t; t2(isnan(t2(:,1))) = t2(~isnan(t2(:,2)),2); % First merge two test blocks 
%x=mooddiff_NoverNeu; y=t2(:,1); % N>Neu
%x=mooddiff_PoverNeu; y=t2(:,1); % P>Neu
subplot(2,3,2);plot(x,y, '.b','MarkerSize',12);
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'RUN 1', 'RUN 2', 'Location', 'NorthEast');



%RUN 1
x=mooddiff_PoverN; x=x(tmp==1); y=moodpref_NegPos_loprob; y=y(tmp==1);
%x=mooddiff_PoverN_zscore'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
% Contrasts vs Neutral can be in either test block
%t2=t; t2(isnan(t2(:,1))) = t2(~isnan(t2(:,2)),2); % First merge two test blocks 
%x=mooddiff_NoverNeu; y=t2(:,1); % N>Neu
%x=mooddiff_PoverNeu; y=t2(:,1); % P>Neu
subplot(2,3,3);plot(x,y, '.r','MarkerSize',12); title('Mood induction against mood bias in binary LO PROB');
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));

%RUN 2
x=mooddiff_PoverN; x=x(tmp==2); y=moodpref_NegPos_loprob; y=y(tmp==2);
%x=mooddiff_PoverN_zscore'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
% Contrasts vs Neutral can be in either test block
%t2=t; t2(isnan(t2(:,1))) = t2(~isnan(t2(:,2)),2); % First merge two test blocks 
%x=mooddiff_NoverNeu; y=t2(:,1); % N>Neu
%x=mooddiff_PoverNeu; y=t2(:,1); % P>Neu
subplot(2,3,3);plot(x,y, '.b','MarkerSize',12);
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'RUN 1', 'RUN 2', 'Location', 'NorthEast');



tmp=vertcat(alldata.taskrun); 
y= vertcat(results.allprefs_PN_conf)*100;y=y(:,2); y=y(tmp==1)';
%x= mooddiff_PoverN_zscore';
x= mooddiff_PoverN(tmp==1);
subplot(2,3,4); plot(x,y, '.r','MarkerSize',12); xlabel('Mood induced (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias pref ratingP>N: r=%.02f, p=%.02f\n',r(2), p(2));


tmp=vertcat(alldata.taskrun);
y= vertcat(results.allprefs_PN_conf)*100; y=y(:,2); y=y(tmp==2)';
x= mooddiff_PoverN(tmp==2);
plot(x,y, '.b','MarkerSize',12); xlabel('Mood induced (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias pref ratingP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'Run 1', 'Run 2');

%RUN 1
x=mooddiff_PoverN; x=x(tmp==1); y=vertcat(results.allprefs_PN_conf_hiprob)*100;y=y(:,2); y=y(tmp==1)';
subplot(2,3,5);plot(x,y, '.r','MarkerSize',12); title('Mood induction against mood bias in binary HI PROB');
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));

%RUN 2
x=mooddiff_PoverN; x=x(tmp==2); y=vertcat(results.allprefs_PN_conf_hiprob)*100;y=y(:,2); y=y(tmp==2)';
%x=mooddiff_PoverN_zscore'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
% Contrasts vs Neutral can be in either test block
%t2=t; t2(isnan(t2(:,1))) = t2(~isnan(t2(:,2)),2); % First merge two test blocks 
%x=mooddiff_NoverNeu; y=t2(:,1); % N>Neu
%x=mooddiff_PoverNeu; y=t2(:,1); % P>Neu
subplot(2,3,5);plot(x,y, '.b','MarkerSize',12); title('Mood induction against mood bias in confidence HI PROB');
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'RUN 1', 'RUN 2', 'Location', 'NorthEast');



%RUN 1
x=mooddiff_PoverN; x=x(tmp==1); y=vertcat(results.allprefs_PN_conf_loprob)*100; y=y(:,2); y=y(tmp==1)';
%x=mooddiff_PoverN_zscore'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
% Contrasts vs Neutral can be in either test block
%t2=t; t2(isnan(t2(:,1))) = t2(~isnan(t2(:,2)),2); % First merge two test blocks 
%x=mooddiff_NoverNeu; y=t2(:,1); % N>Neu
%x=mooddiff_PoverNeu; y=t2(:,1); % P>Neu
subplot(2,3,6); plot(x,y, '.r','MarkerSize',12); title('Mood induction against mood bias in confidence LO PROB');
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));

%RUN 2
x=mooddiff_PoverN; x=x(tmp==2); y=vertcat(results.allprefs_PN_conf_loprob)*100; y=y(:,2); y=y(tmp==2)';
%x=mooddiff_PoverN_zscore'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
% Contrasts vs Neutral can be in either test block
%t2=t; t2(isnan(t2(:,1))) = t2(~isnan(t2(:,2)),2); % First merge two test blocks 
%x=mooddiff_NoverNeu; y=t2(:,1); % N>Neu
%x=mooddiff_PoverNeu; y=t2(:,1); % P>Neu
subplot(2,3,6);plot(x,y, '.b','MarkerSize',12); title('Mood induction against mood bias in confidence LO PROB');
xlabel('Mood induction (Positive > Negative)');ylabel('Preference (%)');set(gca, 'YLim', [0 100]);
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'RUN 1', 'RUN 2', 'Location', 'NorthEast');


%% Plot mood induction against mood effect separately for low and high HPS by RUN

%RUN 1
tmp=vertcat(alldata.taskrun);
x= mooddiff_PoverN';x=x(tmp==1);t=find(isnan(x)); x=x(~isnan(x));
y= vertcat(results.allprefs_PN_conf)*100;y=y(:,2);y=y(tmp==1); y(t)=[]; % P>N
%x= mooddiff_PoverN_zscore';
z=hihps; z=z(tmp==1);
figure('color', [1 1 1]); subplot(1,2,1); plot(x(z),y(z), 'r.','MarkerSize',12); title('Mood induction against mood effect RUN 1');
hold on; plot(x(~z),y(~z), 'b.','MarkerSize',12)
xlabel('Mood induced (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);legend('hi hps', 'lo hps')
% Get fitted values
coeffs = polyfit(x(z), y(z), 1); fittedX = linspace(min(x(z)), max(x(z)), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
coeffs = polyfit(x(~z), y(~z), 1); fittedX = linspace(min(x(~z)), max(x(~z)), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
%[r, p]=corrcoef(x,y)


% RUN 2
x= mooddiff_PoverN';x=x(tmp==2);t=find(isnan(x)); x=x(~isnan(x));
y= vertcat(results.allprefs_PN_conf)*100;y=y(:,2);y=y(tmp==2); y(t)=[]; % P>N
z=hihps; z=z(tmp==2);
%x= mooddiff_PoverN_zscore';
subplot(1,2,2); plot(x(z),y(z), 'r.','MarkerSize',12); title('Mood induction against mood effect RUN 2');
hold on; plot(x(~z),y(~z), 'b.','MarkerSize',12)
xlabel('Mood induced (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);legend('hi hps', 'lo hps')
% Get fitted values
coeffs = polyfit(x(z), y(z), 1); fittedX = linspace(min(x(z)), max(x(z)), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
coeffs = polyfit(x(~z), y(~z), 1); fittedX = linspace(min(x(~z)), max(x(~z)), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);

%[r, p]=corrcoef(x,y)

%% Plot HPS against mood manipulation, mood bias (binary and continuous) by TASK VERSION

% VERSION 20 60
tmp=vertcat(alldata.taskver);
x=hps; x=x(tmp==1); t=find(isnan(x)); x=x(~isnan(x))';
y=mooddiff_PoverN; y=y(tmp==1); y(t)=[]; % P>N
figure('color',[1 1 1]); subplot(2,3,1); 
plot(x,y, '.r','MarkerSize',12); hold on;
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Mood bias (Pos>Neg)');%set(gca, 'YLim', [0 100]);
set(gca, 'XLim', [20 45])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('HPS-mooddiffP>N: r=%.02f, p=%.02f\n',r(2), p(2));
%VERSION 25 75
x=hps; x=x(tmp==2); t=find(isnan(x)); x=x(~isnan(x))';
y=mooddiff_PoverN; y=y(tmp==2); y(t)=[]; % P>N
hold on; subplot(2,3,1); plot(x,y, '.b','MarkerSize',12); hold on;
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Mood bias (Pos>Neg)');%set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('HPS-mooddiffP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, taskverstring{1}, taskverstring{2}, 'Location', 'SouthWest');


%VERSION 20 60
subplot(2,3,2);
x=hps; x=x(tmp==1); t=find(isnan(x)); x=x(~isnan(x))';
y=vertcat(results.allprobs_pos_over_neg); y=y(:,2); y=y(tmp==1)'; y(t)=[];
%y=moodpref_NegPos_loprob' % HPS is correlated with mood bias in lo prob!
%y=moodpref_NegPos_hiprob' % HPS not correlated with mood bias in hi prob 
plot(x,y, '.r','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
xlabel('HPS'); ylabel('Preference (-4 to +4)');set(gca, 'YLim', [-4 4]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('HPS-choseP>N: r=%.02f, p=%.02f\n',r(2), p(2));

%VERSION 25 75
hold on; subplot(2,3,2);
x=hps; x=x(tmp==2); t=find(isnan(x)); x=x(~isnan(x))';
y=vertcat(results.allprobs_pos_over_neg); y=y(:,2); y=y(tmp==2)'; y(t)=[];
%y=moodpref_NegPos_loprob' % HPS is correlated with mood bias in lo prob!
%y=moodpref_NegPos_hiprob' % HPS not correlated with mood bias in hi prob 
plot(x,y, '.b','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
xlabel('HPS'); ylabel('Preference (-4 to +4)');set(gca, 'YLim', [-4 4]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('HPS-choseP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, taskverstring{1}, taskverstring{2}, 'Location', 'SouthWest');

%VERSION 20 60
subplot(2,3,3);
x=hps; x=x(tmp==1); t=find(isnan(x)); x=x(~isnan(x))';
y= vertcat(results.allprefs_PN_conf)*100;y=y(:,2); y=y(tmp==1)'; y(t)=[];
plot(x,y, '.r','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('HPS-prefstrengthP>N: r=%.02f, p=%.02f\n',r(2), p(2));

%VERSION 25 75
subplot(2,3,3);
x=hps; x=x(tmp==2); t=find(isnan(x)); x=x(~isnan(x))';
y= vertcat(results.allprefs_PN_conf)*100;y=y(:,2); y=y(tmp==2)'; y(t)=[];
plot(x,y, '.b','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);xlabel('HPS'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('HPS-prefstrengthP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, taskverstring{1}, taskverstring{2}, 'Location', 'SouthWest');

%VERSION 20 60
subplot(2,3,4); title('HPS against mood bias (lo prob only)')
x=hps; x=x(tmp==1); t=find(isnan(x)); x=x(~isnan(x))';
y=moodpref_NegPos_loprob'; % HPS is correlated with mood bias in lo prob!
%y=moodpref_NegPos_hiprob' % HPS not correlated with mood bias in hi prob 
y=y(tmp==1)'; y(t)=[]; 
plot(x,y, '.r','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
xlabel('HPS'); ylabel('Preference (-4 to +4)');set(gca, 'YLim', [-4 4]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
%[r, p]=corrcoef(x,y); fprintf('HPS-choseP>N: r=%.02f, p=%.02f\n',r(2), p(2));

%VERSION 25 75
subplot(2,3,4); title('HPS against mood bias (lo prob only)')
x=hps; x=x(tmp==2); t=find(isnan(x)); x=x(~isnan(x))';
y=moodpref_NegPos_loprob'; % HPS is correlated with mood bias in lo prob!
%y=moodpref_NegPos_hiprob' % HPS not correlated with mood bias in hi prob 
y=y(tmp==2)'; y(t)=[]; 
plot(x,y, '.b','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
xlabel('HPS'); ylabel('Preference (-4 to +4)');set(gca, 'YLim', [-4 4]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
%[r, p]=corrcoef(x,y); fprintf('HPS-choseP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, taskverstring{1}, taskverstring{2}, 'Location', 'SouthWest');


%VERSION 20 60
subplot(2,3,5);title('HPS against mood bias (hi prob only)')
x=hps; x=x(tmp==1); t=find(isnan(x)); x=x(~isnan(x))';
%y=moodpref_NegPos_loprob'; % HPS is correlated with mood bias in lo prob!
y=moodpref_NegPos_hiprob'; % HPS not correlated with mood bias in hi prob 
y=y(tmp==1)'; y(t)=[]; 
plot(x,y, '.r','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
xlabel('HPS'); ylabel('Preference (-4 to +4)');set(gca, 'YLim', [-4 4]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
%[r, p]=corrcoef(x,y); fprintf('HPS-choseP>N: r=%.02f, p=%.02f\n',r(2), p(2));

%VERSION 25 75
subplot(2,3,5);title('HPS against mood bias (hi prob only)')
x=hps; x=x(tmp==2); t=find(isnan(x)); x=x(~isnan(x))';
%y=moodpref_NegPos_loprob'; % HPS is correlated with mood bias in lo prob!
y=moodpref_NegPos_hiprob'; % HPS not correlated with mood bias in hi prob 
y=y(tmp==2)'; y(t)=[]; 
plot(x,y, '.b','MarkerSize',12); hold on
% Get fitted values
coeffs = polyfit(x, y, 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
xlabel('HPS'); ylabel('Preference (-4 to +4)');set(gca, 'YLim', [-4 4]);
%set(gca, 'XLim', [0 1])
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
%[r, p]=corrcoef(x,y); fprintf('HPS-choseP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, taskverstring{1}, taskverstring{2}, 'Location', 'SouthWest');
hold on; suptitle('HPS against mood manipulation & bias by TASK VERSION');

%% Plot mood induction against mood bias Split by TASK VERSION

% P vs N contrast only possible in 2nd test block
% Separating plot of mood induction against bias by run
% VERSION 20 60

tmp=vertcat(alldata.taskver);
% *****
% Use this line if data from two tasks have been concatenated rather than
% blended as above *****
%tmp = nan(size(alldata)); %tmp(1:end/2)=1;tmp(1+end/2:end)=2; 


x=mooddiff_PoverN(tmp==1)'; y=vertcat(results.allprobs_pos_over_neg); y=y(:,2); y=y(tmp==1); % P>N
figure('color',[1 1 1]); subplot(2,3,1); plot(x,y, '.r','MarkerSize',12); 
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4]);
hold on
% Get fitted values: Run 1
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));


% VERSION 25 75
%tmp=vertcat(alldata.taskver);
x=mooddiff_PoverN(tmp==2)'; y=vertcat(results.allprobs_pos_over_neg); y=y(:,2); y=y(tmp==2); % P>N
hold on; plot(x,y, '.b','MarkerSize',12);
% Get fitted values: Run 2
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, taskverstring{1}, taskverstring{2}, 'Location', 'SouthEast');


% VERSION 20 60
x=mooddiff_PoverN; x=x(tmp==1); y=moodpref_NegPos_hiprob; y=y(tmp==1);
subplot(2,3,2);plot(x,y, '.r','MarkerSize',12); title('Mood induction against mood bias in binary HI PROB');
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));

% VERSION 25 75
x=mooddiff_PoverN; x=x(tmp==2); y=moodpref_NegPos_hiprob; y=y(tmp==2);
subplot(2,3,2);plot(x,y, '.b','MarkerSize',12);
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, taskverstring{1}, taskverstring{2}, 'Location', 'SouthEast');



% VERSION 20 60
x=mooddiff_PoverN; x=x(tmp==1); y=moodpref_NegPos_loprob; y=y(tmp==1);
subplot(2,3,3);plot(x,y, '.r','MarkerSize',12); title('Mood induction against mood bias in binary LO PROB');
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));

% VERSION 25 75
x=mooddiff_PoverN; x=x(tmp==2); y=moodpref_NegPos_loprob; y=y(tmp==2);
subplot(2,3,3);plot(x,y, '.b','MarkerSize',12);
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, taskverstring{1}, taskverstring{2}, 'Location', 'SouthEast');


% ALL PROBS
% VERSION 20 60
%tmp=vertcat(alldata.taskver); 
y= vertcat(results.allprefs_PN_conf)*100;y=y(:,2); y=y(tmp==1)';
%x= mooddiff_PoverN_zscore';
x= mooddiff_PoverN(tmp==1);
subplot(2,3,4); plot(x,y, '.r','MarkerSize',12); xlabel('Mood induced (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias pref ratingP>N: r=%.02f, p=%.02f\n',r(2), p(2));

% VERSION 25 75
%tmp=vertcat(alldata.taskver);
y= vertcat(results.allprefs_PN_conf)*100; y=y(:,2); y=y(tmp==2)';
%x= mooddiff_PoverN_zscore';
x= mooddiff_PoverN(tmp==2);
plot(x,y, '.b','MarkerSize',12); xlabel('Mood induced (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100])
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias pref ratingP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, taskverstring{1}, taskverstring{2}, 'Location', 'SouthEast');

% HI PROB
% VERSION 20 60
x=mooddiff_PoverN; x=x(tmp==1); y=vertcat(results.allprefs_PN_conf_hiprob)*100;y=y(:,2); y=y(tmp==1)';
subplot(2,3,5);plot(x,y, '.r','MarkerSize',12); title('Mood induction against mood bias in binary HI PROB');
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));

% VERSION 25 75
x=mooddiff_PoverN; x=x(tmp==2); y=vertcat(results.allprefs_PN_conf_hiprob)*100;y=y(:,2); y=y(tmp==2)';
subplot(2,3,5);plot(x,y, '.b','MarkerSize',12); title('Mood induction against mood bias in confidence HI PROB');
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, taskverstring{1}, taskverstring{2}, 'Location', 'SouthEast');


% LO PROB
% VERSION 20 60
x=mooddiff_PoverN; x=x(tmp==1); y=vertcat(results.allprefs_PN_conf_loprob)*100; y=y(:,2); y=y(tmp==1)';
subplot(2,3,6); plot(x,y, '.r','MarkerSize',12); title('Mood induction against mood bias in confidence LO PROB');
xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));

% VERSION 25 75
x=mooddiff_PoverN; x=x(tmp==2); y=vertcat(results.allprefs_PN_conf_loprob)*100; y=y(:,2); y=y(tmp==2)';
subplot(2,3,6);plot(x,y, '.b','MarkerSize',12); title('Mood induction against mood bias in confidence LO PROB');
xlabel('Mood induction (Positive > Negative)');ylabel('Preference (%)');set(gca, 'YLim', [0 100]);
hold on
% Get fitted values
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
[r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, taskverstring{1}, taskverstring{2}, 'Location', 'SouthEast');
hold on; suptitle('Mood induction against mood bias by TASK VERSION'); 

%% Plot mood induction against mood effect separately for low and high HPS by TASK VERSION

%
tmp=vertcat(alldata.taskver);
x= mooddiff_PoverN';x=x(tmp==1);t=find(isnan(x)); x=x(~isnan(x));
y= vertcat(results.allprefs_PN_conf)*100;y=y(:,2);y=y(tmp==1); y(t)=[]; % P>N
%x= mooddiff_PoverN_zscore';
z=hihps; z=z(tmp==1);
figure('color', [1 1 1]); subplot(1,2,1); plot(x(z),y(z), 'r.','MarkerSize',12); title('Mood induction against mood effect VERSION 20 60');
hold on; plot(x(~z),y(~z), 'b.','MarkerSize',12)
xlabel('Mood induced (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);legend('hi hps', 'lo hps')
% Get fitted values
coeffs = polyfit(x(z), y(z), 1); fittedX = linspace(min(x(z)), max(x(z)), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
coeffs = polyfit(x(~z), y(~z), 1); fittedX = linspace(min(x(~z)), max(x(~z)), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
%[r, p]=corrcoef(x,y)


% RUN 2
x= mooddiff_PoverN';x=x(tmp==2);t=find(isnan(x)); x=x(~isnan(x));
y= vertcat(results.allprefs_PN_conf)*100;y=y(:,2);y=y(tmp==2); y(t)=[]; % P>N
z=hihps; z=z(tmp==2);
%x= mooddiff_PoverN_zscore';
subplot(1,2,2); plot(x(z),y(z), 'r.','MarkerSize',12); title('Mood induction against mood effect VERSION 25 75');
hold on; plot(x(~z),y(~z), 'b.','MarkerSize',12)
xlabel('Mood induced (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]);legend('hi hps', 'lo hps')
% Get fitted values
coeffs = polyfit(x(z), y(z), 1); fittedX = linspace(min(x(z)), max(x(z)), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
coeffs = polyfit(x(~z), y(~z), 1); fittedX = linspace(min(x(~z)), max(x(~z)), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);

%[r, p]=corrcoef(x,y)

%% HPS vs Avg Mood

x=hps; x=x(~isnan(x)); tmp=vertcat(alldata.taskver);
y=vertcat(alldata.Moodhappyavg); y=y(~isnan(x),:);

%3 plots: one for each mood happiness average
figure('color',[1 1 1]); subplot(1,3,1); 
plot(x, y(:,1), 'k.', 'MarkerSize', 10); hold on; 
xlabel('HPS'); ylabel('Average Happiness'); title('Neutral Mood block');
coeffs = polyfit(x, y(:,1), 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'k-', 'LineWidth', 2);


% Positive
subplot(1,3,2); plot(x, y(:,2), 'g.', 'MarkerSize', 10); hold on; 
xlabel('HPS'); ylabel('Average Happiness'); title('Positive Mood block');
coeffs = polyfit(x, y(:,2), 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'g-', 'LineWidth', 2);

% Negative
subplot(1,3,3); plot(x, y(:,3), 'r.', 'MarkerSize', 10); hold on; 
xlabel('HPS'); ylabel('Average Happiness'); title('Negative Mood block');
coeffs = polyfit(x, y(:,3), 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);suptitle('Mood happiness x HPS');

%% Mood happiness x HPS RUN
tmp=vertcat(alldata.taskrun);
x=hps; x=x(tmp==1); t=find(isnan(x)); x=x(~isnan(x));
y=vertcat(alldata.Moodhappyavg); y=y(tmp==1,:); y(t,:)=[]; 

%Neutral RUN 1
figure('color',[1 1 1]); subplot(1,3,1); 
plot(x, y(:,1), 'r.', 'MarkerSize', 10); hold on; 
xlabel('HPS'); ylabel('Average Happiness'); title('Neutral Mood block');
coeffs = polyfit(x, y(:,1), 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);

%Neutral RUN 2
x=hps; x=x(tmp==2); t=find(isnan(x)); x=x(~isnan(x));
y=vertcat(alldata.Moodhappyavg); y=y(tmp==2,:); y(t,:)=[]; 
subplot(1,3,1); hold on; plot(x, y(:,1), 'b.', 'MarkerSize', 10); hold on; 
coeffs = polyfit(x, y(:,1), 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'Run 1', 'Run 2', 'Location', 'SouthWest');


% Positive RUN 1
x=hps; x=x(tmp==1); t=find(isnan(x)); x=x(~isnan(x));
y=vertcat(alldata.Moodhappyavg); y=y(tmp==1,:); y(t,:)=[]; 
subplot(1,3,2); plot(x, y(:,2), 'r.', 'MarkerSize', 10); hold on; 
xlabel('HPS'); ylabel('Average Happiness'); title('Positive Mood block');
coeffs = polyfit(x, y(:,2), 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);

%Positive RUN 2
x=hps; x=x(tmp==2); t=find(isnan(x)); x=x(~isnan(x));
y=vertcat(alldata.Moodhappyavg); y=y(tmp==2,:); y(t,:)=[]; 
subplot(1,3,2); hold on; plot(x, y(:,1), 'b.', 'MarkerSize', 10); 
hold on; coeffs = polyfit(x, y(:,2), 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'Run 1', 'Run 2', 'Location', 'SouthWest');


% Negative RUN 1
x=hps; x=x(tmp==1); t=find(isnan(x)); x=x(~isnan(x));
y=vertcat(alldata.Moodhappyavg); y=y(tmp==1,:); y(t,:)=[]; 
subplot(1,3,3); plot(x, y(:,3), 'r.', 'MarkerSize', 10); hold on; 
xlabel('HPS'); ylabel('Average Happiness'); title('Negative Mood block');
coeffs = polyfit(x, y(:,3), 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
hold on; plot(fittedX, fittedY, 'r-', 'LineWidth', 2);

%Negative RUN 2
x=hps; x=x(tmp==2); t=find(isnan(x)); x=x(~isnan(x));
y=vertcat(alldata.Moodhappyavg); y=y(tmp==2,:); y(t,:)=[]; 
subplot(1,3,3); hold on; plot(x, y(:,3), 'b.', 'MarkerSize', 10); 
hold on; coeffs = polyfit(x, y(:,3), 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
hold on; plot(fittedX, fittedY, 'b-', 'LineWidth', 2); 
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'Run 1', 'Run 2', 'Location', 'SouthWest');
suptitle('Mood happiness x HPS');

%% Mood happiness x HPS TASK VERSION
tmp=vertcat(alldata.taskver);
x=hps; x=x(tmp==1); t=find(isnan(x)); x=x(~isnan(x));
y=vertcat(alldata.Moodhappyavg); y=y(tmp==1,:); y(t,:)=[]; 
%Neutral Version 20 60
figure('color',[1 1 1]); subplot(1,3,1); 
plot(x, y(:,1), 'r.', 'MarkerSize', 10); axis([15 50 0 1]); hold on; 
xlabel('HPS'); ylabel('Average Happiness'); title('Neutral Mood block');
coeffs = polyfit(x, y(:,1), 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);

%Neutral Version 25 75
x=hps; x=x(tmp==2); t=find(isnan(x)); x=x(~isnan(x));
y=vertcat(alldata.Moodhappyavg); y=y(tmp==2,:); y(t,:)=[]; 
subplot(1,3,1); hold on; plot(x, y(:,1), 'b.', 'MarkerSize', 10); hold on; 
coeffs = polyfit(x, y(:,1), 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
hold on; plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, taskverstring{1}, taskverstring{2}, 'Location', 'SouthEast');

% Positive Version 20 60
x=hps; x=x(tmp==1); t=find(isnan(x)); x=x(~isnan(x));
y=vertcat(alldata.Moodhappyavg); y=y(tmp==1,:); y(t,:)=[]; 
subplot(1,3,2); plot(x, y(:,2), 'r.', 'MarkerSize', 10); axis([15 50 0 1]); hold on; 
xlabel('HPS'); ylabel('Average Happiness'); title('Positive Mood block');
coeffs = polyfit(x, y(:,2), 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);

%Positive Version 25 75
x=hps; x=x(tmp==2); t=find(isnan(x)); x=x(~isnan(x));
y=vertcat(alldata.Moodhappyavg); y=y(tmp==2,:); y(t,:)=[]; 
subplot(1,3,2); hold on; plot(x, y(:,1), 'b.', 'MarkerSize', 10); 
hold on; coeffs = polyfit(x, y(:,2), 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
hold on; plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, taskverstring{1}, taskverstring{2}, 'Location', 'SouthEast');

% Negative Version 20 60
x=hps; x=x(tmp==1); t=find(isnan(x)); x=x(~isnan(x));
y=vertcat(alldata.Moodhappyavg); y=y(tmp==1,:); y(t,:)=[]; 
subplot(1,3,3); plot(x, y(:,3), 'r.', 'MarkerSize', 10); hold on; 
xlabel('HPS'); ylabel('Average Happiness'); title('Negative Mood block');
coeffs = polyfit(x, y(:,3), 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
hold on; plot(fittedX, fittedY, 'r-', 'LineWidth', 2);

%Negative Version 25 75
x=hps; x=x(tmp==2); t=find(isnan(x)); x=x(~isnan(x));
y=vertcat(alldata.Moodhappyavg); y=y(tmp==2,:); y(t,:)=[]; 
subplot(1,3,3); hold on; plot(x, y(:,3), 'b.', 'MarkerSize', 10); axis([15 50 0 1]);
hold on; coeffs = polyfit(x, y(:,3), 1); fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
hold on; plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
suptitle('Mood happiness x HPS by Task Version');
h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); hold on;  legend(h, taskverstring{1}, taskverstring{2}, 'Location', 'SouthEast');



%% 3D plot: Perf, Happiness, Contrast

% x= vertcat(alldata.Moodhappyavg); x=x(:,2);
% y= vertcat(alldata.Moodperfavg); y=y(:,2);
% z= vertcat(results.allprobs_pos_over_neg); z=z(:,2);
% 
% [r,p] = corr(x,y)
% [r,p] = corr(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(1), p(1));
% 
% [r,p] = partialcorri(x,y,z)
% 
% 
% tbl = table(x,y,z,'VariableNames',{'Happy_average','Performance_average','Pos_Neg',});
% 
% lm = fitlm(tbl,'Pos_Neg~Performance_average')
% Res = table2array(lm.Residuals); plotResiduals(lm);
% 
% lm2 = fitlm(tbl,'Happy_average~Performance_average')
% Res2 = table2array(lm2.Residuals); plotResiduals(lm2);
% 
% [r,p] = corr(Res, Res2)
% 
% % scatter(Res(:,1), Res2(:,1))
% 


% 
% > X = c(2,4,15,20)
% > Y = c(1,2,3,4)
% > Z = c(0,0,1,1)
% > mm1 = lm(X~Z)
% > res1 = mm1$residuals
% > mm2 = lm(Y~Z)
% > res2 = mm2$residuals
% > cor(res1,res2)

% %% Mood manipulation and mood bias OLD vs NEW data
% 
% % % % DATASET 1
% % P vs N contrast only possible in 2nd test block
% x=mooddiff_PoverN'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
% %x=mooddiff_PoverN_zscore'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
% % Contrasts vs Neutral can be in either test block
% %t2=t; t2(isnan(t2(:,1))) = t2(~isnan(t2(:,2)),2); % First merge two test blocks 
% %x=mooddiff_NoverNeu; y=t2(:,1); % N>Neu
% %x=mooddiff_PoverNeu; y=t2(:,1); % P>Neu
% figure('color',[1 1 1]); subplot(2,3,1); plot(x,y, 'r.','MarkerSize',12); title('Mood induction against mood bias in binary ALL');
% xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4]); hold on
% % Get fitted values
% coeffs = polyfit(x, y, 1);
% fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
% plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
% [r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
% 
% % % % DATASET 2
% x=mooddiff_PoverNOLD'; y=vertcat(Oldresults.allprobs_pos_over_neg);y=y(:,2); % P>N
% subplot(2,3,1); plot(x,y, 'b.','MarkerSize',12);set(gca, 'YLim', [-4 4])
% hold on
% % Get fitted values
% coeffs = polyfit(x, y, 1);
% fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
% plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
% [r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
% hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'New Data', 'Old Data', 'Location', 'SouthEast');
% 
% 
% % % % BY RUN
% 
% % RUN 1
% % DATASET 1
% tmp=vertcat(alldata.taskrun);
% x=mooddiff_PoverN(tmp==1)'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); y=y(tmp==1); % P>N
% subplot(2,3,2); plot(x,y, 'r.','MarkerSize',12); title('Mood induction against mood bias in binary RUN 1');
% xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4]); hold on
% % Get fitted values
% coeffs = polyfit(x, y, 1);
% fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
% plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
% [r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
% 
% % DATASET 2
% tmp=vertcat(Olddata.taskrun);
% x=mooddiff_PoverNOLD(tmp==1)'; y=vertcat(Oldresults.allprobs_pos_over_neg);y=y(:,2); y=y(tmp==1);% P>N
% subplot(2,3,2); plot(x,y, 'b.','MarkerSize',12);set(gca, 'YLim', [-4 4])
% hold on
% % Get fitted values
% coeffs = polyfit(x, y, 1);
% fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
% plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
% [r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
% hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'New Data', 'Old Data', 'Location', 'SouthEast');
% 
% 
% % % % BY RUN
% 
% % RUN 2
% % DATASET 1
% tmp=vertcat(alldata.taskrun);
% x=mooddiff_PoverN(tmp==2)'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); y=y(tmp==2); % P>N
% subplot(2,3,3); plot(x,y, 'r.','MarkerSize',12); title('Mood induction against mood bias in binary RUN 2');
% xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4]); hold on
% % Get fitted values
% coeffs = polyfit(x, y, 1);
% fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
% plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
% [r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
% 
% % DATASET 2
% tmp=vertcat(Olddata.taskrun);
% x=mooddiff_PoverNOLD(tmp==2)'; y=vertcat(Oldresults.allprobs_pos_over_neg);y=y(:,2); y=y(tmp==2);% P>N
% subplot(2,3,3); plot(x,y, 'b.','MarkerSize',12);set(gca, 'YLim', [-4 4])
% hold on
% % Get fitted values
% coeffs = polyfit(x, y, 1);
% fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
% plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
% [r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
% hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'New Data', 'Old Data', 'Location', 'SouthEast');
% 
% 
% % % % CONFIDENCE
% 
% % DATASET 1
% x=mooddiff_PoverN'; y=vertcat(results.allprefs_PN_conf)*100;y=y(:,2); % P>N
% subplot(2,3,4); plot(x,y, 'r.','MarkerSize',12); title('Mood induction against mood bias in confidence ALL');
% xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]); hold on
% % Get fitted values
% coeffs = polyfit(x, y, 1);
% fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
% plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
% [r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
% 
% % DATASET 2
% x=mooddiff_PoverNOLD'; y=vertcat(Oldresults.allprefs_PN_conf)*100;y=y(:,2); % P>N
% subplot(2,3,4); plot(x,y, 'b.','MarkerSize',12); set(gca, 'YLim', [0 100]); hold on
% % Get fitted values
% coeffs = polyfit(x, y, 1);
% fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
% plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
% [r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias pref rating P>N: r=%.02f, p=%.02f\n',r(2), p(2));
% hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'New Data', 'Old Data', 'Location', 'SouthEast');
% 
% 
% % % % BY RUN
% 
% % RUN 1
% % DATASET 1
% tmp=vertcat(alldata.taskrun);
% x=mooddiff_PoverN(tmp==1)'; y=vertcat(results.allprefs_PN_conf)*100;y=y(:,2); y=y(tmp==1);
% subplot(2,3,5); plot(x,y, 'r.','MarkerSize',12); title('Mood induction against mood bias in confidence RUN 1');
% xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]); hold on
% % Get fitted values
% coeffs = polyfit(x, y, 1);
% fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
% plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
% [r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias pref rating P>N: r=%.02f, p=%.02f\n',r(2), p(2));
% 
% 
% % DATASET 2
% tmp=vertcat(Olddata.taskrun);
% x=mooddiff_PoverNOLD(tmp==1)'; y=vertcat(Oldresults.allprefs_PN_conf)*100;y=y(:,2); y=y(tmp==1);% P>N
% subplot(2,3,5); plot(x,y, 'b.','MarkerSize',12); set(gca, 'YLim', [0 100]); hold on
% % Get fitted values
% coeffs = polyfit(x, y, 1);
% fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
% plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
% [r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias pref rating P>N: r=%.02f, p=%.02f\n',r(2), p(2));
% hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'New Data', 'Old Data', 'Location', 'SouthEast');
% 
% 
% % RUN 2
% % DATASET 1
% tmp=vertcat(alldata.taskrun);
% x=mooddiff_PoverN(tmp==2)'; y=vertcat(results.allprefs_PN_conf)*100;y=y(:,2); y=y(tmp==2);
% subplot(2,3,6); plot(x,y, 'r.','MarkerSize',12); title('Mood induction against mood bias in confidence RUN 2');
% xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100]); hold on
% % Get fitted values
% coeffs = polyfit(x, y, 1);
% fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
% plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
% [r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias pref rating P>N: r=%.02f, p=%.02f\n',r(2), p(2));
% 
% 
% % DATASET 2
% tmp=vertcat(Olddata.taskrun);
% x=mooddiff_PoverNOLD(tmp==2)'; y=vertcat(Oldresults.allprefs_PN_conf)*100;y=y(:,2); y=y(tmp==2);% P>N
% subplot(2,3,6); plot(x,y, 'b.','MarkerSize',12); set(gca, 'YLim', [0 100]); hold on
% % Get fitted values
% coeffs = polyfit(x, y, 1);
% fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
% plot(fittedX, fittedY, 'b-', 'LineWidth', 2);
% [r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias pref rating P>N: r=%.02f, p=%.02f\n',r(2), p(2));
% hold on; h=NaN(2,1); h(1) = plot(NaN,'.r'); h(2)=plot(NaN,'.b'); legend(h, 'New Data', 'Old Data', 'Location', 'SouthEast');
% 
% 



%% Hi and Lo probs for mood bias vs mood manipulation -- NOT COMPLETED YET

% x=mooddiff_PoverN; y=moodpref_NegPos_hiprob;
% %x=mooddiff_PoverN_zscore'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
% % Contrasts vs Neutral can be in either test block
% %t2=t; t2(isnan(t2(:,1))) = t2(~isnan(t2(:,2)),2); % First merge two test blocks 
% %x=mooddiff_NoverNeu; y=t2(:,1); % N>Neu
% %x=mooddiff_PoverNeu; y=t2(:,1); % P>Neu
% subplot(2,3,2);plot(x,y, '.','MarkerSize',12); title('Mood induction against mood bias in binary HI PROB');
% xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4])
% hold on
% % Get fitted values
% coeffs = polyfit(x, y, 1);
% fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
% plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
% [r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
% 
% x=mooddiff_PoverN; y=moodpref_NegPos_loprob;
% %x=mooddiff_PoverN_zscore'; y=vertcat(results.allprobs_pos_over_neg);y=y(:,2); % P>N
% % Contrasts vs Neutral can be in either test block
% %t2=t; t2(isnan(t2(:,1))) = t2(~isnan(t2(:,2)),2); % First merge two test blocks 
% %x=mooddiff_NoverNeu; y=t2(:,1); % N>Neu
% %x=mooddiff_PoverNeu; y=t2(:,1); % P>Neu
% subplot(2,3,3);plot(x,y, '.','MarkerSize',12); title('Mood induction against mood bias in binary LO PROB');
% xlabel('Mood induction (Positive > Negative)'); ylabel('Preference (-4 to 4)');set(gca, 'YLim', [-4 4])
% hold on
% % Get fitted values
% coeffs = polyfit(x, y, 1);
% fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
% plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
% [r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias choiceP>N: r=%.02f, p=%.02f\n',r(2), p(2));
% 
% 
% y= vertcat(results.allprefs_PN_conf_hiprob)*100;y=y(:,2);
% %x= mooddiff_PoverN_zscore';
% x= mooddiff_PoverN';
% subplot(2,3,5); plot(x,y, '.','MarkerSize',12); title('Mood induction against mood bias in confidence HI PROB');
% xlabel('Mood induced (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100])
% hold on
% % Get fitted values
% coeffs = polyfit(x, y, 1);
% fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
% plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
% [r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias pref ratingP>N: r=%.02f, p=%.02f\n',r(2), p(2));
% 
% y= vertcat(results.allprefs_PN_conf_loprob)*100;y=y(:,2);
% %x= mooddiff_PoverN_zscore';
% x= mooddiff_PoverN';
% subplot(2,3,6); plot(x,y, '.','MarkerSize',12); title('Mood induction against mood bias in confidence LO PROB');
% xlabel('Mood induced (Positive > Negative)'); ylabel('Preference (%)');set(gca, 'YLim', [0 100])
% hold on
% % Get fitted values
% coeffs = polyfit(x, y, 1);
% fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
% plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
% [r, p]=corrcoef(x,y); fprintf('Mood diffP>N - Mood bias pref ratingP>N: r=%.02f, p=%.02f\n',r(2), p(2));


%%



