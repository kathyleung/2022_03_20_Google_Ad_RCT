clear all;
close all;

data_file_name = ['Data/Google Ad data from CK 6 Jan.xlsx'];
AA = readtable(data_file_name);

CLK = 1;
IMP = 2;
CTR = 3;

campaign_groups = unique(AA.Campaign);
age_groups = unique(AA.Age);
gender_text = unique(AA.Gender);
day_to_min = 60*24;
startDate = '2021-07-16';

time_partition = [0 14 28 50];
ad_metric = zeros(numel(age_groups)*numel(gender_text), numel(time_partition)-1, 3,2);
for ii=1:size(AA,1)
    AA.time_visit(ii) = (datenum(AA.Day(ii))-datenum(startDate));
    camp_str = char(AA.Campaign(ii));
    if (contains(camp_str, 'Traditional'))
        campaign = 1;
    elseif (contains(camp_str, 'English'))
        campaign = 2;
    else
        campaign = 3;
    end
    
    age_str = char(AA.Age(ii));
    
    for age=1:numel(age_groups)
        if (strcmp(age_str,age_groups(age)))
            break;
        end
    end
    gender_str = char(AA.Gender(ii));
    for gender=1:numel(gender_text)
        if (strcmp(gender_str, gender_text{gender}))
            break;
        end
    end
    
    for jj=1:numel(time_partition)-1
        if (AA.time_visit(ii)>=time_partition(jj) && AA.time_visit(ii)<time_partition(jj+1))
            stratum = age+(gender-1)*numel(age_groups);
            ad_metric(stratum, jj, campaign, CLK) = ad_metric(stratum, jj, campaign, CLK)+AA.Clicks(ii);
            ad_metric(stratum, jj, campaign, IMP) = ad_metric(stratum, jj, campaign, IMP)+AA.Impr_(ii);
        end
    end
end

ad_metric(:,:,:, CTR) = ad_metric(:,:,:, CLK)./ad_metric(:,:,:, IMP);
 
for qq=1:2
    for campaign=1:3
        for  jj=1:numel(time_partition)-1
            ad_metric(:, jj, campaign,qq) = ad_metric(:, jj, campaign,qq)/sum(ad_metric(:, jj, campaign,qq));
        end
    end
end

for gender=1:3
    for age=1:numel(age_groups)
        stratum = age+(gender-1)*numel(age_groups);
        stratum_text{stratum} = [age_groups{age} ' ' gender_text{gender}];
    end
end

week_text = {'Weeks 1-2', 'Weeks 3-4', 'Weeks 5-6'};
campaign_text = {'Traditional Chinese', 'English', 'Simplified Chinese'};
xv = 1:numel(age_groups)*numel(gender_text);

figNum = 0;
for qq=[CLK IMP CTR]
    
    figNum = figNum+1;
    figure(figNum);
    clf;
    hold on;
    for campaign=1:numel(campaign_text)
        subplot(numel(campaign_text),1,campaign);
        bar(ad_metric(:,:,campaign,qq));
        if (qq~=CTR)
            ylim([0 0.3]);
        end
        if (campaign==3)
            xlabel('Age-gender group');
            xtickangle(45);
            set(gca, 'XTick', xv, 'XTickLabel', stratum_text);
        else
            set(gca, 'XTick', xv, 'XTickLabel', []);
        end
        ylabel('Proportion');
        legend(week_text, 'Location', 'Best');
        legend boxoff;
        title(campaign_text{campaign});
    end
    figNum = figNum+1;
    figure(figNum);
    
    clf;
    hold on;
    for wk=1:numel(week_text)
        subplot(numel(week_text),1,wk);
        bar(squeeze(ad_metric(:,wk,:,qq)));
        if (qq~=CTR)
            ylim([0 0.3]);
        end
        if (wk==numel(week_text))
            xlabel('Age-gender group');
            xtickangle(45);
            set(gca, 'XTick', xv, 'XTickLabel', stratum_text);
        else
            set(gca, 'XTick', xv, 'XTickLabel', []);
        end
        ylabel('Proportion');
        legend(campaign_text, 'Location', 'Best');
        legend boxoff
        title(week_text{wk});
    end
end