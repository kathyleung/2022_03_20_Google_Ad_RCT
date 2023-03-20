clear all;
close all;

TRAD_CHINESE = 1;
ENGLISH = 2;
SIMP_CHINESE = 3;
language_text = {'Traditional Chinese', 'English', 'Simplified Chinese'};

category_full_text = {'Exemption', 'Comparison', 'Family', 'Mortality', 'Government', 'Reopen', 'Lottery'};
category_short_text = {'E', 'C', 'F', 'M', 'G', 'R', 'L'};
num_categories = numel(category_full_text);
category_display_order = [5 2 1 3 7 4 6];
reference_cat = 5;   % the control arm

%%%% time 0 is 2021-07-16 00:00:01
time_partition = [0 15 30 31+30];  % Days since time 0; study period is partitioned into Jul, Aug, Sep+Oct
time_partition = cumsum(time_partition);

age_groups = [
    12 29
    30 54
    55 100
    10 100
    ];

num_age_groups = size(age_groups,1);
age_txt = [];
for ii=1:size(age_groups,1)
    age_txt{end+1} = [num2str(age_groups(ii,1)) '-'  num2str(age_groups(ii,2))];;
end

age_txt{end} = 'All';
if (numel(time_partition)==4)
    tp_txt = {'Jul', 'Aug', 'Sep-Oct'};
else
    tp_txt = {'Jul', 'Aug', 'Sep', 'Oct'};
end

load('Results/script_2_output');
ds_idx_vec = 2;
lang_vec = [TRAD_CHINESE SIMP_CHINESE ENGLISH];
for num_bookings_ceiling = [100]
    for margin_for_session_split = [0.5 1 2]*day_to_min  % split into sessions if inter-visit time is largeer than margin_for_session_split
        for margin_for_merging = [2 3 6]
            % for num_bookings_ceiling = [100]
            %     for margin_for_session_split = [1]*day_to_min  % split into sessions if inter-visit time is larger than margin_for_session_split
            %         for margin_for_merging = [3]
            clear processed_linelist
            for lang=1:numel(linelist_large)
                clear linelist_m data_for_inference
                ii_m = 1;
                
                tv_diff_all = [];
                ip_lang_diff = [];
                ip_clash = [];
                linelist = linelist_large{lang};
                for ii=1:numel(linelist)
                    ii
                    xx = linelist(ii);
                    session_idx = 1;
                    
                    xx.time_visit_diff = 0;
                    xx.category_diff = 0;
                    xx.language_diff = 0;
                    xx.time_visit_split{session_idx} = xx.time_visit(1);
                    xx.category_split{session_idx} = xx.category(1);
                    xx.language_split{session_idx} = xx.language(1);
                    xx.click_button_split{session_idx} = xx.click_button(1);
                    
                    if (numel(xx.time_visit)>1)
                        xx.time_visit_diff = xx.time_visit(2:end)-xx.time_visit(1:end-1);
                        xx.category_diff = xx.category(2:end)-xx.category(1:end-1);
                        xx.language_diff = xx.language(2:end)-xx.language(1:end-1);
                        
                        for iii=1:2
                            if (iii==1)
                                fI = find(xx.time_visit_diff<margin_for_merging & xx.category_diff==0 & xx.language_diff==0);
                                fI = fI+1;                                
                            else
                                fI = find(xx.time_visit_diff<margin_for_session_split);
                                if (~isempty(fI))
                                    fI = unique([fI fI+1]);
                                end
                            end
                            
                            xx.time_visit(fI) = [];
                            xx.category(fI) = [];
                            xx.language(fI) = [];
                            xx.click_button(fI) = [];
                            xx.time_visit_diff = xx.time_visit(2:end)-xx.time_visit(1:end-1);
                            xx.category_diff = xx.category(2:end)-xx.category(1:end-1);
                            xx.language_diff = xx.language(2:end)-xx.language(1:end-1);
                        end
                        
                        for jj=2:numel(xx.time_visit)
                            if (xx.time_visit_diff(jj-1)>margin_for_session_split)
                                session_idx = session_idx+1;
                                next_idx = 1;
                            else
                                next_idx = next_idx+1;
                            end
                            xx.time_visit_split{session_idx}(next_idx) = xx.time_visit(jj);
                            xx.category_split{session_idx}(next_idx) = xx.category(jj);
                            xx.language_split{session_idx}(next_idx) = xx.language(jj);
                            xx.click_button_split{session_idx}(next_idx) = xx.click_button(jj);
                        end
                    end
                    if (~isempty(xx.time_visit))
                        linelist_m(ii_m) = xx;
                        ii_m = ii_m+1;
                        tv_diff_all = [tv_diff_all xx.time_visit_diff];
                    end
                end
                
                for ds_idx=ds_idx_vec
                    if (numel(linelist_large)==1)
                        episode_idx = [1 1 1];
                    else
                        episode_idx = 1;
                    end
                    for ii=1:numel(linelist_m)
                        ii
                        xx = linelist_m(ii);
                        if (ds_idx==1)
                            nn = 1;   % consider only the 1st session from each IP
                        else
                            nn = numel(xx.time_visit_split);    % consider all split sessions from each IP
                        end
                        
                        for session_idx=1:nn
                            if (numel(xx.time_visit_split{session_idx})==1)   %%%% choose only split sessions with one visit
                                lang_this = xx.language_split{session_idx};
                                if (numel(episode_idx)==1)
                                    episode_idx_this = episode_idx;
                                else
                                    episode_idx_this = episode_idx(lang_this);
                                end
                                processed_linelist{lang_this, ds_idx}.time_visit(episode_idx_this) = xx.time_visit_split{session_idx};
                                processed_linelist{lang_this, ds_idx}.category(episode_idx_this) = xx.category_split{session_idx};
                                processed_linelist{lang_this, ds_idx}.language(episode_idx_this) = xx.language_split{session_idx};
                                %                     processed_linelist{lang, ds_idx}.click_button(episode_idx) = xx.click_button_split{session_idx};
                                
                                bk_time = xx.booking_time-xx.time_visit_split{session_idx};
                                
                                fI = find(abs(bk_time)<margin_for_session_split);
                                mm = numel(fI);
                                processed_linelist{lang_this, ds_idx}.booking_time(1,episode_idx_this) = NaN;
                                processed_linelist{lang_this, ds_idx}.booking_vax(1,episode_idx_this) = -1;
                                processed_linelist{lang_this, ds_idx}.booking_age(1,episode_idx_this) = -1;
                                processed_linelist{lang_this, ds_idx}.vax_time(1,:,episode_idx_this) = [-1 -1];
                                if (mm>0)
                                    bk_time = bk_time(fI);
                                    vax = xx.booking_vax(fI);
                                    age = xx.booking_age(fI);
                                    vax_time = xx.vax_time(:,fI);
                                    
                                    processed_linelist{lang_this, ds_idx}.booking_time(1:mm,episode_idx_this) = bk_time';
                                    processed_linelist{lang_this, ds_idx}.booking_vax(1:mm,episode_idx_this) = vax';
                                    processed_linelist{lang_this, ds_idx}.booking_age(1:mm,episode_idx_this) = age';
                                    processed_linelist{lang_this, ds_idx}.vax_time(1:mm,:,episode_idx_this) = vax_time';
                                end
                                if (numel(episode_idx)==1)
                                    episode_idx = episode_idx+1;
                                else
                                    episode_idx(lang_this) = episode_idx(lang_this)+1;
                                end
                            else
                                disp('error');
                                return;
                                
                            end
                        end
                    end
                    if (numel(linelist_large)==1)
                        lang_vec_2 = 1:3;
                    else
                        lang_vec_2 = lang_this;
                    end
                    for lang_this_2=lang_vec_2
                        processed_linelist{lang_this_2, ds_idx}.booking_time = processed_linelist{lang_this_2, ds_idx}.booking_time';
                        processed_linelist{lang_this_2, ds_idx}.booking_time(processed_linelist{lang_this_2, ds_idx}.booking_time==0) = NaN;
                        processed_linelist{lang_this_2, ds_idx}.booking_vax = processed_linelist{lang_this_2, ds_idx}.booking_vax';
                        processed_linelist{lang_this_2, ds_idx}.booking_age = processed_linelist{lang_this_2, ds_idx}.booking_age';
                    end
                end
            end
            clear linelist_m
            
            %             for time_window_this = min(margin_for_session_split, [30 60 120 180 360 day_to_min])
                                    for time_window_this = [180 360 720]
%             for time_window_this = 360
                close all;
                clear  data_for_inference
                set(0, 'defaulttextfontname', 'arial');
                set(0, 'defaultaxesfontname', 'arial');
                set(0, 'defaulttextfontsize', 14)
                set(0, 'defaultaxesfontsize', 14);
                set(0, 'defaultlinelinewidth', 1);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                load('Results/script_1_output');
                
                timeMax = -1;
                timeMin = 1e10;
                age_vec = 2021-booking_data.BirthYear;
                for age_idx=1:num_age_groups
                    for vax=[SINOVAC BNT]
                        XX = booking_data.BookingCreationDateTime(fI{vax});
                        fI_age = age_vec(fI{vax})>=age_groups(age_idx,1) & age_vec(fI{vax})<=age_groups(age_idx,2);
                        XX = XX(fI_age & XX>=startDate_2 & XX<=endDate_2);
                        XX = (datenum(XX)-datenum(startDate))*day_to_min;
                        booking_time{age_idx,vax} = XX;
                        if (~isempty(XX))
                            timeMax = max(timeMax, max(XX));
                            timeMin = min(timeMin, min(XX));
                        end
                    end
                    booking_time{age_idx,ALL} = [booking_time{age_idx,SINOVAC}; booking_time{age_idx,BNT}];
                    booking_time{age_idx,ALL} = sort(booking_time{age_idx,ALL});
                end
                
                dt = 10;
                max_time_window = day_to_min;
                time_window_vec = dt:dt:max_time_window;
                time_window_indices(time_window_vec) = 1:numel(time_window_vec);
                tv = floor(timeMin):ceil(timeMax);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                tic
                yy_max = -1;
                for ds_idx=ds_idx_vec
                    for lang = lang_vec
                        
                        close all;
                        ds = processed_linelist{lang, ds_idx};
                        fI = find(ds.time_visit>0);
                        this_time_visit = ds.time_visit(fI);
                        
                        time_visit = ds.time_visit/day_to_min;
                        for cat_idx = 1:num_categories
                            for t_idx = 1:(numel(time_partition)-1)
                                
                                %%%%%%%%% 1st step: Select based on the characeteristics of the visit
                                %%%%%%%%% fI is the set of visits with the
                                %%%%%%%%% specified (cat, lang, click, time_period)
                                %%%%%%%%%
                                
                                fI_cl = find(ds.category==cat_idx);
                                if (t_idx==numel(time_partition))
                                    fI_t = find(time_visit>=time_partition(1) & time_visit<time_partition(end));
                                else
                                    fI_t = find(time_visit>=time_partition(t_idx) & time_visit<time_partition(t_idx+1));
                                end
                                fI_clt = intersect(fI_cl,fI_t);
                                
                                %%%%%%%%% 2nd step: Stratify
                                %%%%%%%%% bookings for each visit
                                %%%%%%%%% stratified by booking age and vax
                                if (isempty(ds.booking_time))
                                    bk_time_age_vax = [];
                                else
                                    XXX = ds.booking_time(fI_clt,:);
                                end
                                
                                fI_xx = find(sum(XXX<0, 2)==0 & sum(XXX>0, 2) < num_bookings_ceiling);
                                fI = fI_clt(fI_xx);
                                bk_time_age_vax_time = ds.booking_time(fI,:);
                                for age_idx = 1:num_age_groups
                                    for vax = ALL
                                        if (age_groups(age_idx,2)<18 && vax==SINOVAC)
                                            continue;
                                        end
                                        if (isempty(bk_time_age_vax_time))
                                            bk_time_age_vax = [];
                                        else
                                            ZZZ = ds.booking_age(fI,:);
                                            fI_a = ZZZ>=age_groups(age_idx,1) & ZZZ<=age_groups(age_idx,2);
                                            bk_time_age_vax = bk_time_age_vax_time.*fI_a;
                                            if (vax~=ALL)
                                                fI_v = ds.booking_vax(fI,:)==vax;
                                                bk_time_age_vax = bk_time_age_vax.*fI_v;
                                            end
                                            bk_time_age_vax(bk_time_age_vax==0) = NaN;
                                        end
                                    end
                                    time_visit_dist{cat_idx, t_idx, vax, age_idx, ds_idx, lang} = sort(time_visit(fI));
                                    
                                    %%%%%%%%% Tally the pre and post no. of age-vax booking
                                    time_window_idx = time_window_indices(time_window_this);
                                    
                                    fI_this = 1:numel(fI);
                                    mm = numel(fI);
                                    
                                    bk_time_age_vax_this = bk_time_age_vax(fI_this,:);
                                    
                                    yy = sum(bk_time_age_vax_this>0 & bk_time_age_vax_this<=time_window_this, 2);
                                    yy_max = max(max(yy(:)),yy_max);
                                    
                                    nn = numel(fI_this);
                                    
                                    data_for_inference.N(cat_idx, t_idx, vax, age_idx, ds_idx, time_window_idx, lang) = nn;
                                    data_for_inference.Y(cat_idx, t_idx, vax, age_idx, ds_idx, time_window_idx, lang) = sum(yy);
                                    data_for_inference.Y_mean(cat_idx, t_idx, vax, age_idx, ds_idx, time_window_idx, lang) = mean(yy);
                                    data_for_inference.Y_var(cat_idx, t_idx, vax, age_idx, ds_idx, time_window_idx, lang) = var(yy);
                                    
                                end
                            end
                        end
                    end
                end
                for ds_idx=ds_idx_vec
                    close all
                    suffix = ['_lang' num2str(numel(linelist_large))...
                        '_sp' num2str(margin_for_session_split) ...
                        '_merge' num2str(margin_for_merging) ...
                        '_tw' num2str(time_window_this) ...
                        '_bkceil' num2str(num_bookings_ceiling) ...
                        '_ds' num2str(ds_idx) ...
                        '_tp' num2str(numel(time_partition))]
                    
                    fig_num = 0;
                    
                    mu_mat = data_for_inference.Y_mean;
                    sigma_2_mat = data_for_inference.Y_var./data_for_inference.N;
                    
                    ntp = numel(time_partition);
                    for cat_idx=1:size(mu_mat,1)
                        for vax=1:size(mu_mat,3)
                            for age_idx=1:size(mu_mat,4)
                                for ds_idx=1:size(mu_mat,5)
                                    for time_window_idx=1:size(mu_mat,6)
                                        for lang=1:size(mu_mat,7)
                                            mu_mat(cat_idx, ntp, vax, age_idx, ds_idx, time_window_idx, lang) = ...
                                                mean(mu_mat(cat_idx, 1:ntp-1, vax, age_idx, ds_idx, time_window_idx, lang));
                                            sigma_2_mat(cat_idx, ntp, vax, age_idx, ds_idx, time_window_idx, lang) = ...
                                                mean(sigma_2_mat(cat_idx, 1:ntp-1, vax, age_idx, ds_idx, time_window_idx, lang));
                                            data_for_inference.N(cat_idx, ntp, vax, age_idx, ds_idx, time_window_idx, lang) = ...
                                                sum(data_for_inference.N(cat_idx, 1:ntp-1, vax, age_idx, ds_idx, time_window_idx, lang));
                                        end
                                    end
                                end
                            end
                        end
                    end
                    sigma_mat = sqrt(sigma_2_mat);
                    
                    
                    fig_num = fig_num+1;
                    fh = figure(fig_num);
                    fh.WindowState = 'maximized';
                    th = tiledlayout(numel(category_display_order)-1, 3,'TileSpacing','Compact','Padding','Compact');
                    zz =  get(th, 'Position');
                    zz(2) = zz(2)+0.06;
                    zz(4) = zz(4)-0.1;
                    set(th, 'Position', zz)
                    annotation(gcf,'textbox',...
                        [0.102083333333334 0.0164587778855482 0.860416666666666 0.0397672162948595],...
                        'String','Jul                      Aug                  Sep-Oct                             Jul                      Aug                  Sep-Oct                               Jul                      Aug                  Sep-Oct',...
                        'LineStyle','none',...
                        'FontSize',11,...
                        'FontName','Arial',...
                        'FontWeight','bold',...
                        'FitBoxToText','off');
                    
                    for vax =ALL
                        for ii = 1:numel(category_display_order(2:end))
                            for lang = lang_vec
                                cat_idx = category_display_order(1+ii);
                                nexttile
                                hold on;
                                
                                age_idx_vec = 1:num_age_groups;
                                t_idx_vec = 1:numel(time_partition)-1;
                                
                                for age_idx = age_idx_vec
                                    for t_idx=t_idx_vec
                                        mu_cat = mu_mat(cat_idx, t_idx, vax, age_idx, ds_idx, time_window_idx, lang);
                                        sigma_cat = sigma_mat(cat_idx, t_idx, vax, age_idx, ds_idx, time_window_idx, lang);
                                        mu_ref = mu_mat(reference_cat, t_idx, vax, age_idx, ds_idx, time_window_idx, lang);
                                        sigma_ref = sigma_mat(reference_cat, t_idx, vax, age_idx, ds_idx, time_window_idx, lang);
                                        
                                        eps_x = 1e-10;
                                        ref_support = norminv([eps_x 1-eps_x], mu_ref, sigma_ref);
                                        xx = linspace(max(1e-10,ref_support(1)), ref_support(end), 2000);
                                        dx = xx(2)-xx(1);
                                        if (normcdf(xx(1), mu_ref, sigma_ref)>0.01)
                                            disp(['Warning: Support contains negative numbers: ' num2str(normcdf(xx(1), mu_ref, sigma_ref))]);
                                            disp([language_text{lang} ' ' category_full_text{cat_idx}])
                                        end
                                        ref_pdf = normpdf(xx, mu_ref, sigma_ref)/(1-normcdf(xx(1), mu_ref, sigma_ref));
                                        
                                        dzz = 0.05;
                                        zz = -3;
                                        m_diff_cdf = 0;
                                        jj = 1;
                                        dx_0 = dx;
                                        xx_0 = xx;
                                        ref_pdf_0 = ref_pdf;
                                        while (m_diff_cdf(end)<0.996)
                                            cat_cdf = normcdf((1+zz(end))*xx, mu_cat, sigma_cat);
                                            m_diff_cdf(jj) = ref_pdf*cat_cdf'*dx;
                                            zz(jj+1) = zz(jj)+dzz;
                                            jj = jj+1;
                                        end
                                        pv = [0.5 2.5 5 50 95 97.5 99.5]/100;
                                        for jj=1:numel(pv)
                                            fI = find(m_diff_cdf>pv(jj), 1, 'first');
                                            if (fI==1)
                                                m_diff_CI(jj) = zz(fI);
                                            else
                                                m_diff_CI(jj) = interp1(m_diff_cdf([fI-1 fI]), zz([fI-1 fI]), pv(jj));
                                            end
                                        end
                                        
                                        if (all(m_diff_CI>0) || all(m_diff_CI<0))
                                            clr = [1 0 0];
                                        elseif (all(m_diff_CI([2 end-1])>0) || all(m_diff_CI([2 end-1])<0))
                                            clr = [1 165/255 0];
                                            %                                         elseif (all(m_diff_CI([3 end-2])>0) || all(m_diff_CI([3 end-2])<0))
                                            %                                             clr = [1 1 0];
                                        else
                                            clr = [0.5 0.5 0.5];
                                        end
                                        
                                        marker_style = 'd^vo';
                                        ms = marker_style(age_idx);
                                        dxx = linspace(-0.3, 0.3, num_age_groups);
                                        
                                        xv = [];
                                        xvl = [];
                                        for iii=1:num_age_groups
                                            xv = [xv t_idx_vec+dxx(iii)];
                                            xvl = [xvl {'12-29', '30-54', '>54', 'Overall'}];
                                        end
                                        xv = sort(xv);
                                        
                                        xx = t_idx+dxx(age_idx);
                                        
                                        if (all(m_diff_CI([2 end-1])>0) || all(m_diff_CI([2 end-1])<0))
                                            disp([language_text{lang} ' ' category_full_text{cat_idx} ' ' xvl{t_idx} ' ' age_txt{age_idx} ': ' num2str(m_diff_CI([2 4 end-1]))])
                                        end
                                        
                                        if (age_idx==1 && t_idx==1)
                                            plot([0 xv(end)+1], [0 0], 'k-');
                                        end
                                        plot([xx xx], m_diff_CI([2 end-1]), 'k-');
                                        plot(xx, m_diff_CI((numel(m_diff_CI)+1)/2), ms, 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', clr);
                                        xlim([xv(1)-0.25 xv(end)+0.25]);
                                        
                                        if (cat_idx==category_display_order(end))
                                            set(gca, 'XTick', xv, 'XTickLabel', xvl);
                                        else
                                            set(gca, 'XTick', xv, 'XTickLabel', []);
                                        end
                                        xtickangle(45);
                                        if (cat_idx==category_display_order(2))
                                            num_subjects = sum(data_for_inference.N(1:num_categories, end, ALL, end, ds_idx, time_window_idx, lang), 'all');
                                            title({language_text{lang}; [' (n = ' num2str(num_subjects) ')']});
                                        end
                                        if (lang==TRAD_CHINESE)
                                            ylabel(category_full_text{cat_idx});
                                        end
                                        
                                        if (lang==TRAD_CHINESE)
                                            yv = -0.5:0.5:1;
                                        elseif (lang==SIMP_CHINESE)
                                            yv = -1:1:2;
                                        else
                                            yv = -1:1:2;
                                        end
                                        set(gca, 'YTick', yv, 'YTickLabel', yv);
                                        ylim(yv([1 end]));
                                    end
                                end
                            end
                        end
                    end
                    saveas(gcf, ['Figures/fig2' suffix]);
                    saveas(gcf, ['Figures/fig2' suffix], 'tiff');
                    
                    x_date = [0 15 31 30 31];
                    x_date = cumsum(x_date);
                    x_date_txt = {'Jul 16', 'Aug 1', 'Sep 1', 'Oct 1', 'Nov 1'};
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    fig_num = fig_num+1;
                    fh = figure(fig_num);
                    fh.WindowState = 'maximized';
                    num_fig_row = 2;
                    num_fig_col = 3;
                    clf;
                    clear yyy
                    tiledlayout(num_fig_row, num_fig_col,'TileSpacing','Compact','Padding','Compact');
                    
                    load('google_ads_stats');
                    for lang = lang_vec
                        nexttile
                        hold on
                        yyaxis left
                        plot(Google_Ads_stats(:,lang));
                        %                         set(gca,'XTick',x_date', 'XTickLabel', x_date_txt);
                        set(gca,'XTick',x_date', 'XTickLabel', []);
                        if (lang==TRAD_CHINESE)
                            yv = 0:1000:5000;
                        elseif (lang==SIMP_CHINESE)
                            yv = 0:200:1000;
                        else
                            yv = 0:100:500;
                        end
                        set(gca,'YTick',  yv, 'YTickLabel', yv);
                        if (lang==TRAD_CHINESE)
                            ylabel('No. of page requests');
                        end
                        ylim(yv([1 end]));
                        
                        yyaxis right
                        yv = 0:0.05:0.25;
                        ylim(yv([1 end]));
                        
                        if (lang==ENGLISH)
                            ylabel('Clickthrough rate');
                            set(gca,'YTick',  yv, 'YTickLabel', yv);
                        else
                            set(gca,'YTick',  yv, 'YTickLabel', []);
                        end
                        plot(Google_Ads_stats(:,3+lang)/100);
                        xlim(x_date([1 end]));
                        %                         title([language_text{lang}]);
                        xtickangle(45);
                        num_subjects = sum(data_for_inference.N(1:num_categories, end, ALL, end, ds_idx, time_window_idx, lang), 'all');
                        title({language_text{lang}; [' (n = ' num2str(num_subjects) ')']});
                        
                    end
                    
                    x_date = [0 15 31 30 31]*day_to_min;
                    x_date = cumsum(x_date);
                    for lang = lang_vec
                        clear mmm;
                        nexttile
                        hold on;
                        total_num_visits = 0;
                        for ii = 1:numel(category_display_order)
                            cat_idx = category_display_order(ii);
                            yyy{cat_idx} = [];
                            for t_idx=1:numel(time_partition)-1
                                yyy{cat_idx} = [yyy{cat_idx} time_visit_dist{cat_idx, t_idx, end, end, ds_idx, lang}*day_to_min];
                            end
                            yyy{cat_idx} = sort(yyy{cat_idx});
                            num_visits = numel(yyy{cat_idx});
                            [hhh xxx] = ecdf(yyy{cat_idx});
                            %                             zz = plot(xxx, hhh*num_visits);
                            zz = plot(xxx, hhh);
                            lt{cat_idx} = [category_full_text{cat_idx} ' (n = ' num2str(num_visits) ')'];
                            color_cat(cat_idx,:) = get(zz, 'Color');
                            total_num_visits = total_num_visits+num_visits;
                        end
                        
                        xlabel([]);
                        set(gca,'XTick',x_date', 'XTickLabel', x_date_txt);
                        if (lang==TRAD_CHINESE)
                            %                             set(gca,'YTickLabel',[]);
                            ylabel({'Cumulative proportion', 'of nudges'});
                        end
                        xlim(x_date([1 end]));
                        lh = legend(lt{category_display_order}, 'Location', 'Southeast', 'FontSize', 8);
                        lh.FontSize = 11;
                        legend boxoff
                        xtickangle(45);
                    end
                    saveas(gcf, ['Figures/figS1' suffix]);
                    saveas(gcf, ['Figures/figS1' suffix], 'tiff');
                    
                    
                    fig_num = fig_num+1;
                    fh = figure(fig_num);
                    fh.WindowState = 'maximized';
                    num_fig_row = 4;
                    num_fig_col = 3;
                    clf;
                    clear yyy
                    
                    th = tiledlayout(num_fig_row, num_fig_col,'TileSpacing','Compact','Padding','Compact');
                    zz =  get(th, 'Position');
                    zz(1) = zz(1)+0.02;
                    zz(3) = zz(3)-0.02;
                    zz(4) = zz(4)-0.08;
                    set(th, 'Position', zz)
                    an_left = repmat(0.013,1, 4);
                    an_bottom = [0.85 0.635 0.42 0.21];
                    an_str = {'A', 'B', 'C', 'D'};
                    for iii=1:numel(an_left)
                        annotation(gcf,'textbox',...
                            [an_left(iii) an_bottom(iii) 0.03203125 0.0455868089233754],...
                            'String',an_str(iii),...
                            'LineStyle','none',...
                            'FontWeight','bold',...
                            'FontSize',14,...
                            'FontName','Arial');
                    end
                    
                    xv = 1:(numel(time_partition)-1);
                    xvl = {'Jul', 'Aug', 'Sep', 'Oct'};
                    cw = 0.3;
                    cat_dx = linspace(-cw, cw, num_categories);
                    
                    for lang = lang_vec
                        nexttile
                        hold on;
                        for ii = 1:numel(category_display_order)
                            cat_idx = category_display_order(ii);
                            YY =   data_for_inference.N(cat_idx, 1:end-1, vax, age_idx, ds_idx, time_window_idx, lang);
                            YY = squeeze(YY);
                            xxx = (1:numel(YY))+cat_dx(ii);
                            bar(xxx, YY,  2*cw/num_categories, 'FaceColor', color_cat(cat_idx,:));
                        end
                        if (lang==1)
                            lh = legend(category_full_text(category_display_order), 'Location', 'North');
                            pos = get(lh, 'Position');
                            lh.NumColumns = 7;
                            lh.FontSize = 12;
                            lh.Orientation = 'Horizontal';
                            
                            legend boxoff
                            set(lh, 'Position',[0.16 1.01 1.2 0.04]);
                            
                        end
                        xvl = [];
                        set(gca, 'XTick', xv, 'XTickLabel', xvl);
                        if (lang==TRAD_CHINESE)
                            yv = [0 3.5 7];
                            ylabel({'No. of nudges'});
                        elseif (lang==SIMP_CHINESE)
                            yv = [0 0.35 0.7];
                        else
                            yv = [0 0.7  1.4];
                        end
                        yv = yv*1e3;
                        set(gca, 'YTick', yv, 'YTickLabel', yv);
                        ylim(yv([1 end]));
                        xlim([xv(1)-0.5 xv(end)+0.5]);
                        xtickangle(45);
                        %                         title([language_text{lang}]);
                        num_subjects = sum(data_for_inference.N(1:num_categories, end, ALL, end, ds_idx, time_window_idx, lang), 'all');
                        title({language_text{lang}; [' (n = ' num2str(num_subjects) ')']});
                        
                    end
                    
                    
                    cw = 0.3;
                    cat_dx = linspace(-cw, cw, num_categories);
                    for lang = lang_vec
                        nexttile
                        hold on;
                        for ii = 1:numel(category_display_order)
                            cat_idx = category_display_order(ii);
                            YY =   data_for_inference.Y(cat_idx, :, vax, age_idx, ds_idx, time_window_idx, lang);
                            YY = squeeze(YY);
                            xxx = (1:numel(YY))+cat_dx(ii);
                            bar(xxx, YY,  2*cw/num_categories, 'FaceColor', color_cat(cat_idx,:));
                        end
                        
                        xvl = [];
                        set(gca, 'XTick', xv, 'XTickLabel', xvl);
                        if (lang==TRAD_CHINESE)
                            ylabel('No. of bookings');
                        end
                        %                         ylim(yv([1 end]));
                        xlim([xv(1)-0.5 xv(end)+0.5]);
                        xtickangle(45);
                    end
                    
                    %%%%%%%% Plot age distribution of bookings
                    for lang = lang_vec
                        nexttile
                        hold on;
                        for ii = 1:numel(category_display_order)
                            cat_idx = category_display_order(ii);
                            YY = data_for_inference.Y(cat_idx, :, vax, 1:end-1, ds_idx, time_window_idx, lang);
                            YY = squeeze(YY);
                            YY = YY./repmat(sum(YY,2),1,size(YY,2));
                            xxx = (1:size(YY,1))+cat_dx(ii);
                            
                            b = bar(xxx, YY, 2*cw/num_categories, 'stacked');
                            for ii3=1:numel(b)
                                b(ii3).FaceColor = color_cat(cat_idx,:);
                                b(ii3).FaceAlpha = ii3/numel(b);
                            end
                        end
                        xvl = [];
                        set(gca, 'XTick', xv, 'XTickLabel', xvl);
                        if (lang==1)
                            ylabel({'Proportion', 'of bookings'});
                        end
                        %                         ylim(yv([1 end]));
                        xlim([xv(1)-0.5 xv(end)+0.5]);
                        xtickangle(45);
                    end
                    
                    %%%%%%%% Plot expected number of bookings per nudge
                    for lang = lang_vec
                        nexttile
                        hold on;
                        
                        for ii = 1:numel(category_display_order)
                            cat_idx = category_display_order(ii);
                            mu_cat = mu_mat(cat_idx, 1:end-1, vax, age_idx, ds_idx, time_window_idx, lang);
                            sigma_cat = sigma_mat(cat_idx, 1:end-1, vax, age_idx, ds_idx, time_window_idx, lang);
                            YY = norminv([0.5 0.025 1-0.025]', mu_cat, sigma_cat);
                            xxx = (1:size(YY,2))+cat_dx(ii);
                            bar(xxx, YY(1,:),  2*cw/num_categories, 'FaceColor', color_cat(cat_idx,:));
                            plot([xxx; xxx], YY(2:3,:), 'k');
                        end
                        xvl = {'July', 'August', 'September', 'October'};
                        set(gca, 'XTick', xv, 'XTickLabel', xvl);
                        yv =0:0.25:0.5;
                        if (lang==1)
                            set(gca, 'YTick', yv, 'YTickLabel', yv);
                            ylabel({'Expected no. of', 'bookings per nudge'});
                        else
                            set(gca, 'YTick', yv, 'YTickLabel', yv);
                        end
                        ylim(yv([1 end]));
                        xlim([xv(1)-0.5 xv(end)+0.5]);
                        %                         xtickangle(45);
                    end
                    saveas(gcf, ['Figures/fig1' suffix]);
                    saveas(gcf, ['Figures/fig1' suffix], 'tiff');
                    
                end
            end
        end
    end
end
% save(['Results/script_3_output_split' num2str(margin_for_session_split)]);
disp('Finished running script 3.');

