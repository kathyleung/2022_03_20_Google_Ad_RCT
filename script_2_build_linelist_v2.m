clear all;
close all;

SINOVAC = 1;
BNT = 2;
ALL = 3;
vax_text = {'Sinovac', 'BioNTech', 'Both'};

startDate = '2021-07-16 00:00:01';
startDate_2 = '2021-07-15 00:00:01';
endDate_2 = '2021-11-05 23:59:59';
endDate = '2021-11-01 00:00:01';

AA = readtable('Data/rev_data_matlab_no_merging');  
master_data = sortrows(AA, [1 2 3 4]);
fI = (master_data.time_visit_rev>=startDate & master_data.time_visit_rev<=endDate);
master_data = master_data(fI,:);
time_start = startDate;
time_end = endDate;
%
fI_click = find(master_data.status_click);
master_data(fI_click,:) = [];   % remove all visits with clicks

languageText = {'tc', 'en', 'sc'};

figure(1);
clf;
for lang=1:3
    subplot(3,1,lang);
    fI = find(strcmp(master_data.landing_page_language, languageText{lang}));
    histogram(master_data.category_idx(fI,:),0.5:1:7.5);
end

day_to_min = 60*24;

%%% Building the line list
AA = master_data;

%%%% v1 = do not stratify by language at this stage
%%%% v2 = stratify by language
for lang=1:3
    clear linelist
    ip_idx = 0;
    prev_ip = '';
    
    fI = find(AA.time_visit_rev > time_start & AA.time_visit_rev < time_end & strcmp(languageText{lang},AA.landing_page_language));
    fI = fI';
    tic
    for ii=fI
        this_time_visit = (datenum(AA.time_visit_rev(ii))-datenum(time_start))*day_to_min;
        this_ip = cell2mat(AA.ip_v4(ii));
        this_booking_id = AA.booking_system_id(ii);
        
        if (~strcmp(this_ip, prev_ip))
            ip_idx = ip_idx+1;
            linelist(ip_idx).ip = this_ip;
            linelist(ip_idx).ip_idx = ip_idx;
            linelist(ip_idx).time_visit = [];
            linelist(ip_idx).category = [];
            linelist(ip_idx).language = [];
            linelist(ip_idx).click_button = [];
            linelist(ip_idx).booking_id = [];
            linelist(ip_idx).booking_time = [];
            linelist(ip_idx).booking_age = [];
            
        end
        if (~ismember(this_time_visit, linelist(ip_idx).time_visit))    % new IP address; these are information about the page visit
            nn = numel(linelist(ip_idx).time_visit)+1;    % next index
            
            linelist(ip_idx).time_visit(nn) = this_time_visit;
            linelist(ip_idx).category(nn) = AA.category_idx(ii);
            linelist(ip_idx).click_button(nn) = AA.status_click(ii);
            
            linelist(ip_idx).language(nn) = -1;
            lang_str = cell2mat(AA.landing_page_language(ii));
            for jj=1:numel(languageText)
                if (strcmp(lang_str, languageText{jj}))
                    linelist(ip_idx).language(nn) = jj;
                end
            end
        end
        if (this_booking_id>0 && ~ismember(this_booking_id, linelist(ip_idx).booking_id))    % new booking; these are information about the bookings associated with the page visit
            nn = numel(linelist(ip_idx).booking_id)+1;
            linelist(ip_idx).booking_id(nn) = this_booking_id;
            linelist(ip_idx).booking_time(nn) =  this_time_visit+AA.time_elapsed(ii);
            if (strcmp(AA.booking_vax(ii), 'BioNTech/Fosun'))
                this_vax = BNT;
            elseif (strcmp(AA.booking_vax(ii), 'Sinovac'))
                this_vax = SINOVAC;
            else
                error('Error in booking vax type.');
            end
            linelist(ip_idx).booking_vax(nn) = this_vax;
            linelist(ip_idx).vax_time(1,nn) =  (datenum(AA.booking_first_dose_date_time(ii))-datenum(time_start))*day_to_min-this_time_visit;
            linelist(ip_idx).vax_time(2,nn) =  (datenum(AA.booking_second_dose_date_time(ii))-datenum(time_start))*day_to_min-this_time_visit;
            linelist(ip_idx).booking_age(nn) = 2021-AA.booking_year_birth(ii);
            
        end
        
        prev_ip = this_ip;
        t2 = toc;
        if (mod(ii,100)==0)
            ii
            disp(['Time elapsed = ' num2str(t2/60) ' min. Time remaining = ' num2str(t2/ii*(numel(fI)-ii)/60) ' min']);
        end
    end
    linelist_large{lang} = linelist;
end

save('Results/script_2_output_v2');
load('Results/script_2_output_v2');

linelist = linelist_large;
for lang=1:3
    linelist_master{lang} = rmfield(linelist{lang}, {'ip', 'booking_id', 'click_button'});
end
clear booking_data AA XX master_data fI this_ip prev_ip this_booking_id linelist_large linelist
save('Results/data_final');

