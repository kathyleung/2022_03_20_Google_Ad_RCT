%%%%%% This is the first script of data processing series
%%%%%% This script load the OGCIO vax booking database

clear all;
close all;

SINOVAC = 1;
BNT = 2;
ALL = 3;
vax_text = {'Sinovac', 'BioNTech', 'Both'};

startDate = '2021-07-16 00:00:01';
startDate_2 = '2021-07-15 00:00:01';
endDate_2 = '2021-11-05 23:59:59';
endDate = '2021-10-31 23:59:59';

booking_data = [];
for ii=1:39
    ii
    if (ii<10)
        data_file_name = ['Data/Booking_Statistic_HKU_0000' num2str(ii) '.csv'];
    else
        data_file_name = ['Data/Booking_Statistic_HKU_000' num2str(ii) '.csv'];
    end
    AA = readtable(data_file_name);
    AA = AA(:,[3 6 7]);
    booking_data = [booking_data; AA];
end

fI{SINOVAC} = find(strcmp(booking_data.VaccineBooked, 'Sinovac'));
fI{BNT} = find(strcmp(booking_data.VaccineBooked, 'BioNTech/Fosun'));
if (numel(fI{BNT})+numel(fI{SINOVAC})~=size(booking_data,1))
    error('There are bookings that are neither Sinovac or BioNTech.');
end

save('Results/script_1_output');
return

day_to_min = 60*24;
x_date = [0 15 31 30 31]*day_to_min;
x_date = cumsum(x_date);
x_date_txt = {'16/7', '1/8', '1/9', '1/10', '1/11'};

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
max_time_window = 60;
time_window_vec = dt:dt:max_time_window;
time_window_indices(time_window_vec) = 1:numel(time_window_vec);
tv = floor(timeMin):ceil(timeMax);
booking_ratio = zeros(numel(tv),num_age_groups,3,numel(time_window_vec));
booking_ratio_v2 = zeros(numel(tv),num_age_groups,3,numel(time_window_vec));

clear booking_data AA
save('Results/script_1_output');

return

max_tw_margin = zeros(num_age_groups, 3, numel(time_window_vec));
for age_idx = 1:num_age_groups
    for vax=[SINOVAC BNT ALL]
        if (vax==ALL)
            num_bookings(:,age_idx,vax) = histcounts(booking_time{age_idx,SINOVAC},tv)+histcounts(booking_time{age_idx,BNT}, tv);
        else
            num_bookings(:,age_idx,vax) = histcounts(booking_time{age_idx,vax}, tv);
        end
        if (all(num_bookings(:,age_idx,vax)==0))
            continue;
        end
        YY = num_bookings(:,age_idx,vax);
        for time_window = fliplr(time_window_vec)
            [age_idx vax time_window]
            time_window_idx = time_window_indices(time_window);
       
            for ii=1:numel(YY)-2*day_to_min
                ii_idx = ii-tv(1)+1;
                tw = time_window;
                
                numer = 1+sum(YY(ii_idx:ii_idx+tw));
                denom = 1+sum(YY(ii_idx-tw:ii_idx));
                booking_ratio(ii,age_idx,vax,time_window_idx) = numer./denom;
                
                denom = sum(YY(ii_idx-tw:ii_idx));
                while (denom==0)
                    tw = tw+1;
                    denom = sum(YY(ii_idx-tw:ii_idx));
                end
                numer = sum(YY(ii_idx:ii_idx+tw));
                max_tw_margin(age_idx,vax,time_window_idx) = max(tw-time_window, max_tw_margin(age_idx,vax,time_window_idx));
                booking_ratio_v2(ii,age_idx,vax,time_window_idx) = numer./denom;
            end
        end
    end
end

clear num_bookings
ylabel_txt = {'No. of bookings', 'Booking ratio'};
for plot_num=1:2
    figure(plot_num);
    clf;
    time_window = 60;
    time_window_idx = time_window_indices(time_window);
    tv = floor(timeMin):time_window:ceil(timeMax);
    s_idx = 1;
    for age_idx = 1:num_age_groups
        for vax=[SINOVAC BNT]
            num_bookings(:,age_idx,vax) = histcounts(booking_time{age_idx,vax}, tv);
            subplot(num_age_groups,2,s_idx);
            if (plot_num==1)
                plot(tv(1:end-1),num_bookings(:,age_idx,vax))
            else
                plot(booking_ratio(:,age_idx,vax,time_window_idx))
            end
            s_idx = s_idx+1;
            title([num2str(age_groups(age_idx,1)) '-' num2str(age_groups(age_idx,2)) ' ' vax_text{vax}]);
            ylabel(ylabel_txt{plot_num});
            set(gca,'XTick',x_date', 'XTickLabel', x_date_txt);
            xlim(x_date([1 end]));
        end
    end
end





