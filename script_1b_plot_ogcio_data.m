clear all;
close all;

load('script_1_output');   

time_window = 60;
time_window_idx = time_window_indices(time_window);

figure(1);
clf
subplot(2,1,1);
plot(tv_num_bookings{time_window_idx}/day_to_min, num_bookings{time_window_idx});
xlabel('Time since July 16 (min)');
ylabel(['No. of bookings (' num2str(time_window) '-min time window)']);
legend('Sinovac', 'BioNTech', 'Total')
subplot(2,1,2);
plot(tv_num_bookings{time_window_idx}/day_to_min, num_bookings{time_window_idx}(1,:)./num_bookings{time_window_idx}(2,:));
xlabel('Time since July 16 (min)');
ylabel('Sinovac-to-BioNTech ratio');
xv = [0 16 16+15 16+15+16];


figure(2)
clf;
subplot(2,1,1);
plot(tv_booking_ratio{time_window_idx}, min(100,booking_ratio{time_window_idx}))
ylabel('Booking ratio');
xlabel('Time since July 16 (min)');
subplot(2,1,2);
hold on;
for ii=1:3
    histogram(booking_ratio{time_window_idx}(ii,:),0:0.1:10);
end
xlabel('Booking ratio');

% for ii=1:size(AA,1)
%     ii
%     GG{ii} = strcat(mat2str(cell2mat(AA.ip_v4(ii))), datestr(AA.time_visit_rev(ii)));
% end
