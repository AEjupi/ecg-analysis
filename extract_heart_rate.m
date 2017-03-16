% Author: Andreas Ejupi (andreas@ejupi.at)

%% data loading

clear; close all; 
path = '/Users/Zero/Desktop/heart_rate_data/young_adults/20170203_MXT_Session1_Sternum_Calibrated_SD.csv'; 
df_orig = dlmread(path, '\t', 5, 0);

%% ecg analysis

% define columns/ sampling rate
Fs = 128; %Hz

XACC = 2; 
YACC = 3; 
ZACC = 4;

VX = 10;
LA_RA = 7;
LL_LA = 8;
LL_RA = 9;


% plot acceleration for selection
x_acc = df_orig(:,XACC); 
y_acc = df_orig(:,YACC);
z_acc = df_orig(:,ZACC);

sv_acc = sqrt(x_acc.^2 + y_acc.^2 + z_acc.^2); 

ts = 1:length(sv_acc);
ts = ts/Fs/60; 
      
figure; hold on; grid;
plot(ts, sv_acc);
title('Acceleration');
ax = gca; 
ax.XTick = round(ts(1:128*30:end),1);


% trim to selected data 
[x,~] = ginput(2);

start_ix = x(1)*60*Fs; 
stop_ix = x(2)*60*Fs;

df = df_orig(start_ix:stop_ix, :);


% select best lead

% smooth data
[b,a] = butter(2,0.5/64,'high'); 

%df(:,LA_RA) = filtfilt( b, a, df(:,LA_RA) );
%df(:,LL_LA) = filtfilt( b, a, df(:,LL_LA) );
%df(:,LL_RA) = filtfilt( b, a, df(:,LL_RA) );

unipolar_ecg = df(:,VX); % unipolar lead Vx
averaged_ecg = (df(:,LA_RA) + df(:,LL_LA) + df(:,LL_RA)) / 3; 


ts = 1:length(averaged_ecg);
ts = ts/Fs/60; 


figure; 

sp(1) = subplot(4,1,1); hold on; grid;
plot(ts, df(:,LA_RA) ); 
title('Lead LA-RA'); 
%ylim([-1 1]);  

sp(2) = subplot(4,1,2); hold on; grid;
plot(ts, df(:,LL_LA) ); 
title('Lead LL-LA'); 
%ylim([-1 1]);  

sp(3) = subplot(4,1,3); hold on; grid;
plot(ts, df(:,LL_RA) ); 
title('Lead LL-RA'); 
%ylim([-1 1]);  

sp(4) = subplot(4,1,4); hold on; grid;
plot(ts, averaged_ecg ); 
title('Averaged Lead 1-3');
%ylim([-1 1]);  

% sp(5) = subplot(5,1,5); hold on; grid;
% plot(ts, unipolar_ecg-mean(unipolar_ecg) ); 
% title('Unipolar Vx');
% ylim([-2 2]); 

[~,~] = ginput(1);
clicked_sp = gca; 

switch (clicked_sp.Title.String)
    case 'Lead LA-RA'
        ecg = df(:,LA_RA); 
    case 'Lead LL-LA'
        ecg = df(:,LL_LA); 
    case 'Lead LL-RA'
        ecg = df(:,LL_RA); 
    case 'Averaged Lead 1-3'
        ecg = averaged_ecg; 
    case 'Unipolar Vx'
        ecg = unipolar_ecg; 
end


% heart rate analysis for selected data

figure; 
sp(1) = subplot(3,1,1); hold on; grid;
plot(ts, ecg);
title(clicked_sp.Title.String);

wt = modwt(ecg,5);
wtrec = zeros(size(wt)); 
wtrec(2:3,:) = wt(2:3,:);
y = imodwt(wtrec,'sym4');

sp(2) = subplot(3,1,2); hold on; grid;
plot(ts, y, 'r');
title('Wavelet-based filtered signal');

sp(3) = subplot(3,1,3); hold on; grid;
plot(ts, y, 'g');
title('R peaks (squared absolute values)');

[~,y_th] = ginput(1);
%y = (y+1).^2;
y(y < 0) = 0; 
[qrspeaks,ix] = findpeaks(y,ts,'MinPeakHeight', y_th,...
    'MinPeakDistance', 0.20/60);

plot(ix,qrspeaks,'ro');
h = brush; 
h.Enable = 'on'; 

linkaxes(sp,'x');

%% heart rate calculation

% remove selected points (brushing)
if exist('remove') && ~isempty(remove) 
    for i=1:length(remove(:,1))
        ix(ix == remove(i,1)) = [];
    end
end
            
heart_rate = length(qrspeaks) / ts(end);
heart_rate_var = std(diff(ix)*60*1000); 
heart_rate_mad = mad(diff(ix)*60*1000); 

fprintf('Duration of activity: %.2f min\n', ts(end));
fprintf('Average Heart Rate: %d BPM\n', round(heart_rate));
fprintf('Heart Rate Variability (std): %d ms\n',round(heart_rate_var));
fprintf('Heart Rate Variability (mad): %d ms\n',round(heart_rate_mad));


