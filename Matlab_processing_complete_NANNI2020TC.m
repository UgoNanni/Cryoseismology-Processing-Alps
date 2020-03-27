
%***********************************************************
% This code is associated to the paper 
% 'Quantification of seasonal and diurnal dynamics of subglacial channels using seismic observations on an Alpine Glacier.' Nanni et al. 2020
% please cite accordingly
% author: Ugo Nanni, 27/03/2020
% ugo.nanni@univ-grenoble-alpes.fr
%***********************************************************
%% Load data
clear all
addpath 'C:\Users\nanni\Documents\PhD\Processing\Matlab_Codes\scripts_MATLAB'
%
close all
beep off

% parameters for the seismic power
f_res = 0.5;
bfsup = 7/f_res+1;
bfinf = 3.5/f_res+1;
bfsupqs = 20/f_res+1;
bfinfqs = 12/f_res+1;

f30sec = 0;
f15min = 1;

if f15min == 1
    type = '15min';
elseif f30sec == 1
    type = '1to50Hz';
end

dofilter = 1; % filter of bin the signals
domean = 0;

%% Load dataset *********************************************************************************************************************************************************

if f15min == 1
    cd 'C:\Users\nanni\Documents\PhD\DATA\Pw_BF_4sec_sorted'% !!!!! Here Pw is in dB !!!!!!!!!!
    load('ARG_Pw_comp_15min')
    load('ARG_Tdoy_comp_15min')
    Pw_comp = ARG_Pw_comp_15min;
    T_pwcomp = ARG_Tdoy_comp_15min;
    
elseif f30sec == 1
    cd 'C:\Users\nanni\Documents\PhD\DATA\Pw_BF_4sec_sorted'% !!!!! Here Pw is in dB !!!!!!!!!!
    load('ARG_Pw_comp_30sec')
    load('ARG_Tdoy_comp_30sec')
    Pw_comp = ARG_Pw_comp_30sec;
    T_pwcomp = ARG_Tdoy_comp_30sec;
end
% !!!!! Here Pw is in dB !!!!!!!!!!
Pwater = median(Pw_comp(bfinf:bfsup,:));
madPwater = mad(Pw_comp(bfinf:bfsup,:),1);
Psed = median(Pw_comp(bfinfqs:bfsupqs,:));
madPsed = mad(Pw_comp(bfinfqs:bfsupqs,:),1);

cd 'C:\Users\nanni\Documents\PhD\DATA\Paper_Data\Discharge'
load('ARG_discharge.mat')
load('ARG_Tdoy_discharge.mat')
load('ARG_qs_raw_m3psec.mat')
load('ARG_Tqs_day.mat')

Qlim = 0.0109*7; % detection limit on the discharge
Qlim = 0.4;

cd 'C:\Users\nanni\Documents\PhD\DATA\DATA_ARG_complementary\ARG_Basal_waterpressure'
load('ARG_pressure.mat')
load('ARG_Tdoy_pressure.mat')

cd 'C:\Users\nanni\Documents\PhD\DATA\DATA_ARG_complementary\ARG_Meteo_Argentiere'
load('ARG_precip.mat')
load('ARG_mto_time.mat')
load('ARG_temp.mat')
ARG_Tdoy_meteo = ARG_mto_time;

cd 'C:\Users\nanni\Documents\PhD\DATA\Paper_Data\Sliding'
load('ARG_sliding_time_5sec_secfrom2017.mat')% !!!!! Here in sec !!!!!!!!!!
ARG_Tdoy_vel= ARG_sliding_time_5sec_secfrom2017/3600/24; % now in DOY
load('ARG_sliding_vel_mmperhour_dt5sec.mat')


cd   'C:\Users\nanni\Documents\PhD\DATA\DATA_ARG_complementary\Vitesse ARG\new_whole2002019'
load('Vslid_whole.mat'); load('dayofyear_whole.mat')
datwhole = datwhole +4;

%% Color schema *********************************************************************************************************************************************************
% figure
RGB = [49 163 84;161 217 155;49 130 189;158 202 225;222 45 38; 252 146 114;117 107 177;188 189 220;99 99 99; 189 189 189];
% for it=1:length(RGB)
%
%     hold on
%     plot(1:10,it*ones(10,1),'-','Color',RGB(it,:)/256,'LineWidth',10)
%
% end
colormap parula
%% Apply similar movmean *********************************************************************************************************************************************************
fs = median(diff(T_pwcomp)); % seismic sampling rate
%dt = 6/24; % in days
dt = 1/96; % in days = 15min
Vel = movmean(ARG_sliding_vel_mmperhour_dt5sec,dt,'omitnan','SamplePoints',ARG_Tdoy_vel);
Velwhole = movmean(velwhole/24*10,dt,'omitnan','SamplePoints',datwhole);
Press = movmean(ARG_pressure,dt,'omitnan','SamplePoints',ARG_Tdoy_pressure);
Temp = movmean(ARG_temp,dt,'omitnan','SamplePoints',ARG_Tdoy_meteo);
Precip = movmean(ARG_precip,dt,'omitnan','SamplePoints',ARG_Tdoy_meteo);
Q = movmean(ARG_discharge,dt,'omitnan','SamplePoints',ARG_Tdoy_discharge);
Qs = movmean(ARG_qs_raw_m3psec,60,'omitnan','SamplePoints',ARG_Tqs_day);
Pwraw = movmean(Pwater,dt,'omitnan','SamplePoints',T_pwcomp);% dB
Pb = movmean(Psed,dt,'omitnan','SamplePoints',T_pwcomp);% dB


%% remove vel peaks *********************************************************************************************************************************************************
med_std = median(movstd(ARG_sliding_vel_mmperhour_dt5sec,dt,'SamplePoints',ARG_Tdoy_vel),'omitnan');
idx = find((movstd(ARG_sliding_vel_mmperhour_dt5sec,dt,'SamplePoints',ARG_Tdoy_vel))>10*med_std);
Vel(idx) = NaN;

tsample = ARG_Tdoy_vel;
day1 = [tsample(1) 268];
day2 = [214 293];
for dt = 1:length(day1)
    id1 = find(abs(tsample-day1(dt))==min(abs(tsample-day1(dt))));
    id2 = find(abs(tsample-day2(dt))==min(abs(tsample-day2(dt))));
    Vel(id1:id2)= NaN;
end



%% Winter noise *********************************************************************************************************************************************************

% Calculate winter mean - anthropogenic contribution
rm_linear = 0;
rm_daily = 1;
decnovoluate = 0;
subtract = 1;
% ------------------------------------------------------
sample = movmean(Pwater,2/24,'omitnan','SamplePoints',T_pwcomp);% dB;% dB
tsample = T_pwcomp;
% ------------------------------------------------------

day1 = 395;
day2 = 460;

for dt = 1:length(day1)
    id1 = find(abs(tsample-day1(dt))==min(abs(tsample-day1(dt))));
    id2 = find(abs(tsample-day2(dt))==min(abs(tsample-day2(dt))));
    Pwinter = sample(id1:id2);
    Twinter = tsample(id1:id2);
end


Twinter_interp = day1:1/96:day2; % regular winter grid with 15 min dt
Pwinter = interp1(Twinter,Pwinter,Twinter_interp); % interp the winter signal

% Go to linear space
Pwinter = (10.^(Pwinter/10));
Pwinter_stat = Pwinter; % interp the winter signal


idx_db = find(mod(Twinter_interp,1) == 0); % day beginnings
idx_db = idx_db(1:end-1);
for it = 0:95 % loop on one day
    
    Pwday_winter(it+1,1) = mean(Pwinter(idx_db(1:end-1)+it),'omitnan');
    Pwday_winter(it+1,2)  = std(Pwinter(idx_db(1:end-1)+it),'omitnan');
    Pwday_winter(it+1,3)  = median(Pwinter(idx_db(1:end-1)+it),'omitnan');
    Pwday_winter(it+1,4)  = mad(Pwinter(idx_db(1:end-1)+it));
    for jt = 1:length(idx_db)-1
        Pwinter_stat(idx_db(jt)+it) = Pwday_winter(it+1,3);
    end
    
end

% interpolare over the whole time serie
Tdaily_winter = floor(tsample(1)):1/96:ceil(tsample(end));
Pdaily_winter = nan(1,length(Tdaily_winter));
idx_db = find(mod(Tdaily_winter,1) == 0); % day beginnings
idx_db = idx_db(1:end-1);

for it = 0:95 % loop on one day
    for jt = 1:length(idx_db)-1
        Pdaily_winter(idx_db(jt)+it) = Pwday_winter(it+1,3);
    end
end

Pdaily_winter = interp1(Tdaily_winter,Pdaily_winter,tsample);

if rm_linear == 1    % Remove in linear the antrhopique contribution
    mw = mean(Pwinter,'omitnan');% + std(Pwinter);% dB
    Pwlin = (10.^(Pwraw/10)-mw);
elseif rm_daily == 1    % Remove the daily the antrhopique contribution
    
    Pwlin = (10.^(Pwraw/10));
    if subtract == 1 % subtraction only
        Pwlin = Pwlin-Pdaily_winter;
    elseif decnovoluate == 1 % deconvolution
        [q,r] = deconv(Pwlin,Pdaily_winter);
        Pwlin = r;
    end
end

Pw = real(10*log10(Pwlin)); % to dB
Pdaily_winter = real(10*log10(Pdaily_winter)); % dB
clear Pwlin

%% Interpolate all variables onto regular grid *********************************************************************************************************************************************************

% the seismic sets the reference time
% ----------------------------------------------
T_interp = floor(120):1/96:ceil(T_pwcomp(end));
Fs = 1/median(diff(T_interp)); % 15 min


% keep NAN in interp
removenan = 1;

Pw_interp = interp1(T_pwcomp,Pw,T_interp,'nearest');% dB
if removenan == 1 idx = find(gradient(Pw_interp)==0); Pw_interp(idx+1)=nan;Pw_interp(idx-1)=nan;  end% remove nan
Pw_interpraw = interp1(T_pwcomp,Pwater,T_interp,'nearest');% dB
if removenan == 1 idx = find(gradient(Pw_interpraw)==0); Pw_interpraw(idx+1)=nan;Pw_interpraw(idx-1)=nan;  end% remove nan
Pb_interp = interp1(T_pwcomp,Pb,T_interp,'nearest');% dB
if removenan == 1 idx = find(gradient(Pw_interp)==0); Pw_interp(idx)=nan; end% remove nan
%Pdaily_winter_interp = interp1(tsample,Pdaily_winter,T_interp);% dB
%Qs_interp = interp1(ARG_Tqs_day,Qs,T_interp);
Vel_interp = interp1(ARG_Tdoy_vel,Vel,T_interp);
Vel_interp = interp1(datwhole,Velwhole,T_interp);
Press_interp = interp1(ARG_Tdoy_pressure,Press,T_interp);
Q_interp = interp1(ARG_Tdoy_discharge,Q,T_interp);
Temp_interp = interp1(ARG_Tdoy_meteo,Temp,T_interp);
Precip_interp = interp1(ARG_Tdoy_meteo,Precip,T_interp);




% night times for T_interp
% ----------------------------------------------

idx_night = find(mod(T_interp,1) >= 0.8 | mod(T_interp,1)  <= 0.2);
idx_night = find(T_interp>0);
idx_Qlim = find(Q_interp>0.01*5);
idx_vel = idx_night;

% loagrithmic values
% ----------------------------------------------
Q_interplog = 10*log10(Q_interp);
Vel_interplog = 10*log10(Vel_interp);

idx = find(isfinite(Vel_interplog)==1);
Vel_interplog = interp1(T_interp(idx),Vel_interplog(idx),T_interp);
idxnanvelN = find(isnan(Vel_interplog)==0);
idxnanvelY = find(isnan(Vel_interplog)==1);

Press_interplog = 10*log10(Press_interp);
idx = find(isfinite(Press_interplog)==1);
Press_interplog = interp1(T_interp(idx),Press_interplog(idx),T_interp);
idxnanpressN = find(isnan(Press_interplog)==0);
idxnanpressY = find(isnan(Press_interplog)==1);

%%  Melt seasons periods *********************************************************************************************************************************************************


% ------------------------------------------------------

day1 = [130 135 120+365 147+365]; %season and day 2017 season and day 2018
day2 = [315 264 315+365 306+365];

day1 = [132 135 110+365 142+365]; %season and day 2017 season and day 2018
day2 = [297 264 315+365 300+365];

idx17s = find(T_interp >= day1(1) & T_interp <= day2(1));
idx17d = find(T_interp >= day1(2) & T_interp <= day2(2));
idx18s = find(T_interp >= day1(3) & T_interp <= day2(3));
idx18d = find(T_interp >= day1(4) & T_interp <= day2(4));

T17s = T_interp(idx17s);
T17d = T_interp(idx17d);
T18s = T_interp(idx18s);
T18d = T_interp(idx18d);

idx_season = [idx17s idx18s];
idx_days = [idx17d idx18d];



%% Spectral filters Low pass *********************************************************************************************************************************************************
% Define references periods

% Low pass for seasons
% ---------------------------------------------------------------------
order = 2;
Fcutlow_lp = [Fs/24 2 1 1/5 1/10 1/30 1/60 24];

clear Qlog & Pwlog & Vlog & Presslog
for it = 1:length(Fcutlow_lp)
    fcutlow = Fcutlow_lp(it);
    for var = 1:8
        
        if var == 1 tmp = Q_interplog; end
        if var == 2 tmp = Pw_interp; end
        if var == 8 tmp = Pw_interpraw; end
        if var == 3 tmp = Vel_interplog; end
        if var == 4 tmp = Press_interplog; end
        if var == 5 tmp = Vel_interp; end
        if var == 6 tmp = Q_interp; end
        if var == 7 tmp = 10.^(Pw_interp/10); end
        
        idx_night = idx_season;
        if var == 5 idx_night = idx_vel; end
        if var == 3 idx_night = idx_vel; end
        
        idxnan = find(isnan(tmp) == 1);
        tmp = tmp(idx_night); % taking only night events
        idx = isfinite(tmp); % compute for not NaN values
        tmp = tmp(idx); m = mean(tmp,'omitnan');
        tmp = (tmp-m); % remove mean, do not detrend
        [b,a] = butter(order,fcutlow/(Fs/2), 'low'); % low pass
        
        if dofilter == 1
            tmp = filtfilt(b,a,tmp);
        elseif domean == 1
            tmp = movmean(tmp,1/fcutlow,'omitnan','SamplePoints',T_interp);
        end
        
        tmp = interp1(T_interp(idx_night(idx)),tmp+m,T_interp);
        tmp(idxnan) = nan;
        % cut first and last day of the season
        if var ~= 5 && var ~= 3
            idx = find(T_interp > day2(1)-1 & T_interp < day1(3)+2 );tmp(idx) = nan;
        end
        
        if var == 1 Qlog(:,it) = tmp; end
        if var == 2 Pwlog(:,it) = tmp; end
        if var == 8 Pwlograw(:,it) = tmp; end
        if var == 7 Pwl(:,it) = tmp; end
        if var == 3 Vlog(:,it) = tmp; end
        if var == 4 Presslog(:,it) = tmp; end
        if var == 6 Ql(:,it) = tmp; end
        if var == 5 Vl(:,it) = tmp; end
    end
end

%% R and S inversion *********************************************************************************************************************************************************

% invert S and R from filtered Q and P to avoid short periods influences on
% the inversions
clear Qlin & Pwlin & Pwlin & Pwref & Qref & R & S & Srothi
order = 2;
idx = idx17s;
for it = 1:length(Fcutlow_lp)
    
    Qlin(:,it) = 10.^(Qlog(:,it)/10);
    Pwlin(:,it) = 10.^(Pwlog(:,it)/10);
    %Pwlin(:,it) = 10.^(Pwlograw(:,it)/10);
    
    Qlin(:,it) = Ql(:,it);
    %Pwlin(:,it) = Pwl(:,it);
    %
    %
    Pwref(it) = nanmin(Pwlin((idx(1)),it));
    Qref(it) = nanmin(Qlin((idx(1)),it));
    
    % inversions
    R(:,it) = ((Pwlin(:,it)/Pwref(it)).^(-9/82)).*((Qlin(:,it)/Qref(it)).^(21/41));
    S(:,it) = ((Pwlin(:,it)/Pwref(it)).^(24/41)).*((Qlin(:,it)/Qref(it)).^(-30/41));
    Srothi(:,it) = ((Qlin(:,it)/Qref(it)).^(2));
    Srothi_eq(:,it) = ((Qlin(:,it)/Qref(it)).^(-2/11));
    
end
% normalize
for it = 1:length(Fcutlow_lp)
    
    id = find(S==min(S(idx)));
    S(:,it) = S(:,it);%./min(S(idx(:),6));
    R(:,it) = R(:,it);%./min(S(idx(:),6));
    Srothi(:,it) = Srothi(:,it)./min(Srothi(idx(1),6));
    
end

%% Spectral filters Band pass *********************************************************************************************************************************************************

% Select the diurnal pattern
% take the season period

Fcutlow_bp = [2/3 2/4 2/3 ];
Fcuthigh_bp = [Fs/48 Fs/3 Fs/3];

% Band pass filter for the diurnal signals
clear Qdayf & Pdayf & Vdayf & Pressdayf

for it = 1:length(Fcutlow_bp)
    for var = 1:7
        
        switch var
            case 1; tmp = Q_interp;
            case 7; tmp = Q_interplog;
            case 2; tmp = Pw_interp;
            case 3; tmp = Vel_interp;
            case 4; tmp = Press_interplog;
            case 5; tmp = R(:,1);
            case 6; tmp = S(:,1);
        end
        
        tmp = tmp(idx_night);
        idx = isfinite(tmp);
        tmp = tmp(idx);
        
        tmp = (tmp-mean(tmp,'omitnan')); % do not detrend
        
        [b,a] = butter(order,[Fcutlow_bp(it) Fcuthigh_bp(it) ]/(Fs/2), 'bandpass');
        tmp = filtfilt(b,a,tmp);
        
        switch var
            case 1; Qdayf(:,it) = interp1(T_interp(idx_night(idx)),tmp,T_interp);
            case 2; Pdayf(:,it) = interp1(T_interp(idx_night(idx)),tmp,T_interp);
            case 3; Vdayf(:,it) = interp1(T_interp(idx_night(idx)),tmp,T_interp);
            case 4; Pressdayf(:,it) = interp1(T_interp(idx_night(idx)),tmp,T_interp);
            case 5; Rdayf(:,it) = interp1(T_interp(idx_night(idx)),tmp,T_interp);
            case 6; Sdayf(:,it) = interp1(T_interp(idx_night(idx)),tmp,T_interp);
            case 7; Qdayflog(:,it) = interp1(T_interp(idx_night(idx)),tmp,T_interp);
                
        end
        
        
        
    end
end






%% *********************************************************************************************************************************************************

% diurnal properties

nhoursday = 18;
nFS = floor(2*Fs);% new sampling rate, Ny*2 thus 30 min
doplot = 0;
nfig = 1;
for dwin = 1:1
    
    n = 1;
    ns = 1;
    spac_qg = 2 ; %multiple of 15 min
    T_day = T_interp(idx_days);
    clear Dayf_var & Pt_dayf & phiQ & phiP & validQ & validP & lagP & lagQ
    % Choose appropriate bandpass
    % ---------------------------------------------------------------------
    Qf = Qdayf(idx_days,dwin);
    Qflog = Qdayflog(idx_days,dwin);
    Pf = Pdayf(idx_days,dwin);
    Vf = Vdayf(idx_days,dwin);
    Pressf = Pressdayf(idx_days,dwin);
    Rf = Rdayf(idx_days,dwin);
    Sf = Sdayf(idx_days,dwin);
    
    % Define hydrological day
    % ---------------------------------------------------------------------
    
    mq = movmin(Qf,0.9,'SamplePoints',T_day); % search min(Q) over 22 hours
    idx = find(Qf == mq);
    % verify that each day is separated by more than n hours
    idx2 = find(diff(T_day(idx))>1/nhoursday)+1; % take the next one because diff shorts the vector by one
    idx = idx(idx2);
    
    
    % Loop through each day
    % ***********************************************************************************************************************
    for it = 1:length(idx)-1
        
        % Cut the time serie and define each day
        % ---------------------------------------------------------------------
        % idx(it) is one minimum, idx(it+1) is where the next minimum is
        Ttmp = T_day(idx(it):idx(it+1));                                    % take the corresponding date between the minimum
        Ttmpreg = (0:1/nFS:1-1/nFS)+floor(Ttmp(1));                         % regular grid from T(1) to T(end)
        Qftmp = Qf(idx(it):idx(it+1));                                      % correspoding Q
        Pftmp = Pf(idx(it):idx(it+1));                                      % corresponding Pw
        Vftmp = Vf(idx(it):idx(it+1));                                      % corresponding Vel
        Pressftmp = Pressf(idx(it):idx(it+1));                              % corresponding water pressure
        Rftmp = Rf(idx(it):idx(it+1));                                      % corresponding R
        Sftmp = Sf(idx(it):idx(it+1));                                      % corresponding S
        
        % interp the signals on the new regular grid,
        % ---------------------------------------------------------------------
        %  then each days starts at the same moment, the hydro day beginning here set to 0
        told = Ttmp-Ttmp(1)+floor(Ttmp(1));                                 % real time shifted to root(day) to make it begins a 0
        Qftmp = interp1(told,Qftmp,Ttmpreg);
        Pftmp = interp1(told,Pftmp,Ttmpreg);
        Vftmp = interp1(told,Vftmp,Ttmpreg);
        Pressftmp = interp1(told,Pressftmp,Ttmpreg);
        Rftmp = interp1(told,Rftmp,Ttmpreg);
        Sftmp = interp1(told,Sftmp,Ttmpreg);
        
        % quantify daily variability with the variance (std)^2
        % ---------------------------------------------------------------------
        Dayf_var(it,1) = std(Qftmp,'omitnan');
        Dayf_var(it,2) = std(Pftmp,'omitnan');
        Dayf_var(it,3) = std(Vftmp,'omitnan');
        Dayf_var(it,4) = std(Pressftmp,'omitnan');
        Dayf_var(it,5) = std(Rftmp,'omitnan');
        Dayf_var(it,6) = std(Sftmp,'omitnan');
        
        Dayf_var(it,1) = min(Qftmp)-max(Qftmp);
        Dayf_var(it,2) = min(Pftmp)-max(Pftmp);
        Dayf_var(it,3) = min(Vftmp)-max(Vftmp);
        Dayf_var(it,4) = min(Pressftmp)-max(Pressftmp);
        Dayf_var(it,5) = min(Rftmp)-max(Rftmp);
        Dayf_var(it,6) = min(Sftmp)-max(Sftmp);
        
        
        % properties of the hydrological day
        % ---------------------------------------------------------------------
        Pt_dayf(it,1) = (Ttmp(1)-floor(Ttmp(1)));           % day beginning
        Pt_dayf(it,2) = (Ttmp(end)-Ttmp(1));              % day real duration
        Pt_dayf(it,3) = median(Ttmpreg);                  % day
        
        % compute the rise and fall idenxing
        % ---------------------------------------------------------------------
        
        % compute with Q
        Qmin = [Qftmp(1) Qftmp(end)];                                       % minimum Q between rising and falling, to limit the bias between rise and fall
        gQ = gradient(Qftmp,spac_qg)./abs(gradient(Qftmp,spac_qg));         % compute the gradient, the value is 1 if rising, -1 if falling
        
        idx_risingQ = find(gQ == 1 & Qftmp >=max(Qmin));
        idx_fallingQ = find(gQ == -1& Qftmp >=max(Qmin));
        
        % compute with pressure
        Pressmin = [Pressftmp(1) Pressftmp(end)];                           % minimum press between rising and falling, to limit the bias between rise and fall
        gP = gradient(Pressftmp,spac_qg)./abs(gradient(Pressftmp,spac_qg)); % compute the gradient, the value is 1 if rising, -1 if falling
        
        idx_risingP = find(gP == 1 & Pressftmp >=max(Pressmin));
        idx_fallingP = find(gP == -1& Pressftmp >=max(Pressmin));
        
        % normalize variables by their minimum
        % ---------------------------------------------------------------------
        Qftmpnorm = Qftmp - min(Qftmp);
        Pftmpnorm = Pftmp - min(Pftmp);
        Vftmpnorm = Vftmp - min(Vftmp);
        Pressftmpnorm = Pressftmp - min(Pressftmp);
        Rftmpnorm = Rftmp - min(Rftmp);
        Sftmpnorm = Sftmp - min(Sftmp);
        
        % compute Q related hysteresis
        % ---------------------------------------------------------------------
        idx_falling = idx_fallingQ;
        idx_rising = idx_risingQ;
        
        % hysteresis
        phiQ(it,1) = NaN;
        phiQ(it,2) = ((mean(Pftmpnorm(idx_rising))-mean(Pftmpnorm(idx_falling)))/mean(Pftmpnorm(idx_falling)));
        phiQ(it,3) = ((mean(Vftmpnorm(idx_rising))-mean(Vftmpnorm(idx_falling)))/mean(Vftmpnorm(idx_falling)));
        phiQ(it,4) = ((mean(Pressftmpnorm(idx_rising))-mean(Pressftmpnorm(idx_falling)))/mean(Pressftmpnorm(idx_falling)));
        phiQ(it,5) = ((mean(Rftmpnorm(idx_rising))-mean(Rftmpnorm(idx_falling)))/mean(Rftmpnorm(idx_falling)));
        phiQ(it,6) = ((mean(Sftmpnorm(idx_rising))-mean(Sftmpnorm(idx_falling)))/mean(Sftmpnorm(idx_falling)));
        
        % valid if tmp(rising(1)<rising(end))
        if isfinite(idx_rising)
            validQ(it,1) =0;
            if Pftmp(1) < Pftmp(idx_rising(end)) validQ(it,2) = 1; else validQ(it,2) = 0; end
            if Vftmp(1) < Vftmp(idx_rising(end)) validQ(it,3) = 1; else validQ(it,3) = 0; end
            if Pressftmp(1) < Pressftmp(idx_rising(end)) validQ(it,4) = 1; else validQ(it,4) = 0; end
            if Rftmp(1) < Rftmp(idx_rising(end)) validQ(it,5) = 1; else validQ(it,5) = 0; end
            if Sftmp(1) < Sftmp(idx_rising(end)) validQ(it,6) = 1; else validQ(it,6) = 0; end
            
            % lag between peaks
            idx1 = find(Qftmp==max(Qftmp));
            if isfinite(idx1)
                lagQ(it,1) = NaN;
                lagQ(it,2) = Ttmpreg(idx1(1)) - mean(Ttmpreg(Pftmp==max(Pftmp(idx_rising))),'omitnan');
                lagQ(it,3) = Ttmpreg(idx1(1)) - mean(Ttmpreg(Vftmp==max(Vftmp(idx_rising))),'omitnan');
                lagQ(it,4) = Ttmpreg(idx1(1)) - mean(Ttmpreg(Pressftmp==max(Pressftmp(idx_rising))),'omitnan');
                lagQ(it,5) = Ttmpreg(idx1(1)) - mean(Ttmpreg(Rftmp==max(Rftmp(idx_rising))),'omitnan');
                lagQ(it,6) = Ttmpreg(idx1(1)) - mean(Ttmpreg(Sftmp==max(Sftmp(idx_rising))),'omitnan');
                
                %idx_reg = max(1,idx1(1)-length(Ttmpreg)/4):min(idx1(1)+length(Ttmpreg)/4,length(Ttmpreg));
                %idx_reg = max(1,idx1(1)-2*length(Ttmpreg)/4):min(idx1(1)+2*length(Ttmpreg)/4,length(Ttmpreg));
                idx_reg = 1:length(Ttmpreg)/4*3; n =5;
                % modified 08/11/19 to capture al maximums
                gP = abs(gradient(Pftmp(idx_reg),spac_qg)); [gP, index] = sort(gP);
                gP = index(1:n);
                gR = abs(gradient(Rftmp(idx_reg),spac_qg)); [gR, index] = sort(gR);
                gR = index(1:n);
                gS = abs(gradient(Sftmp(idx_reg),spac_qg)); [gS, index] = sort(gS);
                gS = index(1:n);
                
                lagQ(it,2) = Ttmpreg(idx1(1)) - mean(Ttmpreg(Pftmp==max(Pftmp(idx_reg(gP)))),'omitnan');
                lagQ(it,3) = Ttmpreg(idx1(1)) - mean(Ttmpreg(Vftmp==max(Vftmp(idx_reg))),'omitnan');
                lagQ(it,4) = Ttmpreg(idx1(1)) - mean(Ttmpreg(Pressftmp==max(Pressftmp(idx_reg))),'omitnan');
                lagQ(it,5) = Ttmpreg(idx1(1)) - mean(Ttmpreg(Rftmp==max(Rftmp(idx_reg(gR)))),'omitnan');
                lagQ(it,6) = Ttmpreg(idx1(1)) - mean(Ttmpreg(Sftmp==max(Sftmp(idx_reg(gS)))),'omitnan');
            end
            
        else
            validQ(it,:) = 0;
            lagQ(it,:) = NaN;
        end
        
        % compute Press related hysteresis
        % ---------------------------------------------------------------------
        idx_falling = idx_fallingP;
        idx_rising = idx_risingP;
        
        phiP(it,1) = ((mean(Qftmpnorm(idx_rising))-mean(Qftmpnorm(idx_falling)))/mean(Qftmpnorm(idx_falling)));
        phiP(it,2) = ((mean(Pftmpnorm(idx_rising))-mean(Pftmpnorm(idx_falling)))/mean(Pftmpnorm(idx_falling)));
        phiP(it,3) = ((mean(Vftmpnorm(idx_rising))-mean(Vftmpnorm(idx_falling)))/mean(Vftmpnorm(idx_falling)));
        phiP(it,4) = NaN;
        phiP(it,5) = ((mean(Rftmpnorm(idx_rising))-mean(Rftmpnorm(idx_falling)))/mean(Rftmpnorm(idx_falling)));
        phiP(it,6) = ((mean(Sftmpnorm(idx_rising))-mean(Sftmpnorm(idx_falling)))/mean(Sftmpnorm(idx_falling)));
        
        % valid if tmp(rising(1)<rising(end))
        if isfinite(idx_rising)
            if Qftmp(1) < Qftmp(idx_rising(end)) validP(it,1) = 1; else validP(it,1) = 0; end
            if Pftmp(1) < Pftmp(idx_rising(end)) validP(it,2) = 1; else validP(it,2) = 0; end
            if Vftmp(1) < Vftmp(idx_rising(end)) validP(it,3) = 1; else validP(it,3) = 0; end
            validP(it,4) = 0;
            if Rftmp(1) < Rftmp(idx_rising(end)) validP(it,5) = 1; else validP(it,5) = 0; end
            if Sftmp(1) < Sftmp(idx_rising(end)) validP(it,6) = 1; else validP(it,6) = 0; end
            
            % lag between peaks
            idx1 = find(Pressftmp==max(Pressftmp));
            if isfinite(idx1)
                lagP(it,1) = Ttmpreg(idx1(1)) - mean(Ttmpreg(Qftmp==max(Qftmp(idx_rising))),'omitnan');
                lagP(it,2) = Ttmpreg(idx1(1)) - mean(Ttmpreg(Pftmp==max(Pftmp(idx_rising))),'omitnan');
                lagP(it,3) = Ttmpreg(idx1(1)) - mean(Ttmpreg(Vftmp==max(Vftmp(idx_rising))),'omitnan');
                lagP(it,4) = NaN;
                lagP(it,5) = Ttmpreg(idx1(1)) - mean(Ttmpreg(Rftmp==max(Rftmp(idx_rising))),'omitnan');
                lagP(it,6) = Ttmpreg(idx1(1)) - mean(Ttmpreg(Sftmp==max(Sftmp(idx_rising))),'omitnan');
                
            end
        else
            validP(it,:) = 0;
            lagP(it,:) = NaN;
        end
        % Plot data
        % ***********************************************************************************************************************
        %         figure(f1)
        %
        %         if ns > 30
        %
        %             print(f1,['Ndaily ' num2str(1/Fcuthigh_bp(dwin)*24) '-' num2str(1/(Fcutlow_bp(dwin))*24) 'h ' num2str(nfig)],'-dpng')
        %             close(f1)
        %             f1 = figure('units','normalized','outerposition',[0 0 1 1],'DefaultAxesFontSize',10);
        %             ns=1;
        %             subplot(6,5,ns)
        %             nfig = nfig +1;
        %
        %         else
        %             subplot(6,5,ns)
        %         end
        %         %
        %         subplot(6,5,ns)
        %         qmin = Qftmp-min(Qftmp);
        %         plot(Qftmp,'b')
        %         hold on
        %         yyaxis right
        %         plot(Pftmp,'r')
        %         ylabel('')
        %         xlabel('')
        %         title(['d' num2str(Ttmpreg(1)) ' it' num2str(it) ' \phi' num2str(floor(phiQ(it,2)*10)/10) ' lag' num2str(floor(lagQ(it,2)*10)/10)])
        %
        %         ns = ns+1;
        %         subplot(6,5,ns)
        %         scatter(qmin,Pftmp-min(Pftmp),10,Ttmpreg,'filled')
        %         ylabel('')
        %         xlabel('')
        %         ylim([0 max(Pftmp-min(Pftmp))])
        %         ns = ns+1;
        
        % ----------
        
        
    end
    
    % Clean time series
    % ***********************************************************************************************************************
    [Tregc,ia,ic] = unique(Pt_dayf(:,3),'sorted');
    idx = idx(ia);
    lagP = lagP(ia,:);
    lagQ = lagQ(ia,:);
    Dayf_var = Dayf_var(ia,:);
    phiQ = phiQ(ia,:);
    phiP = phiP(ia,:);
    validQ = validQ(ia,:);
    validP = validP(ia,:);
    Pt_dayf = Pt_dayf(ia,:);
    
end

% *********************************************************************************************************************************************************
% *********************************************************************************************************************************************************








