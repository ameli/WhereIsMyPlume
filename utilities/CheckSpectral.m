clear all

load InputLog.mat

U = Data.Ocean.AbsoluteVel_u;
V = Data.Ocean.AbsoluteVel_v;
T = Data.Time / (24*3600);    % In days
DepthIndex = 1;

StartTime = datenum([2018 08 16 05 00 00]);
EndTime = datenum([2018 08 17 11 00 00]);

StartTimeIndex = find(T >= StartTime,1,'first');
EndTimeIndex = find(EndTime <= T,1,'first');

u_raw = U(DepthIndex,StartTimeIndex:EndTimeIndex);
v_raw = V(DepthIndex,StartTimeIndex:EndTimeIndex);
t = T(StartTimeIndex:EndTimeIndex);

% Moving average
dt = t(2) - t(1);                         % one timestep of data (in unit of days)
MeanWindow = datenum([0 0 0 0 10 0]);     % 30 minutes duration (in unit of days)
MeanWindowLength = floor(MeanWindow / dt);
u_mean = movmean(u_raw,MeanWindowLength);
v_mean = movmean(v_raw,MeanWindowLength);
u = u_raw - u_mean;
v = v_raw - v_mean;
vel = [u(:),v(:)];

% Statistics of filtered velocities
Mean_u = mean(u(:));
Mean_v = mean(v(:));
Std_u = std(u(:));
Std_v = std(v(:));

fprintf('u: mean: %f, std: %f\n',Mean_u,Std_u);
fprintf('v: mean: %f, std: %f\n\n',Mean_v,Std_v);

% Ljung-Box Q-test
h_u = lbqtest(u);
h_v = lbqtest(v);
fprintf('Ljung-Box Q-test for u: %d\n',h_u)
fprintf('Ljung-Box Q-test for v: %d\n\n',h_v)

% Plot
EastColor = [0, 0.4470, 0.7410];
NorthColor = [0.8500, 0.3250, 0.0980];
NumberOfXTicks = 13;
figure()
subplot(2,2,1)
hold on
plot(t,u_raw,'DisplayName','raw','color',EastColor)
plot(t,u_mean,'DisplayName','mean','linewidth',1.5,'color','black')
xticks(linspace(t(1),t(end),NumberOfXTicks))
datetick('x','mm/dd HH:MM','keepticks','keeplimits')
xtickangle(30)
xlabel('Time (mm/dd HH:MM)');
ylabel('Velocity (m)');
xlim([t(1),t(end)]);
title('East Velocity')
grid()
legend()

subplot(2,2,3)
plot(t,u,'DisplayName','filtered','color',EastColor)
xticks(linspace(t(1),t(end),NumberOfXTicks))
datetick('x','mm/dd HH:MM','keepticks','keeplimits')
xtickangle(30)
xlabel('Time (mm/dd HH:MM)');
ylabel('Velocity (m)');
xlim([t(1),t(end)]);
title('East Velocity')
grid()
legend()

subplot(2,2,2)
hold on
plot(t,v_raw,'DisplayName','raw','color',NorthColor)
plot(t,v_mean,'DisplayName','mean','linewidth',1.5,'color','black')
xticks(linspace(t(1),t(end),NumberOfXTicks))
datetick('x','mm/dd HH:MM','keepticks','keeplimits')
xtickangle(30)
xlabel('Time (mm/dd HH:MM)');
ylabel('Velocity (m)');
xlim([t(1),t(end)]);
title('North Velocity')
grid()
legend()

subplot(2,2,4)
plot(t,v,'DisplayName','filtered','color',NorthColor)
xticks(linspace(t(1),t(end),NumberOfXTicks))
datetick('x','mm/dd HH:MM','keepticks','keeplimits')
xtickangle(30)
xlabel('Time (mm/dd HH:MM)');
ylabel('Velocity (m)');
xlim([t(1),t(end)]);
title('North Velocity')
grid()
legend()

% Probability density function
[f_u,xi_u] = ksdensity(u(:));
[f_v,xi_v] = ksdensity(v(:));

NormalDistribution= @(xi,mu,sigma) 1/sqrt(2 * 3.14159265 * sigma^2) * exp(-(xi-mu).^2/(2*sigma^2));
NormalDist_u = NormalDistribution(xi_u,Mean_u,Std_u);
NormalDist_v = NormalDistribution(xi_v,Mean_v,Std_v);

figure();
subplot(1,2,1)
plot(xi_u,f_u,'color','black')
hold on
plot(xi_u,NormalDist_u,'--','color','black')
xlim([xi_u(1),xi_u(end)])
xlabel('u (m/s)')
ylabel('f_u')
title('Probability density function of u')
grid()

subplot(1,2,2)
plot(xi_v,f_v,'color','black')
hold on
plot(xi_v,NormalDist_v,'--','color','black')
xlim([xi_v(1),xi_v(end)])
xlabel('v (m/s)')
ylabel('f_v')
title('Probability density function of v')
grid()

%% Spectrogram

% [p_u,w_u] = periodogram(u);
% [p_v,w_v] = periodogram(v);
N = 256;
window = ones(N,1)/N;
[p_uu,w_uu] = cpsd(u,u,window);
[p_vv,w_vv] = cpsd(v,v,window);
[p_uv,w_uv] = cpsd(u,v,window);
[p_vu,w_vu] = cpsd(v,u,window);

s(:,1,1) = p_uu;
s(:,1,2) = p_uv;
s(:,2,1) = p_vu;
s(:,2,2) = p_vv;

figure()
hold on
plot(w_uu,10*log10(s(:,1,1)),'DisplayName','s_{11}')
plot(w_uu,10*log10(s(:,1,2)),'DisplayName','s_{12}')
plot(w_uu,10*log10(s(:,2,1)),'DisplayName','s_{21}')
plot(w_uu,10*log10(s(:,2,2)),'DisplayName','s_{22}')
legend()
grid()
title('Spectogram')

%% Autocorrelation

[r,lags] = xcorr(vel);

ZeroLagIndex = find(lags == 0);
r = r(ZeroLagIndex:end,:);
lags = lags(ZeroLagIndex:end);

r = reshape(r,size(r,1),2,2);

r0 = squeeze(r(1,:,:));
c = inv(sqrtm(r0));
c2 = inv(r0);

% Autocorrelation
rho = zeros(size(r,1),2,2);
for i = 1:size(r,1)
    % rho(i,:,:) = c * squeeze(r(i,:,:)) * c;
    rho(i,:,:) = c2 * squeeze(r(i,:,:));
    % rho(i,:,:) = squeeze(r(i,:,:));
end

figure()
hold on
plot(lags,rho(:,1,1),'DisplayName','r_{11}')
plot(lags,rho(:,1,2),'DisplayName','r_{12}')
plot(lags,rho(:,2,1),'DisplayName','r_{21}')
plot(lags,rho(:,2,2),'DisplayName','r_{22}')
grid on
legend()
title('rho')

%% Log Autocorrelation

% LogRho = zeros(size(rho,1),2,2);
% for i = 1:size(LogRho,1)
%     LogRho(i,:,:) = logm(squeeze(rho(i,:,:))) / lags(i);
% end
% 
% figure()
% hold on
% plot(lags,LogRho(:,1,1),'DisplayName','r_{11}')
% plot(lags,LogRho(:,1,2),'DisplayName','r_{12}')
% plot(lags,LogRho(:,2,1),'DisplayName','r_{21}')
% plot(lags,LogRho(:,2,2),'DisplayName','r_{22}')
% grid on
% legend()
% title('Log rho')
