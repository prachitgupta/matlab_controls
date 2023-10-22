% Load model and parameters
open_system('ACC.mdl');

T = 0.5;
load('ACCparams.mat'); % Load ACCparams file with K_P, K_I, etc.
Fl = 0;
Kp = 2;
Ki = 1.5;
n = 0.9;
C = 1.5;

% % % Retrieve simulation results
t_sim = v.time;
v_sim = v.signals.values; 
e_sim =  e.signals.values; 
alpha_des_sim = alpha_des.signals.values;
alpha_ac_sim = alpha.signals.values; 
alpha_Sat = alpha_Sat.signals.values;

S = stepinfo(v_sim, t_sim,'SettlingTimeThreshold', 0.05);
disp(S);
% Calculate settling time and maximum overshoot (you need to define these metrics)
t_settling = S.SettlingTime;
% Plot results
figure;
subplot(2,2,1);
plot(t_sim, v_sim);
ylabel('velocity');
xlabel('Time (s)');
% % % Add grid to subplot 1 with grid size 0.02
%  grid on;
%  ax = gca;
%  ax.XTick = 0:0.02:max(t_sim);
%  ax.YTick = 15:0.02:max(v_sim);
subplot(2,2,2);
plot(t_sim, e_sim);
ylabel('error');
xlabel('Time (s)');

subplot(2,2,3);
plot(t_sim, alpha_des_sim);
ylabel('alphaDes');
xlabel('Time (s)');

subplot(2,2,4);
plot(t_sim, alpha_ac_sim);
ylabel('alphaAC');

xlabel('Time (s)');
