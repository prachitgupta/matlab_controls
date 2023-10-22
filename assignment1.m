%%problem1
% Define the system parameters symbolically
A1 = [-1 0; 0 -2];
B1 = [1; sqrt(2)];
C1 = [1 -1*sqrt(2)/2];
D1 = 0;

% Define initial condition x0
x0_1 = [1; -1];
t_values = 0:0.01:10;

[G,h_t,step_t,zero_input_t, Adcf,x,tau,n] = calculate_ans(A1,B1,C1,D1,x0_1,t_values);

display_ans(G,h_t,step_t,zero_input_t, Adcf);

plot_ans(h_t, step_t, zero_input_t, x, t_values, tau,n); 

%%problem2
% Define the system parameters symbolically
A2 = [0 1 0; 0 0 1; -52 -30 -4];
B2 = [0;0;1];
C2 = [20 1 0];
D2 = 0;

% Define initial condition x0
x0_2 = [1; 2; 3];
t_values = 0:0.01:10;

[G,h_t,step_t,zero_input_t, Adcf,x,tau,n] = calculate_ans(A2,B2,C2,D2,x0_2,t_values);

display_ans(G,h_t,step_t,zero_input_t, Adcf);


plot_ans(h_t, step_t, zero_input_t, x, t_values, tau,n); 

%%problem3
% Define the system parameters symbolically
A3 = [0 1 0 0; 0 0 1 0; 0 0 0 1; -962 -126 -67 -4];
B3 = [0;0;0;1];
C3 = [300 0 0 0];
D3 = 0;

% Define initial condition x0
x0_3 = [4;3;2;1];
t_values = 0:0.01:10;

[G,h_t,step_t,zero_input_t, Adcf,x,tau,n] = calculate_ans(A3,B3,C3,D3,x0_3,t_values);

display_ans(G,h_t,step_t,zero_input_t, Adcf);

plot_ans(h_t, step_t, zero_input_t, x, t_values, tau,n); 

%%problem4
% Define the system parameters symbolically
A4 = [0 1 0 0; 0 0 1 0; 0 0 0 1; -680 -176 -86 -6];
B4 = [0;0;0;1];
C4 = [100 20 10 0];
D4 = 0;

% Define initial condition x0
x0_4 = [1; 2; 3; 4];
t_values = 0:0.01:10;

[G,h_t,step_t,zero_input_t, Adcf,x,tau,n] = calculate_ans(A4,B4,C4,D4,x0_4,t_values);

display_ans(G,h_t,step_t,zero_input_t, Adcf);

plot_ans(h_t, step_t, zero_input_t, x, t_values, tau,n); 


%% function definations
function [G,h_t,step_t,zero_input_t, Adcf,x, tau,n] = calculate_ans(A,B,C,D,x0,tspan)
syms s t; % Define symbolic variables
n = size(A,1);

%%Transfer functions
H = C * inv(s*eye(n) - A) * B;
G = H + D;

% Impulse response as a symbolic expression of t
h_t = ilaplace(H, s, t);

% Unit step response as a symbolic expression of t L(1) = 1/s
step_t = ilaplace(H / s, s, t);

% % Simulate zero input response
zero_input_response = ilaplace(C*inv(s*eye(n) - A)*x0,s,t);
% % Simplify the zero input response
zero_input_t = simplify(zero_input_response);

%%diagnol canonical form E = diagnol matrix consisting of eigenValues
[T, E] = eig(A);
Adcf = inv(T)*A*T;
Bdcf = inv(T)*B;
Cdcf = C*inv(T);
Ddcf = D;

%%Lyapunov analysis 

% Calculate the eigenvalues of the system matrix A
eigenvalues = eig(A);

% Check stability using eigenvalue analysis
if all(real(eigenvalues) < 0)
    disp('Eigenvalue Analysis: System is stable (all eigenvalues have negative real parts)');
else
    disp('Eigenvalue Analysis: System is unstable (at least one eigenvalue has a non-negative real part)');
end

% % Lyapunov stability analysis
% Q = -eye(2); % any Q
% P = lyap(A, Q);
% 
% % Check stability using Lyapunov analysis
% if all(real(eig(P)) > 0)
%     disp('Lyapunov Analysis: System is stable (all eigen values of P are strictly positive hence P is positive definate)');
% else
%     disp('Lyapunov Analysis: System is unstable (at least eigen value of P is non-positive)');
% end

%%trajectory
dxdt = @(T, x) A * x;

% Solve the system using ode45
[tau, x] = ode45(dxdt, tspan, x0);
end

function display_ans(G,h_t,step_t,zero_input_t, Adcf)
% Display symbolic expressions
disp("open loop transfer function");
disp(G);
disp('Impulse Response :');
disp(h_t);
disp('Unit Step Response :');
disp(step_t);
disp('Zero Input Response with Initial Condition :');
disp(zero_input_t);
disp('Diagnol Canonical Form:');
disp(Adcf);
end

function plot_ans(impulse_response, step_response,zero_input_response_simplified,x,t_values,T,n)
syms t; 
colorCodes = {'r', 'g', 'b', 'm'};
figure;
subplot(2, 2, 1);
% Plot impulse response (numerical)
hold on;
impulse_response_numerical = double(subs(impulse_response, t, t_values));
plot(t_values, impulse_response_numerical);
title('Impulse Response');
xlabel('Time');
ylabel('Amplitude');
hold on;
subplot(2, 2, 2);
% Plot unit step response (numerical)
step_response_numerical = double(subs(step_response, t, t_values));
plot(t_values, step_response_numerical);
title('Unit Step Response)');
xlabel('Time');
ylabel('Amplitude');
hold on;
subplot(2, 2, 3);
% % Plot zero input response (numerical)
zero_input_response_numerical = double(subs(zero_input_response_simplified, t, t_values));
plot(t_values, zero_input_response_numerical);
title('Zero Input Response');
xlabel('Time');
ylabel('Amplitude');
hold on;
subplot(2, 2, 4);
for c = 1:n
plot(T, x(:, c), colorCodes{c}, 'LineWidth', 2); % Plot x1 in blue
hold on;
end
hold off;
title('Trajectories');
xlabel('Time');
ylabel('Amplitude');
legend('x1', 'x2', 'x3 ', 'x4');
% Save the current figure as a PNG file
index = num2str(n);
formattedString = sprintf('Problem %.2f.png', index);
saveas(gcf, formattedString);
end