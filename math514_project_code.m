clc; clear; close all; format long

N=1e6; S0 = N-1; I0 = 1;
beta = 0.3; gamma = 0.8;

% y = [S I R]'
dS = @(t,y) -beta * y(1) * y(2) / N;
dI = @(t,y) beta * y(1) * y(2) / N - gamma * y(2);
dR = @(t,y) gamma * y(2);
f = @(t,y) [dS(t,y); dI(t,y); dR(t,y)];
y0 = [S0 I0 0]'; Nsteps = 5000; 
a = 0; b = 1; h=(b-a)/Nsteps;

%% RK4 Loop
t = zeros(1, Nsteps+1);
y = zeros(3, Nsteps+1);
t(1) = a; y(:, 1) = y0;
for i = 1:N
    k_1 = f(t(i),y(:,i));
    k_2 = f(t(i)+0.5*h,y(:,i)+0.5*h*k_1);
    k_3 = f((t(i)+0.5*h),(y(:,i)+0.5*h*k_2));
    k_4 = f((t(i)+h),(y(:,i)+k_3*h));
    y(:,i+1) = y(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    t(i+1) = t(i) + h;
end

plot(t,y,'LineWidth',1.5)
legend('S', 'I', 'R')
xlabel('Time t');
ylabel('Population');
title('SIR Model Approximation');
grid on;

%% Stiffness

eig_diff = zeros(size(t));

for k = 1:length(t)
    J = [ -(beta*y(2,k))/N, -(beta*y(1,k))/N,      0;
           (beta*y(2,k))/N, (beta*y(1,k))/N-gamma, 0;
            0,              gamma,                 0];
    lambda = eig(J);
    lambda = sort(lambda, 'ComparisonMethod', 'real'); 
    eig_diff(k) = lambda(2) - lambda(1);  
end

% Plot result
figure;
plot(t, eig_diff, 'LineWidth', 1.5);
xlabel('Time t');
ylabel('Eigenvalue difference |\lambda_2 - \lambda_1|');
title('Difference Between Nonzero Eigenvalues Over Time');
grid on;