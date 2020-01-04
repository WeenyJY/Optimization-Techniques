% Konstantinos Letros 8851
% Optimization Techniques
% Project - Genetic Algorithm

%% Clean the screen

clc
clear
close all;
format long;

%% Problem Definition
% Rate of Incoming Vehicles
V = 100;

% Road Capacity c_i (Row Vector)
c = [59.85, 43.05, 53.55, 26.25, ...
    44.10, 64.05, 36.75, 35.70, ...
    13.65, 21.00, 54.60, 63.00, ...
    33.60, 54.60, 46.20, 31.50];

% Constant a_i
a = ones(1,16);

% Minimum time t_i
t = c/4;

% Initial Objective Function
f = @(x) sum(x.*t + a.*(x.^2)./(1-x./c));

% Equality Constraints
h = @(x) [ ...
    % 9 Equalities
    x(1)+x(2)+x(3)-V;
    x(6)+x(7)-x(1);
    x(8)+x(9)-x(7);
    x(15)-x(9)-x(10);
    x(5)+x(6)+x(8)-x(10)-x(11)-x(16);
    x(3)+x(4)-x(5)-x(12);
    x(2)-x(4)-x(13);
    x(11)+x(12)+x(13)-x(14);
    x(14)+x(15)+x(16)-V];

Aeq = [1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
    -1 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0  -1 1 1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0  0 -1 -1 0 0 0 0 1 0;
    0 0 0 0 1 1 0  1 0 -1 -1 0 0 0 0 -1;
    0 0 1 1 -1 0 0 0 0 0 0 -1 0 0 0 0;
    0 1 0 -1 0 0 0 0 0 0 0 0 -1 0 0 0;
    0 0 0 0 0 0 0 0 0 0 1 1 1 -1 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1];

beq = [V;zeros(7,1);V];

% Inequality Constraints
lb = zeros(1,16); % x >= 0
ub = c; % x <= c , lambda < Infinity

r = 1e3;

% Augmented Lagrangian
func = @(params) f(params(1:16));
% func = @(params) f(params(1:16)) + params(17:25)*h(params(1:16))+r*h(params(1:16))'*h(params(1:16));
% func = @(params) r/(f(params(1:16)) + r*sum(abs(h(params(1:16)))));

%% Testing

options = optimoptions('ga','FunctionTolerance',0,'PlotFcn', @gaplotbestf,'MigrationFraction',0.8);

[x,fval,exitflag,output,population,scores]  = ga(func,16,[],[],Aeq,beq,lb,ub,[],options);
f(x(1:16))