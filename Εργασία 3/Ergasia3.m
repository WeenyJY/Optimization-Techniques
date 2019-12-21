% Konstantinos Letros 8851
% Optimization Techniques
% Project 3 -
%


%% Clean the screen

clc
clear
close all;
format long;

%% Parameters
e = 1e-4;
plotNum = 0;

% Inital Conditions
x_1 = 5;
y_1 = -10;

% Boundaries
a1 = 3;
b1 = 30;
a2 = -25;
b2 = -5;

%% Steepest Descent
fprintf("####### STEEPEST DESCENT #######\n\n")

% Count Number of Plots
plotNum = plotFunction(plotNum,a1,b1,a2,b2);

fprintf("Initial Conditions [x,y] =  [%f,%f]\n",x_1,y_1);

% Steepest Descent
[x,y,k] = steepestDescent(x_1,y_1,e,a1,b1,a2,b2);

if k < 101
    fprintf("Min(f) = %f at [x,y] = [%f,%f] after %d repetitions\n\n", ...
        f(x(end),y(end)),x(end),y(end),k);
end

% Plot Trace
tracePlot (x,y,k,plotNum-1)
tracePlot (x,y,k,plotNum)


figure(plotNum-1)
title('3D Plot - Steepest Descent - Optimized Gamma')
xlabel("x")
ylabel("y")
zlabel("f(x,y)")

figure(plotNum)
title('2D Plot - Steepest Descent - Optimized Gamma')
xlabel("x")
ylabel("y")
zlabel("f(x,y)")

%% Save Plots

% for i = 1 : plotNum
%     figure(i)
%       if (mod(i,2)==1)
%          view(-5,-20)
%       end
%     savePlot([mfilename,'_',num2str(i)])
% end

%% Functions

% Objective Function
function res = f(x,y)

res = x.*y+2*(x-y).^2;

end

% Beta function
function res = Beta(x,y,a1,b1,a2,b2)

h{1} = @(x,y) a1-x;
h{2} = @(x,y) x-b1;
h{3} = @(x,y) a2-y;
h{4} = @(x,y) y-b2;

F = @(x,y,h) -1/h; 

res = 0;
for i = 1:4
    res = res + F(x,y,h{i}(x,y));
end

end

% Auxiliary Function 
function res = phi(x,y,a1,a2,b1,b2)

r_k = 0.5;
res = f(x,y)+r_k*Beta(x,y,a1,b1,a2,b2);

end

% Derivative of the Auxiliary Function
function res = gradphi(x,y)

res = 

end

% Plot given Function
function plotNum = plotFunction(plotNum,a1,b1,a2,b2)

x =a1-3:0.1:b1+3;
y =a2-3:0.1:b2+3;

[X,Y] = meshgrid(x,y);

func = [];
for i = 1:length(y)
    for j = 1:length(x)
        func(i,j) = f(X(i,j),Y(i,j));
    end
end

figure(plotNum + 1)
surf(X,Y,func)
view(-10,25)
colorbar

t = -1e+3:10:1e+4;
line(a1*ones(size(t)),a2*ones(size(t)),t,'linewidth',3,'Color','#EDB120','LineStyle','--');
line(a1*ones(size(t)),b1*ones(size(t)),t,'linewidth',3,'Color','#EDB120','LineStyle','--');
line(b1*ones(size(t)),b2*ones(size(t)),t,'linewidth',3,'Color','#EDB120','LineStyle','--');
line(b1*ones(size(t)),a2*ones(size(t)),t,'linewidth',3,'Color','#EDB120','LineStyle','--');

figure(plotNum + 2)
contour(X,Y,func,20)
colorbar
rectangle('Position',[a1 a2 b1-a1 b2-a2],'linewidth',3,'EdgeColor','#EDB120','LineStyle','--');

plotNum = plotNum + 3;

end

% Steepest Descent Method
function [x,y,k] = steepestDescent(x,y,e,a1,b1,a2,b2)
k = 1;
d = [];
gamma = [];

while( norm( gradphi(x(k),y(k)) ) >= e)
    
    d(:,k) = - gradphi(x(k),y(k));
    gamma(k) = calcGamma(x(k),y(k),d(:,k),a1,b1,a2,b2);
    x(k+1) = x(k) + gamma(k)*d(1,k);
    y(k+1) = y(k) + gamma(k)*d(2,k);
    
    k = k + 1;
    if k>50000
        fprintf("Algorithm did not converge\n")
        break
    end
end

end
function gamma = calcGamma(x_k, y_k, d_k,a1,b1,a2,b2)

% function to be minimized with respect to gamma
func = @(gamma) phi(x_k+gamma*d_k(1) , y_k+gamma*d_k(2) ,a1,b1,a2,b2);

% Golden Section Parameters
l = 1e-3;
a = 1e-4;
b = 2;

% Minimization using Golden Section Method
gamma = goldenSectionMethod(func,l,a,b);

end

function minX = goldenSectionMethod(func,l,a,b)

gamma = (-1+sqrt(5))/2;
k = 1;
x1(k)=a(k)+(1-gamma)*(b(k)-a(k));
x2(k)=a(k)+gamma*(b(k)-a(k));

while b(k)-a(k)>=l
    
    if(func(x1(k))<=func(x2(k)))
        a(k+1)= a(k);
        b(k+1) = x2(k);
        x1(k+1) = a(k+1) + (1-gamma)*(b(k+1)-a(k+1));
        x2(k+1) = x1(k);
    else
        a(k+1) = x1(k);
        b(k+1)= b(k);
        x1(k+1) = x2(k);
        x2(k+1) = a(k+1) + gamma*(b(k+1)-a(k+1));
    end
    k = k + 1;
end

% Minimum point in the center of the final interval
minX = (a(end)+b(end))/2;

end

% Plot trace
function tracePlot (x,y,k,plotNum)

trace_f = [];

for i = 1:k
    trace_f(i) = f(x(i),y(i));
end

figure(plotNum)
hold on;
plot3(x,y,trace_f,"-r+",'linewidth',1)

% Plot minimum point
hold on;
plot3(x(end),y(end),f(x(end),y(end)),"-r*",'linewidth',7)

end

% Function to automatically save plots in high resolution
function savePlot(name)

% Resize current figure to fullscreen for higher resolution image
set(gcf, 'Position', get(0, 'Screensize'));

% Save current figure with the specified name
saveas(gcf, join([name,'.jpg']));

% Resize current figure back to normal
set(gcf,'position',get(0,'defaultfigureposition'));

end