% Konstantinos Letros 8851
% Optimization Techniques
% Project 3 - Part A
% Objective Function g(x,y)

%% Clean the screen

clc
clear
close all;
format long;

%% Parameters
e = 1e-3;
plotNum = 0;

% Inital Conditions
x_1 = [-20,-5,-29,-3];
y_1 = [-3,-27,-1.5,-5.5];

% Boundaries
global b1 b2

b1 = -1;
b2 = -1;

%% Steepest Descent - f(x,y)
fprintf("####### STEEPEST DESCENT - f(x,y) #######\n\n")

% Count Number of Plots
plotNum = plotFunction(plotNum);
% For all Initial Conditions
for i=1:length(x_1)
        
fprintf("Initial Conditions [x,y] =  [%f,%f]\n",x_1(i),y_1(i));

% Steepest Descent
[x,y,k] = steepestDescent(x_1(i),y_1(i),e);

if k < 50001
    fprintf("Min(f) = %f at [x,y] = [%f,%f] after %d repetitions\n\n", ...
        g(x(end),y(end)),x(end),y(end),k);
end

% Plot Trace
color = rand(1,3);
tracePlot (x,y,k,plotNum-1,color)
tracePlot (x,y,k,plotNum,color)
end

figure(plotNum-1)
title('3D Plot - Steepest Descent - Optimized Gamma')
xlabel("x")
ylabel("y")
zlabel("g(x,y)")

figure(plotNum)
title('2D Plot - Steepest Descent - Optimized Gamma')
xlabel("x")
ylabel("y")
zlabel("g(x,y)")

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
function res = g(x,y)

res = (x-y).^2;

end

function res = gradg(x,y)

res = [2.*x-2.*y ; 2.*y-2.*x];

end


% Auxiliary Function and the derivative
function res = phiFun(n)
global b1 b2

h{1} = @(x,y) (x-b1);
h{2} = @(x,y) (y-b2);

gradh{1} = @(x,y) [1;0];
gradh{2} = @(x,y) [0;1];

F = @(h) (-1/h ).* (h < 0) + 1e+25.*(h>=0) ; 
gradF = @(h) (1./h.^2).* (h < 0)  + 1e+25.*(h>=0) ;

Beta = @(x,y) F(h{1}(x,y))+F(h{2}(x,y));
gradBeta = @(x,y) gradF(h{1}(x,y))*gradh{1}(x,y) +gradF(h{2}(x,y))*gradh{2}(x,y);

alpha = 2;

if n==1
    res = @(x,y,r) (gradg(x,y)+r*gradBeta(x,y));
else
    res = @(x,y) alpha * g(x,y)/Beta(x,y);
end

end


function res = gradphi(x,y,r)
func = phiFun(1);
res = func(x,y,r);
end

function res = r0(x,y)
func = phiFun(2);
res = func(x,y);
end

% Plot given Function
function plotNum = plotFunction(plotNum)
global b1 b2

x =-30:0.1:b1+3;
y =-30:0.1:b2+3;

[X,Y] = meshgrid(x,y);

func = [];
for i = 1:length(y)
    for j = 1:length(x)
        func(i,j) = g(X(i,j),Y(i,j));
    end
end

figure(plotNum + 1)
surf(X,Y,func)
view(-10,25)
colorbar

t = -1e+3:10:1e+4;
line(b1*ones(size(t)),b2*ones(size(t)),t,'linewidth',3,'Color','#EDB120','LineStyle','--');

figure(plotNum+2)
contour(X,Y,func,20)
colorbar
hold on
plot([b1,-30],[b1,b1],'linewidth',3,'color','#EDB120','LineStyle','--')
hold on
plot([b1,b1],[-30,b1],'linewidth',3,'color','#EDB120','LineStyle','--')
hold on
plot(-30:0.01:-1,-30:0.01:-1,'color','r','LineStyle','--')

plotNum = plotNum + 2;

end

% Steepest Descent Method
function [x,y,k] = steepestDescent(x,y,e)
k = 1;
d = [];
gamma = 0.005;
r = r0(x,y);

c = 0.9;

while( norm( gradphi(x(k),y(k),r(k)) ) >= e)
    
    f_prev=g(x(end),y(end));    
    
    d(:,k) = - gradphi(x(k),y(k),r(k));
    x(k+1) = x(k) + gamma*d(1,k);
    y(k+1) = y(k) + gamma*d(2,k);
    r(k+1) = c *r(k);
    
    k = k + 1;
    f_xk=g(x(end),y(end));
       
    if abs((f_xk-f_prev)/f_xk)<= e
        break
    elseif k>50000
        fprintf("Algorithm did not converge\n")
        break
    end
end

end

% Plot trace
function tracePlot (x,y,k,plotNum,num)

trace_g = [];

for i = 1:k
    trace_g(i) = g(x(i),y(i));
end

figure(plotNum)
hold on;
plot3(x,y,trace_g,'-+','color',num,'linewidth',1)

% Plot minimum point
hold on;
plot3(x(end),y(end),g(x(end),y(end)),"-r*",'linewidth',4)

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