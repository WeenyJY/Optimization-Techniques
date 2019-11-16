clc
clear
close all;
format long;

%% Parameters


plotNum = 0;
gammaStrings = ["Constant Gamma","Optimized Gamma","Changing Gamma"];


% Inital Conditions
x_1 = [-0.5,0.5,-2,-3, -0.1];
y_1 = [-0.5,0.5,-1, 3, -0.9];

% Gamma Parameters
gamma = [ 0.1,1,2,10,...
        0.1,0.3,0.01 ];

% Accuracy Parameters
e = [0.01,0.01,0.01,0.01, 0.01,0.02,0.01];

%% Steepest Descent - Part A
fprintf("####### STEEPEST DESCENT #######\n\n")

% 3 Different Methods of Calculating Gamma
for method = 1 : 4
    
    % Count Number of Plots
    plotNum = plotFunction(plotNum);
    fprintf("******* Gamma Method %d *******\n\n",method)
    
    % For all Initial Conditions
    for i=1:length(x_1)
        
        fprintf("Initial Conditions [x,y] =  [%f,%f]\n",x_1(i),y_1(i));
        
        fprintf("Gamma = %f , Accuracy e = %f \n", gamma(method) ,  e(method) )
        
        % Steepest Descent
        [x,y,k] = steepestDescent(x_1(i),y_1(i),e(method),gamma(method));
        
        fprintf("Min(f) = %f at [x,y] = [%f,%f] after %d repetitions\n\n", ...
            f(x(end),y(end)),x(end),y(end),k);
        
        % Plot Trace
        tracePlot (x,y,k,plotNum-1)
        tracePlot (x,y,k,plotNum)
    end
    
    figure(plotNum-1)
    title(['3D Plot - Steepest Descent - Gamma = ',num2str( gamma(method) ),' - Accuracy e = ',num2str(e(method))])
    xlabel("x")
    ylabel("y")
    zlabel("f(x,y)")
    
    figure(plotNum)
    title(['2D Plot - Steepest Descent - Gamma = ',num2str( gamma(method) ),' - Accuracy e = ',num2str(e(method))])
    xlabel("x")
    ylabel("y")
    zlabel("f(x,y)")
end


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

res = 1/2*(x.^2+y.^2);

end

% Derivative of the Objective Function

function res = gradf(x,y)

res = [x ; y];

end

% Plot given Function
function plotNum = plotFunction(plotNum)


x = linspace(-3.5, 3.5, 100);
y = linspace(-3.5, 3.5, 100);

[X,Y] = meshgrid(x,y);

func = [];
for i = 1:length(x)
    for j = 1:length(y)
        func(i,j) = f(X(i,j),Y(i,j));
    end
end



figure(plotNum+1)
surf(X,Y,func)
view(-10,25)
colorbar

figure(plotNum + 2)
contour(X,Y,func,20)
colorbar

plotNum = plotNum + 2;

end

function [x,y,k] = steepestDescent(x,y,e,gammaMethod)
k = 1;
d = [];
gamma = [];


while( norm( gradf(x(k),y(k)) ) >= e)
    
    d(:,k) = - gradf(x(k),y(k));
    gamma(k) = gammaMethod;
    x(k+1) = x(k) + gamma(k)*d(1,k);
    y(k+1) = y(k) + gamma(k)*d(2,k);
    
    k = k + 1;
    
    if(k>100)
       fprintf("Diverging...");
       break;
    end
    
end

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