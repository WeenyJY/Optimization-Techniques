clc
clear
close all;
format long;

%% Parameters

plotNum = 0;

% Inital Conditions
x_1 = [8,-5, 11];
y_1 = [3, 7,  3];

% Gamma Parameters
gamma_A = [0.1,1,2,10];
gamma_B = [0.1,0.3,0.01];

% Accuracy Parameters
e_A = [0.01,0.01,0.01,0.01];
e_B = [0.01,0.02,0.01];

% Projection Parameters
s = [15,1,0.1];
a1 = -20;
b1 =  10;
a2 = -12;
b2 =  15;

%% Steepest Descent - Part A (a)
fprintf("####### STEEPEST DESCENT #######\n\n")
partB = 0;

for i = 1 : 4
    
    % Count Number of Plots
    plotNum = plotFunction(plotNum,partB,a1,b1,a2,b2);
    
    % For all Initial Conditions
    for i=1:length(x_1)
        
        fprintf("Initial Conditions [x,y] =  [%f,%f]\n",x_1(i),y_1(i));
        
        fprintf("Gamma = %f , Accuracy e = %f \n", gamma_A(i) ,  e_A(i) )
        
        % Steepest Descent
        [x,y,k] = steepestDescent(x_1(i),y_1(i),e_A(i),gamma_A(i));
        
        if k < 101
            fprintf("Min(f) = %f at [x,y] = [%f,%f] after %d repetitions\n\n", ...
                f(x(end),y(end)),x(end),y(end),k);
        end
        
        % Plot Trace
        tracePlot (x,y,k,plotNum-1)
        tracePlot (x,y,k,plotNum)
    end
    
    figure(plotNum-1)
    title(['3D Plot - Steepest Descent - Gamma = ',num2str( gamma_A(i) ),' - Accuracy e = ',num2str(e_A(i))])
    xlabel("x")
    ylabel("y")
    zlabel("f(x,y)")
    
    figure(plotNum)
    title(['2D Plot - Steepest Descent - Gamma = ',num2str( gamma_A(i) ),' - Accuracy e = ',num2str(e_A(i))])
    xlabel("x")
    ylabel("y")
    zlabel("f(x,y)")
end

%% Steepest Descent with Projection - Part B (b,c,d)
fprintf("\n\n####### STEEPEST DESCENT with Projection #######\n\n")
partB = 1;

for i = 1 : 3
    
    % Count Number of Plots
    plotNum = plotFunction(plotNum,partB,a1,b1,a2,b2);
    
    fprintf("Initial Conditions [x,y] =  [%f,%f]\n",x_1(i),y_1(i));
    
    fprintf("Gamma = %f , Accuracy e = %f , s = %f \n", gamma_B(i),  e_B(i), s(i) )
    
    % Steepest Descent
    [x,y,k] = steepestDescentProj(x_1(i),y_1(i),e_B(i),gamma_B(i),s(i),a1,b1,a2,b2);
    
    if k < 101
        fprintf("Min(f) = %f at [x,y] = [%f,%f] after %d repetitions\n\n", ...
            f(x(end),y(end)),x(end),y(end),k);
    end
    
    % Plot Trace
    tracePlot (x,y,k,plotNum-1)
    tracePlot (x,y,k,plotNum)
    
    
    figure(plotNum-1)
    title(['3D Plot - Steepest Descent with Projection - Gamma = ',num2str( gamma_B(i) ),...
        ' - Accuracy e = ',num2str(e_B(i)),' - s_k = ',num2str(s(i))])
    
    xlabel("x")
    ylabel("y")
    zlabel("f(x,y)")
    
    figure(plotNum)
    title(['2D Plot - Steepest Descent with Projection - Gamma = ',num2str( gamma_B(i) ),...
        ' - Accuracy e = ',num2str(e_B(i)),' - s_k = ',num2str(s(i))])
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
function plotNum = plotFunction(plotNum,partB,a1,b1,a2,b2)


x = linspace(-21, 21, 100);
y = linspace(-21, 21, 100);

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
if partB == 1
     t = linspace(-50,500,1000);
     line(a1*ones(size(t)),a2*ones(size(t)),t,'linewidth',3,'Color','#EDB120','LineStyle','--');
     line(a1*ones(size(t)),b1*ones(size(t)),t,'linewidth',3,'Color','#EDB120','LineStyle','--');
     line(b1*ones(size(t)),b2*ones(size(t)),t,'linewidth',3,'Color','#EDB120','LineStyle','--');
     line(b1*ones(size(t)),a2*ones(size(t)),t,'linewidth',3,'Color','#EDB120','LineStyle','--');
end

figure(plotNum + 2)
contour(X,Y,func,20)
colorbar
if partB == 1
    rectangle('Position',[a1 a2 b1-a1 b2-a2],'linewidth',3,'EdgeColor','#EDB120','LineStyle','--');
end

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
    
    if(k>5000)
        fprintf("Diverging...\n\n");
        break;
    end
    
end

end


function [x,y,k] = steepestDescentProj(x,y,e,gammaMethod,s,a1,b1,a2,b2)
k = 1;
x_bar = [];
y_bar = [];
gamma = [];

while( norm( gradf(x(k),y(k)) ) >= e)
    
    gamma(k) = gammaMethod;
    
    x_bar(k) = calcProj(s,x(k),y(k),a1,b1,1);
    y_bar(k) = calcProj(s,x(k),y(k),a2,b2,2);
    
    x(k+1) = x(k) + gamma(k)*(x_bar(k)-x(k));
    y(k+1) = y(k) + gamma(k)*(y_bar(k)-y(k));
    
    k = k + 1;
    
      if k>5000
        fprintf("Diverging...\n\n");
        break;
      end
    
end

end

function proj = calcProj(s,x_k,y_k,a,b,index)

if index == 1
    vec_k = x_k;
elseif index == 2
    vec_k = y_k;
else
    fprintf("ERROR")
end

if vec_k <= a
    proj = a;
elseif vec_k >= b
    proj = b;
else
    grad = gradf(x_k,y_k);
    proj = vec_k-s*grad(index);
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