% Konstantinos Letros 8851
% Optimization Techniques
% Project 2 - Modified Gradient Descent Methods - Part B
% Conjugate Gradients - Quasi Newton

%% Clean the screen

clc
clear
close all;
format long;

%% Parameters

e = 1e-4;
plotNum = 0;
gammaStrings = ["Constant Gamma","Optimized Gamma","Armijo's Gamma"];
dVecStrings = ["Conjugate Gradients Method","Quasi Newton Method"];

% Inital Conditions
x_1 = [0,-0.6,1];
y_1 = [0,-0.6,1];

%% Modified Gradient Descent
fprintf("####### MODIFIED GRADIENT DESCENT METHODS #######\n\n")


% 2 Different Methods of Calculating d_k Vector
for dVecMethod = 1 : 2
    
    fprintf("******* %s *******\n\n",dVecStrings(dVecMethod))
    
    % 3 Different Methods of Calculating Gamma
    for gammaMethod = 1 : 3
        
        % Count Number of Plots
        plotNum = plotFunction(plotNum);
        fprintf("******* %s *******\n\n",gammaStrings(gammaMethod))
        
        % For all Initial Conditions
        for i=1:length(x_1)
            
            fprintf("Initial Conditions [x,y] =  [%f,%f]\n",x_1(i),y_1(i));
            
            if gammaMethod == 1
                fprintf("Gamma is constant \n")
            elseif gammaMethod == 2
                fprintf("Gamma is optimized using Golden Section Method \n")
            else
                fprintf("Gamma is optimized using Armijo's Condition \n")
            end
            
            % Gradient Descent
            [x,y,k] = gradientDescent(x_1(i),y_1(i),e,gammaMethod,dVecMethod);
            
            fprintf("Min(g) = %f at [x,y] = [%f,%f] after %d repetitions\n\n", ...
                g(x(end),y(end)),x(end),y(end),k);
            
            % Plot Trace
            tracePlot (x,y,k,plotNum-2)
            tracePlot (x,y,k,plotNum-1)
            
            % Plot Objective function through iterations
            objFuncReps(x,y,k,plotNum,i)
            
        end
        
        figure(plotNum-2)
        title(['3D Plot - ',num2str(gammaStrings(gammaMethod)),' - ',num2str(dVecStrings(dVecMethod))])
        xlabel("x")
        ylabel("y")
        zlabel("g(x,y)")
        
        figure(plotNum-1)
        title(['2D Plot - ',num2str(gammaStrings(gammaMethod)),' - ',num2str(dVecStrings(dVecMethod))])
        xlabel("x")
        ylabel("y")
        
        figure(plotNum)
        sgtitle(['Objective Function through Iterations - ',num2str(gammaStrings(gammaMethod)),' - ',num2str(dVecStrings(dVecMethod))])
        
    end
end

%% Save Plots

% for i = 1 : plotNum
%     figure(i)
%     if (mod(i,3)==1)
%         view(-20,65)
%     end
%     savePlot([mfilename,'_',num2str(i)])
% end

%% Functions

% Plot given Function
function plotNum = plotFunction(plotNum)


x = linspace(-1.1, 1.1, 200);
y = linspace(-1.1, 1.1, 200);

[X,Y] = meshgrid(x,y);

func = [];
for i = 1:length(x)
    for j = 1:length(y)
        func(i,j) = g(X(i,j),Y(i,j));
    end
end

figure(plotNum+1)
surf(X,Y,func)
view(-10,25)
colorbar

figure(plotNum + 2)
contour(X,Y,func,20)
colorbar

plotNum = plotNum + 3;

end

function [x,y,k] = gradientDescent(x,y,e,gammaMethod,dVectorMethod)
k = 1;
d = [];
gamma = [];


while( norm( gradg(x(k),y(k)) ) >= e)
    
    if dVectorMethod == 1
        if k ~= 1
            d(:,k) = conjGradients(x(k),y(k),x(k-1),y(k-1),d(:,k-1));
        else
            d(:,1) = - gradg(x(k),y(k));
        end
    else
        if k ~= 1
            [d(:,k),Delta_k] = quasiNewton(x(k),y(k),x(k-1),y(k-1),Delta_k);
        else
            Delta_k = 4*eye(2);
            d(:,1) = -Delta_k*gradg(x(k),y(k));
        end
    end
    
    gamma(k) = calcGamma(x(k),y(k),d(:,k),gammaMethod);
    x(k+1) = x(k) + gamma(k)*d(1,k);
    y(k+1) = y(k) + gamma(k)*d(2,k);
    
    k = k + 1;
end

end

function d_k = conjGradients(x_k,y_k,x_prev,y_prev,d_prev)

beta_k = gradg(x_k,y_k)'*( gradg(x_k,y_k)-gradg(x_prev,y_prev) )/(gradg(x_prev,y_prev)'*gradg(x_prev,y_prev));
d_k = -gradg(x_k,y_k) + beta_k*d_prev;

end

function [d_k,Delta_k] = quasiNewton(x_k,y_k,x_next,y_next,Delta_prev)
xi = 0.5;


q_k = gradg(x_next,y_next)-gradg(x_k,y_k);
p_k = [x_next-x_k;y_next-y_k];

t_k = q_k'*Delta_prev*q_k;
v_k = p_k/(p_k' * q_k) - Delta_prev*q_k/t_k;
Delta_k = Delta_prev + p_k*p_k'/(p_k' * q_k) - ...
    Delta_prev*q_k*q_k'*Delta_prev/(q_k'*Delta_prev*q_k) + xi*t_k*v_k*v_k';

d_k = -Delta_k*gradg(x_k,y_k);
end

function gamma = calcGamma(x_k, y_k, d_k, gammaMethod)

if gammaMethod == 1
    gamma = 0.04;
elseif gammaMethod == 2
    
    % function to be minimized with respect to gamma
    func = @(gamma) g(x_k+gamma*d_k(1) , y_k+gamma*d_k(2) );
    
    % Golden Section Parameters
    l = 1e-3;
    a = 1e-4;
    b = 2;
    
    % Minimization using Golden Section Method
    gamma = goldenSectionMethod(func,l,a,b);
    
else
    
    % Armijo's Condition
    
    s = 0.9;
    beta = 0.8;
    alpha = 0.05;
    
    m_k = 0;
    gamma = s*beta^m_k;
    
    while g(x_k,y_k)-g(x_k+gamma*d_k(1) , y_k+gamma*d_k(2) ) < -alpha*gamma*d_k'*gradg(x_k,y_k)
        m_k = m_k + 1;
        gamma = s*beta^m_k;
    end
    
end

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
    trace_f(i) = g(x(i),y(i));
end

figure(plotNum)
hold on;
plot3(x,y,trace_f,"-r+",'linewidth',1)

% Plot minimum point
hold on;
plot3(x(end),y(end),g(x(end),y(end)),"-r*",'linewidth',7)

end

% Plot Objetive Function through Iterations
function objFuncReps(x,y,k,plotNum,subNum)

trace_g = [];
k = min(k,50);

for i = 1:k
    trace_g(i) = g(x(i),y(i));
end

figure(plotNum)
subplot(3,1,subNum)
plot(0:k-1,trace_g,"-r+",'linewidth',1)
title(["Initial Point: $$(x,y)$$=("+num2str(x(1))+","+num2str(y(1))+")"],"interpreter","latex")
xlabel("Iterations k")
ylabel("$$g(x,y)$$","interpreter","latex")

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