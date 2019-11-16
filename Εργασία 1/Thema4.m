% Konstantinos Letros 8851
% Optimization Techniques
% Project 1 - Bisection Derivative Method

%% Clean the screen

clc
clear
close all;
format long;

%% Functions to be tested

f = cell(3,1);

f{1} = @(x) sqrt(x+1)+(x^2-2)*log(x+1);
f{2} = @(x) (x-2)^2+cos(x)^2;
f{3} = @(x) exp(-3*x)-(sin(x-2)-2)^2;

%% Derivatives of the functions to be tested

dfdx = cell(length(f),1);

syms x

for i = 1:length(f)
    dfdx{i} = matlabFunction( diff(f{i}(x)) );
end

%% Parameters

% Interval [a,b]
a_init = 0;
b_init = 4;

% l values to be tested
l_val = 3e-3:4e-2:5e-1;

%% Bisection Derivative Method

% Count Plots
plotCounter = 1;
c = [];
for i = 1 : length(f)
    k=[];
    %%% Correction to interval for function 3 due to convexity issues
    if i == 3
        b_init = 3;
    end
    %%%
    
    counter = plotCounter;
    for l = l_val
        [a,b,k(end+1)] = bisectionDerivativeMethod(dfdx{i},l,a_init,b_init);
        
        figure(counter)
        sgtitle('Bisection Derivative Method (Intervals)')
        
        subplot(3,1,i)
        plot(1:k(end),a,"-x")
        
        hold on;
        plot(1:k(end),b,"-o")
        title(['$$ f_',num2str(i),'(x) = $$',funcToLatex(f{i}),'  at l = ',num2str(l)],'Interpreter','Latex')
        xlabel("Repetitions until convergence - k")
        ylabel('k_{th} Interval [a_k,b_k]')
        
        counter = counter + 1;
    end
    
    
    figure(counter)
    
    sgtitle('Bisection Derivative Method (Accuracy)')
    subplot(3,1,i)
    plot(l_val,k,'-x')
    title(['$$ f_',num2str(i),'(x) = $$',funcToLatex(f{i})],'Interpreter','Latex')
    xlabel("Accuracy - l")
    ylabel("Calculations until convergence")
    
     c(:,i) = k;
end

plotCounter = counter;

%% Save Plots

% for i = 1 : plotCounter
%     figure(i)
%     savePlot([mfilename,'_',num2str(i)])
% end

%% Functions

function [a,b,k] =   bisectionDerivativeMethod(dfdx,l,a,b)
k = 1;
n = ceil(log2((b(k)-a(k))/l));

for i=1:n
    x(k) = (a(k)+b(k))/2;
    if dfdx(x(k)) == 0
        break;
    elseif dfdx(x(k))>0
        a(k+1) = a(k);
        b(k+1) = x(k);
    else
        a(k+1) = x(k);
        b(k+1) = b(k);
    end
    k = k + 1;
end

end

function expr = funcToLatex(func)
% Remove arguments, element-wise operators
fstr = regexprep(func2str(func), '^@\(.*?\)|\.', '');
% Convert to symbolic expression
fsym = feval(symengine, 'hold', fstr);
expr = ['$$ ',latex(fsym),' $$'];
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