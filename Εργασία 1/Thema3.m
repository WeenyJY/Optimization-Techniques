% Konstantinos Letros 8851
% Optimization Techniques
% Project 1 - Fibonacci Method

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

%% Parameters

% Interval [a,b]
a_init = 0;
b_init = 4;

% l values to be tested
l_val = 3e-3:4e-2:5e-1;

%% Fibonacci Method

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
        [a,b,k(end+1)] = fibonacciMethod(f{i},l,a_init,b_init);
        
        figure(counter)
        sgtitle('Fibonacci Method (Intervals)')
        
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
    % Plot [a,b](k)
    sgtitle('Fibonacci Method (Accuracy)')
    subplot(3,1,i)
    plot(l_val,k-1,'-x')
    title(['$$ f_',num2str(i),'(x) = $$',funcToLatex(f{i})],'Interpreter','Latex')
    xlabel("Accuracy - l")
    ylabel("Calculations until convergence")
     c(:,i) = k-1;
end

plotCounter = counter;

%% Save Plots

% for i = 1 : plotCounter
%     figure(i)
%     savePlot([mfilename,'_',num2str(i)])
% end

%% Functions

function [a,b,k] =  fibonacciMethod(func,l,a,b)

k = 1;
n = 1;

while fibonacci(n) <= (b(k)-a(k))/l
    n = n + 1;
end

x1(k)=a(k)+fibonacci(n-2)/fibonacci(n)*(b(k)-a(k));
x2(k)=a(k)+fibonacci(n-1)/fibonacci(n)*(b(k)-a(k));

for i=1:n-1
    
    if(func(x1(k))<=func(x2(k)))
        a(k+1) = a(k);
        b(k+1) = x2(k);
        x2(k+1) = x1(k);
        x1(k+1) = a(k+1) + fibonacci(n-k-2)/fibonacci(n-k)*(b(k+1)-a(k+1));
    else
        a(k+1) = x1(k);
        b(k+1) = b(k);
        x1(k+1) = x2(k);
        x2(k+1) = a(k+1) + fibonacci(n-k-1)/fibonacci(n-k)*(b(k+1)-a(k+1));
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