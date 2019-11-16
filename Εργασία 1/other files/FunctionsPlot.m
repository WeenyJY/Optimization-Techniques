clc
clear
close all;
format long;

f = cell(3,1);

f{1} = @(x) sqrt(x+1)+(x^2-2)*log(x+1);
f{2} = @(x) (x-2)^2+cos(x)^2;
f{3} = @(x) exp(-3*x)-(sin(x-2)-2)^2;

dfdx = cell(length(f),1);

syms x

for i = 1:length(f)
    dfdx{i} = matlabFunction( diff(f{i}(x)) );
end


x_min = fminbnd(f{1},0,4);
y_min = f{1}(x_min);
t1=-0.9:1e-4:5;
figure(1);
plot(t1,arrayfun(f{1},t1))
title(['Minimum of f_1(x) : f_1(x*) = ', num2str(y_min) ,' at x* = ',num2str(x_min)])
xlabel("x")
ylabel(['$$ f_1(x) = $$',funcToLatex(f{1})],'Interpreter','Latex')
hold on;
xline(0, 'r');
xline(4, 'b');
hold on;
plot(t1,arrayfun(dfdx{1},t1))
hold on;
plot(x_min,y_min,'*r')
legend('Function f_1(x)','Derivative of f_1(x)','Boundary b','Boundary a')


x_min = fminbnd(f{2},0,4);
y_min = f{2}(x_min);
t2=-3:1e-4:7;
figure(2);
plot(t2,arrayfun(f{2},t2))
title(['Minimum of f_2(x) : f_2(x*) = ', num2str(y_min) ,' at x* = ',num2str(x_min)])
xlabel("x")
ylabel(['$$ f_2(x) = $$',funcToLatex(f{2})],'Interpreter','Latex')
hold on;
xline(0, 'r');
xline(4, 'b');
hold on;
plot(t2,arrayfun(dfdx{2},t2))
hold on;
plot(x_min,y_min,'*r')
legend('Function f_2(x)','Derivative of f_2(x)','Boundary b','Boundary a')


x_min = fminbnd(f{3},0,3);
y_min = f{3}(x_min);
t3=-1:1e-4:8;
figure(3);
plot(t3,arrayfun(f{3},t3))
title(['Minimum of f_3(x) : f_3(x*) = ', num2str(y_min) ,' at x* = ',num2str(x_min)])
xlabel("x")
ylabel(['$$ f_3(x) = $$',funcToLatex(f{3})],'Interpreter','Latex')
hold on;
xline(0, 'r');
xline(3, 'b');
hold on;
plot(t3,arrayfun(dfdx{3},t3))
hold on;
plot(x_min,y_min,'*r')
legend('Function f_3(x)','Derivative of f_3(x)','Boundary b','Boundary a')


for i = 1 : 3
    figure(i)
    savePlot([mfilename,'_',num2str(i)])
end


function expr = funcToLatex(func)
% Remove arguments, element-wise operators
fstr = regexprep(func2str(func), '^@\(.*?\)|\.', '');
% Convert to symbolic expression
fsym = feval(symengine, 'hold', fstr);
expr = ['$$ ',latex(fsym),' $$'];
end

function savePlot(name)

% Resize current figure to fullscreen for higher resolution image
set(gcf, 'Position', get(0, 'Screensize'));

% Save current figure with the specified name
saveas(gcf, join([name,'.jpg']));

% Resize current figure back to normal
set(gcf,'position',get(0,'defaultfigureposition'));

end