clc
clear
format long;

func1 = @(x) sqrt(x+1)+(x^2-2)*log(x+1);
func2 = @(x) (x-2)^2+cos(x)^2;
func3 = @(x) exp(-3*x)-(sin(x-2)-2)^2;
syms x
dfdx1 = matlabFunction( diff(func1(x)) );
dfdx2 = matlabFunction( diff(func1(x)) );
dfdx3 = matlabFunction( diff(func1(x)) );

e = 1e-4;
l = 0.007;
a_init = 0;
b_init = 4;

[a,b,k] = bisectionMethod(func1,e,l,a_init,b_init);
fprintf("Bisection Method\n")
fprintf("Min in interval : [%f,%f] after %d repetitions\n\n",a(end),b(end),k);

[a,b,k] = goldenSectionMethod(func1,l,a_init,b_init);
fprintf("Golden Section Method\n")
fprintf("Min in interval : [%f,%f] after %d repetitions\n\n",a(end),b(end),k);

[a,b,k] = fibonacciMethod(func1,l,a_init,b_init);
fprintf("Fibonacci Method\n")
fprintf("Min in interval : [%f,%f] after %d repetitions\n\n",a(end),b(end),k);

[a,b,k] = BisectionDerivativeMethod(dfdx1,l,a_init,b_init);
fprintf("Bisection Derivative Method\n")
fprintf("Min in interval : [%f,%f] after %d repetitions\n\n",a(end),b(end),k);


function [a,b,k] =  bisectionMethod(func,e,l,a,b)

k = 1;

while b(k)-a(k)>=l
    x1(k) = (a(k) + b(k))/2 - e;
    x2(k) = (a(k) + b(k))/2 + e;
    
    if func(x1(k))<=func(x2(k))
        a(k+1) = a(k);
        b(k+1) = x2(k);
    else
        a(k+1) = x1(k);
        b(k+1) = b(k);
    end
    k = k + 1;
end

end


function [a,b,k] =  goldenSectionMethod(func,l,a,b)

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

end

function [a,b,k] =  fibonacciMethod(func,l,a,b)

k = 1;
n = 1;

while fibonacci(n+1) <= (b(k)-a(k))/l
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
%     if k == n-2
%         x1(k+1) = x1(k);
%         x2(k+1) = x1(k) + e;
%         if(func(x1(k))<=func(x2(k)))
%             a(k+1)=a(k);
%             b(k+1)=x2(k);
%         else
%             a(k+1)=x1(k);
%             b(k+1)=b(k);
%         end
%         break;
%     end
    k = k + 1;
end

end

function [a,b,k] =   BisectionDerivativeMethod(dfdx,l,a,b)
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