function [f_all,gnorm_all] = gradmeth(fun, x0, tol, maxit)
% code for gradient method, including backtracking line search, goes here,
% instead of the following dummy code, which just takes one gradient step
% and does not even check whether fun(x) < fun(x0).
% Note that the input parameter 'fun' is an anonymous function.
%
[f0,g0] = fun(x0);  % call the anonymous function to get f0 and g0
x = x0 - g0;        % a gradient step
f_all(1) = fun(x);         % evaluate fun at new x
gnorm_all(1) = norm(g0);   % we took one iteration with no line search
