function f = fitness(x)
f = (exp(x(1) - 1) - 1)^2 + (exp(2*x(1) - 2) - 1)^2 + x(2)^2;