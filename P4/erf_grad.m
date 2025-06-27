% Return the gradient of erf function 

function [y] = erf_grad(x)
    y = (2/sqrt(pi)) * exp(-x^2);
end