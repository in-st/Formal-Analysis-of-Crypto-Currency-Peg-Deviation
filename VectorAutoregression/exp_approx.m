function ex=exp_approx(x,a,h)
% Taylor approximation of exponential function
% Fernando Pérez Forero

ex=exp(a);

for j=1:h
    ex=ex+(1/factorial(j))*(x-a)^j;
end

end