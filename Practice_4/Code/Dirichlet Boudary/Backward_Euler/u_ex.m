function u_exact = u_ex(x,t)
% u_exact = ((x^2)/2 - (x^3)/3 - 1/12)*exp(-2*t); 
u_exact = exp(-1/4*pi^2*t)*sin(2*pi*x);
end