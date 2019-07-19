function p = henyey_greenstein_two(g, c, theta)
 
% function p = henyey_greenstein_two(g, c, theta)
% -------------------------------------------------------------------------
% Two-lobe Henyey-Greenstein phase function
 
x = henyey_greenstein(g, theta);
y = henyey_greenstein(-g, theta);
 
p = 0.5 * ( (1-c).*x + (1+c).*y );
