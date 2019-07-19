function p = henyey_greenstein( g, theta )
 
% function p = henyey_greenstein( g, theta )
% -------------------------------------------------------------------------
% Henyey-Greenstein phase function for asymmetry parameter g and phase
% angle theta
 
p = (1-g.*g)./(1+g.*g-2*g.*cos(theta)).^(1.5);
