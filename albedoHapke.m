function A = albedoHapke(wmix, inc, emi, theta, gmix) 

% function A = albedoHapke(w, B0, h, b, c, inc, emi, theta)
% -------------------------------------------------------------------------
% Hapke's bidirectional reflectance function. q, w, g, n, and a are the
% extinction efficiency, single-scattering albedo, asymmetry parameter,
% number density, and radius of the particles. inc and emi 
% are the incidence and emergence angles (radians). theta is the phase 
% angle, in radians (zero for backscatter).

%variables to change based on parameters:
%size of particle
%a = 1^-5;
%n = number density


%------------------------------------------------------------------%
% Free parameter describing the opposition effect (Hapke, 1981)
B0 = 0.9;

mu0 = cos(inc);
mu = cos(emi);

%h = (3/8) * q.^(3/2);
%h = (1/2)*(n*pi*(10^-5)*q.^(3/2));
h = .003;

B = B0 .* (1 + tan(abs(theta/2))./h).^-1;

% phase function
%p = phaseHapke(theta);
c = .4;
p = henyey_greenstein_two(gmix, c, theta);

    
    
% Chandra H-functions
Hmu0 = H_twostream(mu0, wmix);
Hmu = H_twostream(mu, wmix);

A = (wmix/(4*pi)).*(mu0./(mu0+mu)).*( (1+B).*p + Hmu0.*Hmu - 1 );
end
