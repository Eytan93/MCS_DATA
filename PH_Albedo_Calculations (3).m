function [qmix, wmix, gmix] = miemix( q, w, g, n, a )

% function [qmix, wmix, gmix] = miemix( q, w, g, n, a )
% -------------------------------------------------------------------------
% Calculate Mie parameters for a mixture of 2+ components (e.g. Hapke,
% 1980)
% q[] = extinction efficiencies (from table)
% w[] = single scattering albedos (from table)
% g[] = asymmetry parameters (from table)
% n[] = number densities (try different ones, must be two column array)
% a[] = radii (look up)

nel = size(q,1);

qmix = zeros(nel,1);
wmix = qmix;
gmix = qmix;

for i=1:nel
 qmix(i) = sum( n.*a.^2.*q((i,:)) ) / sum( n.*a.^2 );
 wmix(i) = sum( n.*a.^2.*q(i,:).*w(i,:) ) / sum( n.*a.^2.*q(i,:) );
 gmix(i) = sum( n.*a.^2.*q(i,:).*w(i,:).*g(i,:) ) / sum( n.*a.^2.*q(i,:).*w(i,:) );
end
end

function p = phaseHapke(theta)

% function p = phaseHapke(theta)
% -------------------------------------------------------------------------
% Lunar single-particle phase function from Hapke (1963)

p = (4*pi/5)*( (sin(theta) + (pi-theta)*cos(theta))/pi + ...
 0.1*(1-cos(theta)).^2 );
end


function H = H_twostream(x, w)

% function H_twostream(x, w)
% -------------------------------------------------------------------------
% Two-stream approximation for Chandrasekhar's H-function. x = cos(theta)
% is the cosine of the photometric angle, and w is the single scattering
% albedo

% "albedo factor"
gamma = sqrt(1-w);
x = cos(emi);

H = (1+2*x)./(1+2*gamma.*x);
end

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
%a = 10^-5;
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
