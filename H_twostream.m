function H = H_twostream(x, w)

% function H_twostream(x, w)
% -------------------------------------------------------------------------
% Two-stream approximation for Chandrasekhar's H-function. x = cos(theta)
% is the cosine of the photometric angle, and w is the single scattering
% albedo

% "albedo factor"
gamma = sqrt(1-w);

H = (1+2*x)./(1+2*gamma.*x);
end
