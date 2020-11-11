function [epsr,mur,sigma] = fdtd_profile(profile, N, b, varargin)

epsr  = ones(N,1);
mur   = ones(N,1);
sigma = zeros(N,1);

if (profile == 0)        
    % eps = eps_o, mu = mu_o, sigma = 0.
    return
end

for i = 1:N
    if (profile == 1) % Permittivity discontinuity
        if i >= b
            epsr(i) = 9;
        end
    end
    if (profile == 2) % Permeability discontinuity
        if i >= b
            mur(i) = 9;
        end
    end    
    if (profile == 3) % Conducting half space
        if i >= b
            sigma(i) = 1.0;
        end
    end
end