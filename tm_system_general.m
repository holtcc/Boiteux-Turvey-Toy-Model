function F = tm_system_general(x,c,c_r,g,r,pbar,epsilon,s,loadcdf)
R = x(1);
G = x(2);

% constant:
% F(1) = g / (1 - loadcdf(epsilon*R + G) ) + c - pbar - s;
% F(2) = r + c_r*epsilon*(1 - loadcdf(epsilon*R)) - c*epsilon*(loadcdf(epsilon*R+G)-loadcdf(epsilon*R))...
%     - (pbar + s)*epsilon*(1 - loadcdf(epsilon*R + G));

% quadratic:
F(1) = g*G + (c - pbar - s)*(1 - loadcdf(epsilon*R + G) );
F(2) = r*R + c_r*epsilon*(1 - loadcdf(epsilon*R)) - c*epsilon*(loadcdf(epsilon*R+G)-loadcdf(epsilon*R))...
    - (pbar + s)*epsilon*(1 - loadcdf(epsilon*R + G));
