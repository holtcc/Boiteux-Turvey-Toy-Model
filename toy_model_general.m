% Toy model for analyzing questions pertaining to efficient investment in
% electricity markets.

clear all

pbar = 2;
v = 3; %2
s_optimal = v - pbar;
c = 1;
% if using a renewable MC > 0, make less than c (to preserve merit order)
c_r = 0;

epsilon = 0.3;

% set optimal R, G (R listed first always)
rtemp = 0.4;
gtemp = 1.2;

% define PDF and CDF here:
% note upper bound on integral is now infinity
loadpdf = @(x) pdf('LogNormal',x,0,1);
loadcdf = @(x) cdf('LogNormal',x,0,1);
% define other functions:
evfun = @(x) x.*loadpdf(x);

% use first order conditions to back out g, r that produce interior
% solution above

% quadratic inv. costs (no intercept term):
g = ( (pbar + s_optimal - c)*[1 - loadcdf(epsilon*rtemp +gtemp) ] ) / gtemp
r = ( (c*epsilon*[loadcdf(epsilon*rtemp + gtemp) - loadcdf(epsilon*rtemp)]) + (pbar + s_optimal)*epsilon*[1 - loadcdf(epsilon*rtemp + gtemp)] - c_r*epsilon*[1 - loadcdf(epsilon*rtemp)] ) / rtemp


% constant marginal inv. costs:
% g = ( (pbar + s_optimal - c)*[1 - loadcdf(epsilon*rtemp +gtemp) ] ) 
% r = ( (c*epsilon*[loadcdf(epsilon*rtemp + gtemp) - loadcdf(epsilon*rtemp)]) + (pbar + s_optimal)*epsilon*[1 - loadcdf(epsilon*rtemp + gtemp)] - c_r*epsilon*[1 - loadcdf(epsilon*rtemp)] )

start = [0.5,0.4];
options = optimoptions('fsolve','Display','iter');

% check optimal solution to make sure same as rtemp, gtemp
xstar_optimal = fsolve(@tm_system_general,start,options,c,c_r,g,r,pbar,epsilon,s_optimal,loadcdf)
[rtemp gtemp]

% this is for testing domain of s
fsolve(@tm_system_general,start,options,c,c_r,g,r,pbar,epsilon,3,loadcdf)
% tm_system_cinvc(start,c,c_r,g,r,pbar,epsilon,2.5)

lb = 0; 
ub = 3;
svec = lb:0.01:ub;

N = length(svec)

Rstar = 999*ones(N,1);
Gstar = 999*ones(N,1);
totcap = 999*ones(N,1);

llp = 999*ones(N,1);
vll = 999*ones(N,1);

utility = 999*ones(N,1);
energy_exp = 999*ones(N,1);
cap_pmt = 999*ones(N,1);
cs = 999*ones(N,1);

gencosts = 999*ones(N,1);
invcosts = 999*ones(N,1);
ps = 999*ones(N,1);

W = 999*ones(N,1);
TC = 999*ones(N,1);


for i = 1:N
    s = svec(i);
    xstar = fsolve(@tm_system_general,start,options,c,c_r,g,r,pbar,epsilon,s,loadcdf);
    
    Rstar(i) = xstar(1);
    Gstar(i) = xstar(2);
    R = xstar(1);
    G = xstar(2);

    totcap(i) = epsilon*R + G;
    llp(i) = (1 - loadcdf(epsilon*R + G));

    vllfun = @(x) v*(x - epsilon*R - G).*loadpdf(x);
    vll(i) = integral(vllfun, (epsilon*R + G) , Inf);    
 
    utility(i) = v*integral(evfun,0,Inf) - vll(i);

    energy_exp(i) = c_r*integral(evfun,0,epsilon*R) + c*integral(evfun,epsilon*R,epsilon*R + G) + pbar*(epsilon*R + G)*(1 - loadcdf(epsilon*R + G));
    
    cap_pmt(i) = s*(epsilon*R + G)*(1 - loadcdf(epsilon*R + G));

    cs(i) = utility(i) - energy_exp(i) - cap_pmt(i);

    midcostfun = @(x) (x - epsilon*R).*loadpdf(x);

    gencosts(i) = c_r*integral(evfun,0,epsilon*R) + c_r*epsilon*R*(1 - loadcdf(epsilon*R)) ...
        + c*integral(midcostfun,epsilon*R,epsilon*R + G) + c*G*(1 - loadcdf(epsilon*R + G));
        
    invcosts(i) = g*G + r*R;
    
    ps(i) = energy_exp(i) + cap_pmt(i) - gencosts(i) - invcosts(i);

    W(i) = cs(i) + ps(i);
    TC(i) = vll(i) + gencosts(i) + invcosts(i);
    
end

plot(svec,round(Rstar,3))

plot(svec,Gstar)


plot(svec,totcap)

% makes sense that this should be zero when constant marginal inv. cost?

plot(svec,cs)
plot(svec,W)
plot(svec,round(ps,5))

max(W)
cs(W==max(W))
svec(W==max(W))

plot(svec,energy_exp)
plot(svec,cap_pmt)
plot(svec,energy_exp+cap_pmt)
plot(svec,utility)
plot(svec,vll)