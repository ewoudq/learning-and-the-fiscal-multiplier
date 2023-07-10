%% Present-value multiplier
% © 2015 Ewoud Quaghebeur

%%
function [multiplier,x_pv,x_pvsum] = PVMultiplier(y,g,r,ratio,horizon)
%PVMultiplier Calculates the present-value y-multiplier of the fiscal
%instrument g.
%   Following Mountford and Uhlig (2009), this function calculates the
%   present-value multiplier
%
%   $$\frac{\sum_{j=0}^k{R^{-j}y_j}}{\sum_{j=0}^k{R^{-j}g_j}}\frac{1}{g/y}$$
%
%   at lag 0 up to _horizon_. 
%
%   * $y_j$ is the response of the variable of interest at period $j$.
%   * $g_j$ is the response of the fiscal variable of interest at period
%   $j$.
%   * $R$ is the interest rate. $R$ can also be time-varying.
%   * $g/y$ (_ratio_) is the ratio of the fiscal variable to the variable
%   y.
%

% check if r is a scalar or an array
if isscalar(r)
    r=ones(horizon,1)*r;
end

x=[y g];
x_pv=nan(horizon,2);
x_pvsum=nan(horizon,2);

for i=1:2
    x_pv(1,i)=x(1,i);
    x_pvsum(1,i)=x(1,i);
    for k=2:horizon
        x_pv(k,i)=x(k,i)*r(k)^(-(k-1));
        x_pvsum(k,i)= x_pvsum(k-1,i)+x_pv(k,i);
    end
end

multiplier=nan(horizon,1);
for k=1:horizon
    multiplier(k)=x_pvsum(k,1)/x_pvsum(k,2)*(1/ratio);
end
end