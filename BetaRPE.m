%% RPE initial beliefs
% © 2015 Ewoud Quaghebeur

%%
function [a_RPE,m_RPE,c_RPE,beta0_RPE,Omega_RPE,Pi_RPE] = BetaRPE(alpha_mat,beta_mat,delta_mat,kappa_mat,phi_mat,m_states,n_endog,k_exog,loc_x_obs,loc_z_obs,itermax,damp,Sigma,constant,draw)
%BetaRPE Initial beliefs from the Restricted Perceptions Equilibrium.
%   I use the approach of Guse (JEDC, 2008) to find the RPE coefficients.

%%
% I use the approach of Guse (JEDC, 2008) to find the RPE coefficients.
%

% It is possible to specify the observed variables for each forward-looking
% variable separately in via loc_z and loc_v. Here I assume it is always
% x_obs and z_obs.
loc_z=zeros(m_states+n_endog);
loc_z(:,loc_x_obs)=1;
loc_v=zeros(m_states+n_endog,k_exog);
loc_v(:,loc_z_obs)=1;

% Find RPE coefficients
[a_RPE,m_RPE,c_RPE,Omega_RPE,Pi_RPE] = RPE(alpha_mat,beta_mat,delta_mat,kappa_mat,phi_mat,m_states,n_endog,k_exog,loc_z,loc_v,itermax,damp,Sigma,draw);

%% Recover coefficients w.r.t. the current-period shocks
% The iterative scheme above provides the RPE coefficients for the lagged
% shock variables $w_{t-1}$. Here, I recover the coefficients for the
% current-period shock variables $w_t$.
%

rho_temp=diag(phi_mat);
if sum(rho_temp==0)>0
    warning('RPE. Some shock variables are not auto-correlated.');    
    logicals=(rho_temp~=0)';
    A = c_RPE(:,logicals);
    v = rho_temp(logicals)';
    B = bsxfun(@rdivide, A, v);
    c_RPE(:,logicals)=B;
else
    c_RPE=c_RPE/phi_mat;
end
%c_RPE=c_RPE/phi_mat;
%% Construct beta0_RPE

if constant==1
    if loc_z_obs==0 % if constant and no observed shocks
        beta0_RPE=[a_RPE'; m_RPE(:,loc_x_obs)'];
    else % if constant and observed shocks
        beta0_RPE=[a_RPE'; m_RPE(:,loc_x_obs)'; c_RPE(:,loc_z_obs)'];        
    end
else
    if loc_z_obs==0 % if no constant and no observed shocks
        beta0_RPE=m_RPE(:,loc_x_obs)'; 
    else % if no constant and observed shocks
        beta0_RPE=[m_RPE(:,loc_x_obs)'; c_RPE(:,loc_z_obs)'];
    end
end