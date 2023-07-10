%% Toolbox for finding the Rational Expectations Equilibrium
% 2014 Ewoud Quaghebeur

%%
function [PP,QQ,RR,SS,NN,Omega,Phi] = REE_Toolbox(AA,BB,CC,DD,FF,GG,HH,JJ,KK,LL,MM,NN,Sigma,VARNAMES)

%% Discription
% This function calls part of the Uhlig Toolbox 4.1. (Author: Harald Uhlig)
% to solve a linear rational expectations model for which the matrices of the 
% first order conditions are given.

%% Calls 
% Calls 'options.m', 'solve.m', 'calc_qrs.m', and 'solve_qz.m' from the Uhlig Toolbox.

%% Inputs
% * The matrices AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, and NN of the dynamical system
%
% $$0 = AA x_t + BB x_{t-1} + CC y_t + DD z_t$$
%
% number of rows = l = "l_equ"
% 
% $$0 = E_t [ FF x_{t+1} + GG x_{t} + HH x_{t-1} + JJ y_{t+1} + KK y_{t} +
%       LL z_{t+1} + MM z_t]$$
%
% number of rows = n+m-l = "q_expectational_equ"
%
% $$z_{t+1} = NN z_) + \epsilon_{t+1} $$ 
% 
% number of rows = k (number of exogenous processes) = "k_exog"
%
% with $E_t[\epsilon_{t+1}] = 0$,
%
% * Sigma: variance-covariance matrix of the exogenous variables
% * VARNAMES: character array of variable names in the following order:
%       * Endogenous state variables "x_t"
%       * Endogenous other variables "y_t"
%       * Exogenous state variables  "z_t"

%% Outputs
% * The matrices PP, QQ, RR, SS, NN of the dynamical system
%
% $$x_t = PP x_{t-1} + QQ z_t$$
%
% $$y_t = RR x_{t-1} + SS z_t$$
%
% $$z_t = NN z_{t-1} + \epsilon_t$$
%
% * The matrices "Omega" and "Phi" of the system
%
% $$X_t = \Omega X_{t-1} + \Phi\epsilon_t$$
%
% with $X_t = [x_t', y_t', z_t']'$ and $\Phi$ =[QQ;SS;I(k)]
%

%% Function 
% Calls part of the Uhlig Toolbox 4.1. (Author: Harald Uhlig)
%

warnings = [];
DISPLAY_IMMEDIATELY = 1;
options;
solve;

% Matrices for X(t) = Omega X(t-1) + Phi Epsilon(t)
II_lag = [ PP, zeros(m_states,n_endog),zeros(m_states,k_exog)
         RR, zeros(n_endog, n_endog),zeros(n_endog, k_exog)
         zeros(k_exog,(m_states+n_endog)), NN                ];
II_contemp = eye(m_states+n_endog+k_exog) + ...
   [ zeros(m_states,(m_states+n_endog)), QQ
     zeros(n_endog, (m_states+n_endog)), SS
     zeros(k_exog,  (m_states+n_endog)), zeros(k_exog,k_exog) ];
 
Omega = II_contemp*II_lag;
Phi(:,1:k_exog) = [QQ;SS;eye(k_exog)];

end

