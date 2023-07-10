%% Fiscal Policy under Adaptive Learning: Optimality conditions
% New Keynesian model with distortionary taxation
% © 2017 Ewoud Quaghebeur

%% Optimality conditions
% The matrices AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, and NN of the dynamical system
%
% BLOCK 1: $\quad 0 = AA x_t + BB x_{t-1} + CC y_t + DD z_t$
%
% BLOCK 2: $\quad 0 = E_t [ FF x_{t+1} + GG x_{t} + HH x_{t-1} + JJ y_{t+1} + KK y_{t} +
%       LL z_{t+1} + MM z_t]$
%
% BLOCK 3: $\quad z_{t+1} = NN z_t + \epsilon_{t+1}$ 
% 
% with $E_t[\epsilon_{t+1}] = 0$,
%
% see the |RE_Toolbox| function and the documentation of Uhlig's toolbox
% for more information.
%

% for b(t),k(t),r(t),y(t)
AA = [ 0,0,0,-1
       0,0,0,-1
       0,0,0,0
       0,0,0,0
       0,-1,0,0
       0,0,0,0
       0,0,0,1
       0,0,-1,0
       by,0,0,0
       0,0,0,0
       0,0,0,0];

% for b(t-1),k(t-1),r(t-1),y(t-1)
BB = [ 0,0,0,0
       0,alpha,0,0
       0,alpha,0,0
       0,(alpha-1),0,0
       0,(1-delta),0,0
       0,0,0,0
       0,-rk_ss*ky,0,0
       0,0,0,0
       -by/beta,tauk*rk_ss*ky,-by/beta,0
       0,-delta*sigmai,0,0
       -rhob,0,0,0];

%Order:	cons        div   inv           mc   labour         pi          q       rk              tr      wage
CC = [	(1-gy-iy),  0,    iy,           0,	0,              0,          0,      0,              0,      0
      	0,          0,    0,            0,	(1-alpha),      0,          0,      0,              0,      0      	  
      	0,          0,    0,            1,	-alpha,         0,          0,      0,              0,      -1         
     	0,          0,    0,            1,	(1-alpha),      0,          0,      -1,             0,      0     
      	0,          0,    delta,        0,	0,              0,          0,      0,              0,      0
        1,          0,    0,            0,	n_ss/(1-n_ss),  0,          0,      0,              0,      -1
        0,          -dy,  0,            0,	-wny,           0,          0,      -rk_ss*ky,      0,      -wny
        0,          0,    0,            0,  0,              rhopi,      0,      0,              0,      0
        tauc*(1-gy-iy),0, 0,            0,  tauw*wny,       by/beta,    0,      tauk*rk_ss*ky,  -ty,    tauw*wny
        0,          0,    delta*sigmai, 0,  0,              0,          -1,     0,              0,      0
        0,          0,    0,            0,  0,              0,          0,      0,              -1,     0];   

%      g(t) ur(t)z(t)
DD = [ gy,  0,  0
       0,   0,  1
       0,   0,  1
       0,   0,  1
       0,   0,  0
       0,   0,  0 
       0,   0,  0
       0,   1,  0
       0,   0,  0
       0,   0,  0
       0,   0,  0];

if financing==2 % labour tax financing
    DD(6,1)=gy/(wny*(1-tauw));
elseif financing ==0 % lump-sum financing
    DD(9,1)=-gy;
end
%      b(t+1),k(t+1),r(t+1),y(t+1)
FF = [0,0,0,0
      0,0,0,0
      0,0,0,0];

%      b(t),k(t),r(t),y(t)
GG = [0,0,1,0
      0,-beta*sigmai*delta^2,-1,0
      0,0,0,0];

HH = [0,0,0,0
      0,0,0,0
      0,0,0,0];

%Order:	cons                div	inv                 mc    labour                              pi          q                 rk                  tr   	
JJ = [  (phi*(1-sigma)-1),  0,  0,                  0,    -(1-sigma)*(1-phi)*n_ss/(1-n_ss),   -1,         0,                0,                  0,      0
        0,                  0,  beta*sigmai*delta^2,0,    0,                                  1,          beta*(1-delta),   beta*rk_ss*(1-tauk),0,      0
        0,                  0,  0,                  0,    0,                                  beta,       0,                0,                  0,      0];

%Order:	cons                div inv   mc                                labour                              pi      q   rk	 tr     wage 
KK = [  -(phi*(1-sigma)-1), 0,  0,    0,                                (1-sigma)*(1-phi)*n_ss/(1-n_ss),    0,      0,  0,	 0,     0
        0,                  0,  0,    0,                                0,                                  0,      -1, 0,   0,     0
        0,                  0,  0,    (1-theta)*(1-beta*theta)/theta,   0,                                  -1,     0,  0,   0,     0];

%       g(t+1) ur(t+1) z(t+1)
LL = [  0,0,0
        0,0,0
        0,0,0];

%       g(t) ur(t) z(t)
MM = [  0,0,0
        0,0,0
        0,0,0];
    
if financing==1 % capital tax financing
    MM(2,1)=-beta*gy/ky*rhog;
end   
% Rho   eg      er      ez
NN = [  rhog,   0,      0
        0,      rhor,   0
        0,      0,      rhoz];
    
% Sdev      eg      er      ez
Sigma = [   sigmag, 0,      0
            0,      sigmar, 0
            0,      0,      sigmaz];

% l_equ = number of equations in BLOCK 1
% m_states = number of state variables x
% n_endog = number of other endogenous variables y
% k_exog = number of stochastic processes
[l_equ,m_states] = size(AA);
[l_equ,n_endog] = size(CC);
[l_equ,k_exog] = size(DD);
[q_expectational_equ,m_states] = size(FF);
m_obs=size(names_x_obs,2);
k_obs=size(names_z_obs,2);

%%
% Matrices correspond to the dynamical system in the book of Evans and
% Honkapohja (2001, p. 237):
%
% $$Y_{t} = \alpha + \beta*E_{t}Y_{t+1} + \delta*Y_{t-1} + \kappa*z_{t}$$
%
% $$z_{t} = \phi*z_{t-1} + \epsilon_{t}$$
%
% Since there is no intercept, alpha drops out. But for generality I inlude
% it here nevertheless.

alpha_mat = zeros(m_states+n_endog,1); 
beta_mat = -[AA CC;GG KK]\[zeros(l_equ,m_states) zeros(l_equ,n_endog); FF JJ];
delta_mat = -[AA CC;GG KK]\[BB zeros(l_equ,n_endog); HH zeros(q_expectational_equ,n_endog)];
kappa_mat = -[AA CC;GG KK]\[DD;MM];
phi_mat = NN;
