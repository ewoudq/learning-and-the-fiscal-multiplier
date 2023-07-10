function Output = NeoClassical(parameters,irf_shock,sigmag,sigmaz,T,S,draw,itermax,damp,constant,gain,model)
%Neoclassical Summary of this function goes here
%   Detailed explanation goes here

%% Parameters
alpha=parameters(1);
beta=parameters(2);
delta =parameters(3);
eps =parameters(4);
rhog =parameters(5);
rhoz =parameters(6);
sigma =parameters(7);
gy =parameters(8);
rk_ss =parameters(9);
iy =parameters(10);
cy =parameters(11);
ky =parameters(12);
wny =parameters(13);
dy =parameters(14);
n_ss =parameters(15);
phi =parameters(16);
sigmai =parameters(17);
by=parameters(21);
ty=gy+(1/beta-1)*by;

%% Variable names
VARNAMES = ['k ';'r ';'y ';'c ';'d ';'i ';'n ';'q ';'rk';'t ';'w ';'g ';'z ']; 
names_x = {'k','r','y'}; 
names_y = {'c','d','i','n','q','rk','t','w'}; 
names_z = {'g','z'}; 
names_disturbances = {'eg','ez'}; 
names_x_obs={'k'};
%names_x_obs={'k','r'}; 
names_z_obs={'z'}; 
%names_z_obs={'g','z'}; 
names_endogenous=[names_x names_y];
model=1;

OrderAndReorder;

%%
NeoClassicalOptimalityConditions;

%% RE solution

[PP,QQ,RR,SS,NN,Omega,Phi] = RE_Toolbox(AA,BB,CC,DD,FF,GG,HH,JJ,KK,LL,MM,NN,Sigma,VARNAMES);

% *Impulse response analysis*

for i=1:size(irf_shock,2)
    shock=zeros(k_exog,T);
    shock(loc_shock(i),1)=Sigma(loc_shock(i),loc_shock(i));
    eval(['Output.REE_IRF_' names_disturbances{loc_shock(i)} '= RE_Simul(Omega,Phi,T,shock,cellstr(VARNAMES),draw);']);
end

%% 
% _RPE initial beliefs_
%
[a_RPE,m_RPE,c_RPE,beta0_RPE,Omega_RPE,Pi_RPE] = BetaRPE(alpha_mat,beta_mat,delta_mat,kappa_mat,phi_mat,m_states,n_endog,k_exog,loc_x_obs,loc_z_obs,itermax,damp,Sigma,constant,draw);

%% 
% _MSV representation_ of the solution under rational expectations.
%
% $$Y_{t} = b_{msv}*x_{t-1} + c_{msv}*z_{t}$$
%
% $$z_{t} = \varphi*z_{t-1} + \epsilon_{t}$$

a_msv = zeros(m_states+n_endog,1);
b_msv = [PP;RR];
c_msv = [QQ;SS];
rho=NN;

% The number of columns of b_msv and c_msv equals the nr of states
% ("m_states") and the nr of shocks ("k_exog") respectively. The matrices
% b_msv_full and c_msv_full have a column for each variable of the system
% (the colums corresponding to other exogenous variables). 
b_msv_full = [b_msv zeros(m_states+n_endog,n_endog+k_exog)];
c_msv_full = [c_msv zeros(m_states+n_endog,m_states+n_endog)];

%%
if constant==1
    if loc_z_obs==0 % if constant and no observed shocks
        beta0_msv=[a_msv'; b_msv(:,loc_x_obs)'];
    else % if constant and observed shocks
        beta0_msv=[a_msv'; b_msv(:,loc_x_obs)'; c_msv(:,loc_z_obs)'];        
    end
else
    if loc_z_obs==0 % if no constant and no observed shocks
        beta0_msv=b_msv(:,loc_x_obs)'; 
    else % if no constant and observed shocks
        beta0_msv=[b_msv(:,loc_x_obs)'; c_msv(:,loc_z_obs)'];
    end
end

%% Rational Expectations Simulation

% Random deviates:random draws for the disturbances from a normal 
% distribution with zero mean and second moments matrix Sigma.
rng(120)% set seed
epsilon=zeros(k_exog,S);
for i=1:size(Sigma,2)
    epsilon(i,:)=normrnd(0,Sigma(i,i),1,S);
end
[REE_Simul] = RE_Simul(Omega,Phi,S,epsilon,cellstr(VARNAMES),draw);
Output.REE_Simul = REE_Simul;
%% 
% _Regression-based initial beliefs_
%
loc_z_obs_xy=loc_z_obs+m_states+n_endog; %location of observed shocks in xy-vector
if constant==1
    if loc_z_obs==0 % if constant and no observed shocks
        beta0_res=[a_msv'; b_msv(:,loc_x_obs)'];
        data_res = [ones(S-1,1) REE_Simul.y(2:end,loc_x_obs)];
    else % if constant and observed shocks
        beta0_res=[a_msv'; b_msv(:,loc_x_obs)'; c_msv(:,loc_z_obs)'];        
        data_res = [ ones(S-1,1) REE_Simul.y(2:end,loc_x_obs) REE_Simul.y(1:end-1,loc_z_obs_xy)];        
    end
else
    if loc_z_obs==0 % if no constant and no observed shocks
        beta0_res=b_msv(:,loc_x_obs)';
        data_res = REE_Simul.y(2:end,loc_x_obs);
    else % if no constant and observed shocks
        beta0_res=[b_msv(:,loc_x_obs)'; c_msv(:,loc_z_obs)'];
        data_res =[REE_Simul.y(2:end,loc_x_obs) REE_Simul.y(1:end-1,loc_z_obs_xy)];
    end
end
R0_res=cov(data_res); % initial second moments vector

%% Learning


%load('R0_res2.mat');
%R0_res=[R0_res R0_res(:,3)];
%R0_res=[R0_res; R0_res(1,:)];

beta0=beta0_RPE;
for x=1:size(irf_shock,2)
    shocks=zeros(k_exog,T);
    shocks(loc_shock(x),1)=Sigma(loc_shock(x),loc_shock(x));
    eval(['[Output.IHL_' names_disturbances{loc_shock(x)} '] = IHLearning(beta0,R0_res,names_x,names_y,names_z,names_x_obs,names_z_obs,shocks,phi_mat,parameters,constant,gain,model);']);
    eval(['[Output.RL_' names_disturbances{loc_shock(x)} '] = RestrictedLearning(beta0,R0_res,names_x,names_y,names_z,names_x_obs,names_z_obs,shocks,T,alpha_mat,beta_mat,delta_mat,kappa_mat,phi_mat,gain,constant,draw);']);
end

%% Plot
if draw==1
PlotComparison;
end

end

