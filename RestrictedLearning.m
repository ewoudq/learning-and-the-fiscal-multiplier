%% Recursive Least Squares Learning with Restricted Perceptions
% © 2015 Ewoud Quaghebeur

%%
function [Returns] = RestrictedLearning(beta0,R0,names_x,names_y,names_z,names_x_obs,names_z_obs,epsilon,S,alpha_mat,beta_mat,delta_mat,kappa_mat,phi_mat,gain,constant,draw)

%% Description
% This function extends the _RLSLearning_ function with the possibility
% that the Perceived Law of Motion is underparameterised, because a subset
% of variables (state variables and/or shock variables) is ommited from the data vector in the learning scheme.
%
% The main differences with the _RLSLearning_ function are put in italics.
%
% The function simulates a given macroeconomic model under a (constant
% gain variant) of recursive least squares (RLS) learning. Expectations are
% formed at time t (see Evans and Honkapohja, p. 236).
%
% The function uses the following model:
%
% $$y_t=\alpha+\beta E^*_ty_{t+1}+\delta y_{t-1}+\kappa w_t,\qquad
% w_t=\varphi w_{t-1}+e_t$$
%
% where   
%
% * m = number of endogenous state variables ($x$)
% * n = number of other endogenous variables ($y$)
% * k = number of shocks ($w$)
%
% The expectations $E_t^*y_{t+1}$ are formed using the forecast rule
% (Perceived Law of Motion, PLM)
%
% $$E_ty_{t+1}=a+by_t+c\varphi w_t$$
%
% where the beliefs $a$, $b$, and $c$ are updated according to the RLS
% algorithm (Evans and Honkapohja, 2001, p. 238).
% 
% The dynamics of learning are described by the following Actual Law of
% Motion (ALM):
%
% $$y_t=a_t+b_ty_{t-1}+c_t\varphi w_t$$
%

%% Function inputs
%
% * Initialisation of the recursive algorithm: "beta0" _[(1+m_obs+k_obs,m+n)
% matrix]
% and "R0" [1+m_obs+k_obs square matrix] (second moments matrix of data vector).
% "beta0" are the initial beliefs and have the following structure:
% $\beta_0=[a; b; c]$, where $a$ is $1\times (m+n)$, $b$ is $m_{obs}\times(m+n)$,
% and $c$ is $k_{obs}\times(m+n)$._
% * Names [arrays]:  "names_x" (state variables), "names_y" (other 
% endogenous variables), and "names_z" (shocks)
% * "S": number of simulation periods
% * Disturbances: "epsilon" [(1,S) matrix]
% * Model parameters: "alpha_mat" ($\alpha$), "beta_mat" ($beta$), 
% "delta_mat" ($\delta$),"kappa_mat" ($\kappa$), and "phi_mat" ($\varphi$)
% * Gain parameter "gain". If gain==0, then the learning mechanism is RLS.
% * "draw": Function returns figures with simulated data and beliefs iff
% draw==1.
% * _"names_x_obs" and "names_z_obs" are the names of the observed state
% and shock variables_
%

%% Function outputs
%
% * "y": matrix of simulated endogenous variables (n+m times S)
% * "w": matrix of simulated shocks (k times S)
% * "beta": matrix of beliefs at time S _(1+m_obs+k_obs times n+m)_
% (without "1" if no constant)
% * "a","b","c": the parameters of the ALM (m+n times S)
% * "names_x", "names_y", "names_z", and "names_xy"
%

%% Function

%%
% *Initialisation*
%%
% _ 1. Number of variables_
n=size(names_y,2);
m=size(names_x,2);
k=size(names_z,2);
m_obs=size(names_x_obs,2); %new
k_obs=size(names_z_obs,2); %new
%%
% _2. Names_
names_xy = [names_x names_y];

%%
% _3.Location of observed variables_
%
%   The variables x_obs and z_obs are a subset of x and z. The vectors
%   loc_x_obs and loc_z_obs give the location of this subset within x and
%   z.
loc_x_obs=zeros(1,m_obs);
for i=1:m_obs
    for j=1:m
        if strcmp(names_x_obs(i),names_x(j))            
            loc_x_obs(i)=j;
        end
    end
end
loc_z_obs=zeros(1,k_obs);
for i=1:k_obs
    for j=1:k
        if strcmp(names_z_obs(i),names_z(j))            
            loc_z_obs(i)=j;
        end
    end
end
%%
% _4. Initial beliefs_
%
% $\beta_0$ is $[a; b; c]$ where a, b, and c are the _restricted_ coefficients.
% These initial coefficients are stored in
%
% * $a_0$ is $(n+m)\times 1$
% * $b_0$ is $(n+m)\times m_{obs}$
% * $c_0$ is $(n+m)\times k_{obs}$
% * $b_{0,full}=[b_0 \quad zeros(n+m,n)]$ is $(n+m)\times (n+m)$ ?????
%
% Notice the transposition. Each period the coefficients are stored in
% $a_{new}$, $b_{new}$, $c_{new}$, and $\beta_{new}$. _(should be adapted when m and/or k
% is larger than 1)._
%
% Restricted initial beliefs:
%
if constant==1
    a0=beta0(1,:)';
else
    constant=0;
    a0=zeros(1,n+m)';
end
b0=beta0(1+constant:m_obs+constant,:)';
c0=beta0(m_obs+constant+1:end,:)';
%b0_full = [b0 zeros(n+m,n)]; %???
a_new(:,1)=a0;
b_new(:,1)=reshape(b0,1,(n+m)*m_obs); %=b0;
if ~isempty(c0)
    c_new(:,1)=reshape(c0,1,(n+m)*k_obs);%c_new(:,1)=c0;
else
    c_new(:,1)=reshape(zeros(n+m,k),1,(n+m)*k);%c_new(:,1)=zeros(n+m,k);
end
beta_new=beta0;
R_new = R0; % initial moments matrix

%% 
% *First period*
%
% The function |T_map| returns the coefficients of the Actual Law of Motion
% (ALM) for the first period, based on the initial beliefs. Given the
% vector of disturbances in the first period ($\epsilon(:,1)$), the vector
% of shocks $w_1$ and of endogenous variables ($y_1$, size $m+n$) are
% calculated.
% 
% _Given that the initial beliefs are (potentially) restricted, we construct
% unrestricted initial beliefs to be able to calculate the
% ALM-coefficients._
%
% _Unrestricted initial beliefs_
%
b0_full=zeros(n+m,m+n);
c0_full=zeros(n+m,k);
b0_full(:,loc_x_obs)=b0;
if ~isempty(c0) 
    c0_full(:,loc_z_obs)=c0;
end

[a_alm,b_alm,c_alm] = T_map(a0,b0_full,c0_full,alpha_mat,beta_mat,delta_mat,kappa_mat,phi_mat);
w(:,1) = phi_mat*zeros(k,1) + epsilon(:,1);
y(:,1) = a_alm + b_alm*zeros(n+m,1)+c_alm*w(:,1);
yf(:,1)= a_new(:,1)+b0_full*y(:,1)+c0_full*phi_mat*w(:,1); % expectations
wf(:,1)= phi_mat*w(:,1); % expectations

%% 
% *Algorithm*
%
% From the second period onwards, the loop with the recursive learning
% algorithm is as follows
% 
% # The beliefs are updated using the formulas of the learning rule below.
% # The new beliefs $\beta_{t}$ are stored in $a_{new}$, $b_{new}$, $c_{new}$, and $\beta_{new}$. _(should be adapted when m and/or k
% is larger than 1)._ (See "Initialisation" for the dimensions).
% # The coefficients of the PLM (a,b,c) are mapped to the coefficients of
% the ALM. 
% # Period-t endogenous variables and shocks are realised
%
% The formulas of the learning rule are
%
% $$ R_{t} = R_{t-1} + \gamma*(Z*Z'-R_{t-1})$$
%
% $$ error_{t} = y_{t-1}-\beta_{t-1}'*Z$$
%
% $$ \beta_{t} = \beta_{t-1} + \gamma*R_t^{-1}*Z*error'$$
%
% where the data vector $Z=[1;y_{t-2};w_{t-1}]$. Parameter $\gamma$ is the
% "gain" parameter. 

for t=2:S
    % 1. UPDATING BELIEFS
    beta_old  = beta_new; % previous period beliefs (i.e. beta(t-1))
    R_old = R_new;

    if t==2 % data vector Z
        if ~isempty(c0)
            z = [zeros(m_obs,1); w(loc_z_obs,t-1)];            
        else
            z = zeros(m_obs,1);
        end
    else
        if ~isempty(c0)
            z = [y(loc_x_obs,t-2); w(loc_z_obs,t-1)];
        else
            z = y(loc_x_obs,t-2);
        end
    end
    if constant == 1
      z = [1; z];
    end
    if gain == 0 % if gain = 0, then RLS (decreasing gain) learning is executed.
        g = t^(-1);
    else
        g = gain;
    end    
    R_new = R_old + g*(z*z'-R_old); % Learning rules
    error = y(:,t-1)-beta_old'*z;
    beta_new = beta_old + g*(R_new\(z*error')); %beta_new = beta_old + g*R_new^(-1)*z*error';  
    
    % 2. NEW BELIEFS ARE STORED IN a, b, and c.
    if constant==1
        a_new(:,t)=beta_new(1,:)';
    else
        a_new(:,t)=zeros(1,n+m);
    end
    %b_new(:,t)=beta_new(1+constant:m_obs+constant,:)'; before extension
    b_new_temp=beta_new(1+constant:m_obs+constant,:)'; %this and next line: extention
    b_new(:,t)=reshape(b_new_temp,1,(n+m)*m_obs);
    %b_new_full = [b_new(:,t) zeros(n+m,n)];
    %c_new(:,t)=beta_new(m_obs+2:end,:)';
    if ~isempty(c0)
        %c_new(:,t)=beta_new(m_obs+constant+1:end,:)'; before extension
        c_new_temp=beta_new(m_obs+constant+1:end,:)';
        c_new(:,t)=reshape(c_new_temp,1,(n+m)*k_obs);
    else
        c_new(:,t)=reshape(zeros(n+m,k),1,(n+m)*k);%c_new(:,t)=zeros(n+m,k);
    end
    % New: full
    b_new_full=zeros(n+m,m+n);
    c_new_full=zeros(n+m,k);
    %b_new_full(:,loc_x_obs)=b_new(:,t); before extension
    b_new_full(:,loc_x_obs)=b_new_temp; % after extension
    %c_new_full(:,loc_z_obs)=c_new(:,t);
    if ~isempty(c0) 
        c_new_full(:,loc_z_obs)=c_new_temp; %after extension
        %c_new_full(:,loc_z_obs)=c_new(:,t); before extension
    end    
    % 3. The MAP from PLM to ALM
    [a_alm,b_alm,c_alm] = T_map(a_new(:,t),b_new_full,c_new_full,alpha_mat,beta_mat,delta_mat,kappa_mat,phi_mat);
    
    % 4. REALISATION of contemporaneous shocks and endogenous variables
    w(:,t) = phi_mat*w(:,t-1) + epsilon(:,t);
    y(:,t) = a_alm + b_alm*y(:,t-1)+c_alm*w(:,t);
    yf(:,t)= a_new(:,t)+b_new_full*y(:,t)+c_new_full*phi_mat*w(:,t); % expectations
    wf(:,t)= phi_mat*w(:,t); % expectations
end

%% 
% *Figures*
% 
% If draw==1, the function creates figures with the simulated variables and
% beliefs.

if draw==1
    i=1;% Endogenous variables
    figure('Name','Restricted Learning simulation: endogenous variables');
    for j=1:n+m               
        subplot(3,3,i);plot(1:S,y(j,:));
        title(names_xy(j)); 
        i=i+1;
        if mod(j,9)==0
            figure('Name','Restricted Learning simulation: endogenous variables');
            i=1;
        end 
    end
    i=1;
    figure('Name','Restricted Learning simulation: stochastic processes');
    for j=1:k % Stochastic processes
        subplot(3,3,i);plot(1:S,w(j,:));
        title(names_z(j));
        i=i+1;
        if mod(j,9)==0
            figure('Name','Restricted Learning simulation: stochastic processes');
            i=1;
        end
    end    
% Beliefs 
    figure('Name','Restricted Learning simulation: Beliefs a');%   beliefs a
    plot(1:S,a_new,'-',1:S,a_new(:,1)*ones(1,S),'--');
    title('Beliefs a')
    % beliefs b
    for xx=1:m_obs
        titl = sprintf('Restricted Learning simulation: Beliefs b: state %s.', names_x_obs{xx});
        figure('Name',titl);
        i=1;
        series=b_new((xx-1)*(n+m)+1:(n+m)*xx,:);
        for j=1:size(series,1)
        str = sprintf('Beliefs b for %s.', names_xy{j});    
        %subplot(3,3,i);plot(1:S,series(j,:),'-',1:S,series(j,1)*ones(1,S),'--');
        subplot(3,3,i);plot(1:S,series(j,:),'-',1:S,b0(j,xx)*ones(1,S),'--');
        title(str)
        i=i+1;
        if mod(j,9)==0
            figure('Name',titl);
            i=1;
        end 
        end
    end
    % beliefs c
    for zz=1:k_obs
        titl = sprintf('Restricted Learning simulation: Beliefs c: shock %s.', names_z_obs{zz});
        figure('Name',titl);
        i=1;
        series=c_new((zz-1)*(n+m)+1:(n+m)*zz,:);
        for j=1:size(series,1)
        str = sprintf('Beliefs c for %s.', names_xy{j});
        %subplot(3,3,i);plot(1:S,series(j,:),'-',1:S,series(j,1)*ones(1,S),'--');
        subplot(3,3,i);plot(1:S,series(j,:),'-',1:S,c0(j,zz)*ones(1,S),'--');
        title(str)
        i=i+1;
        if mod(j,9)==0
            figure('Name',titl);
            i=1;
        end 
        end
    end
end

%% 
% *Returns*
Returns.y=y;
Returns.w=w;
Returns.yf=yf;
Returns.wf=wf;
Returns.beta=beta_new;
Returns.a=a_new;
Returns.b=b_new;
Returns.c=c_new;
Returns.names_x=names_x;
Returns.names_y=names_y;
Returns.names_z=names_z;
Returns.names_xy=names_xy;
Returns.names_x_obs=names_x_obs;
Returns.names_z_obs=names_z_obs;
Returns.epsilon=epsilon;
end