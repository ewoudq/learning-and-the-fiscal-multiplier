function [Returns] = IHLearning(beta0,R0,names_x,names_y,names_z,names_x_obs,names_z_obs,epsilon,phi_mat,parameters,constant,gain,model)
%IHLearning Summary of this function goes here
%   Detailed explanation goes here

%%
% *Initialisation*
%%

option = optimoptions('fsolve','Display','off');

% _ 1. Number of variables_
n=size(names_y,2);
m=size(names_x,2);
k=size(names_z,2);
m_obs=size(names_x_obs,2); %new
k_obs=size(names_z_obs,2); %new
T=size(epsilon,2);
nobs=m_obs+k_obs;

w=nan(k,T);
y=nan(n+m,T);
yf=nan(n+m,T);
beliefs=beta0';
x_initial=ones(1,n+m);
%shocks=zeros(k,T);
states=zeros(m+k,1);%[0;0;0;0]; % endogenous states ; exogenous states

% consumption function
c_consfunc=nan(1,T);
SW=nan(1,T);
SG=nan(1,T);
SR=nan(1,T);
SPi=nan(1,T);
Srk=nan(1,T);
SD=nan(1,T);

% investment function
Invrk=nan(1,T);
InvR=nan(1,T);
InvG=nan(1,T);
InvPI=nan(1,T);
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
R_new=R0;

%% * First period *
epsi=epsilon(:,1)';
switch model
    case 1
        FUN=@(input) NeoClassicalSystemEq(input,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
    case 2
        FUN=@(input) NewKeynesianSystemEq(input,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
    case 3
        FUN=@(input) DistortionarySystemEq(input,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
    case 4
        FUN=@(input) DebtSystemEq(input,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
end
[x_res] = fsolve(FUN,x_initial);
%[x_res] = fminsearch(FUN,x_initial);
y(:,1)=x_res';
w(:,1) = phi_mat*zeros(k,1) + epsilon(:,1);

% save terms of the consumption function
switch model
    case 1
        CFunc=NeoClassicalCFunc(x_res,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
        CFunc.InvG=NaN;
    case 2
        CFunc=NewKeynesianCFunc(x_res,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
        CFunc.InvG=NaN;
    case 3
        CFunc=DistortionaryCFunc(x_res,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
        %FUN=@(input) DistortionarySystemEq(input,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
    case 4
        CFunc=DebtCFunc(x_res,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
end
c_consfunc(1)=CFunc.c;
G1=CFunc.G1;
G2=CFunc.G2;
G3=CFunc.G3;
G4=CFunc.G4;
if isfield(CFunc,'G5')
    Returns.G5=CFunc.G5;
end    
if isfield(CFunc,'G6')
    Returns.G6=CFunc.G6;
end 
SW(1)=CFunc.SW;
SG(1)=CFunc.SG;
SR(1)=CFunc.SR;
SPi(1)=CFunc.SPi;
Srk(1)=CFunc.Srk;
SD(1)=CFunc.SD;

% save terms of the investment function
Invrk(1)=CFunc.Invrk;
InvR(1)=CFunc.InvR;
InvPI(1)=CFunc.InvPI;
InvG(1)=CFunc.InvG;
%% Loop

for it=2:T
    % 1. UPDATING BELIEFS
    beta_old  = beta_new; % previous period beliefs (i.e. beta(t-1))
    R_old = R_new;
    
    if it==2 % data vector Z
        if ~isempty(c0)
            z = [zeros(m_obs,1); w(loc_z_obs,it-1)];
        else
            z = zeros(m_obs,1);
        end
    else
        if ~isempty(c0)
            z = [y(loc_x_obs,it-2); w(loc_z_obs,it-1)];
        else
            z = y(loc_x_obs,it-2);
        end
    end
    if constant == 1
        z = [1; z];
    end
    if gain == 0 % if gain = 0, then RLS (decreasing gain) learning is executed.
        g = it^(-1);
    else
        g = gain;
    end
    R_new = R_old + g*(z*z'-R_old); % Learning rules
    err = y(:,it-1)-beta_old'*z;
    yf(:,it-1)=beta_old'*z;
    beta_new = beta_old + g*(R_new\(z*err'));
    
    % 2. NEW BELIEFS ARE STORED IN a, b, and c.
    if constant==1
        a_new(:,it)=beta_new(1,:)';
    else
        a_new(:,it)=zeros(1,n+m);
    end
    %b_new(:,t)=beta_new(1+constant:m_obs+constant,:)'; before extension
    b_new_temp=beta_new(1+constant:m_obs+constant,:)'; %this and next line: extention
    b_new(:,it)=reshape(b_new_temp,1,(n+m)*m_obs);
    %b_new_full = [b_new(:,t) zeros(n+m,n)];
    %c_new(:,t)=beta_new(m_obs+2:end,:)';
    if ~isempty(c0)
        %c_new(:,t)=beta_new(m_obs+constant+1:end,:)'; before extension
        c_new_temp=beta_new(m_obs+constant+1:end,:)';
        c_new(:,it)=reshape(c_new_temp,1,(n+m)*k_obs);
    else
        c_new(:,it)=reshape(zeros(n+m,k),1,(n+m)*k);%c_new(:,t)=zeros(n+m,k);
    end
    
    % 3. REALISATION OF VARIABLES
    beliefs=beta_new';
    epsi=epsilon(:,it)';
    w(:,it) = phi_mat*w(:,it-1) + epsilon(:,it);
    states=[y(1:m,it-1);w(:,it-1)];
    switch model
        case 1
            FUN=@(input) NeoClassicalSystemEq(input,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
    case 2
            FUN=@(input) NewKeynesianSystemEq(input,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
    case 3
            FUN=@(input) DistortionarySystemEq(input,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
    case 4
            FUN=@(input) DebtSystemEq(input,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
    end
    [x_res,fval,exitflag,output] = fsolve(FUN,x_res,option);
    %[x_res,fval,exitflag,output] = fminsearch(FUN,x_res);
    y(:,it)=x_res';
    if exitflag<1
        error('Error. Equation not solved.')
    end
    
    % save terms of the consumption function
    switch model
        case 1
            %FUN=@(input) NeoClassicalSystemEq(input,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
            CFunc=NeoClassicalCFunc(x_res,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
            CFunc.InvG=NaN;
        case 2
            CFunc=NewKeynesianCFunc(x_res,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
            CFunc.InvG=NaN;
        case 3
            CFunc=DistortionaryCFunc(x_res,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
            %FUN=@(input) DistortionarySystemEq(input,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
        case 4
            CFunc=DebtCFunc(x_res,epsi,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters);
    end
    c_consfunc(it)=CFunc.c;
    SW(it)=CFunc.SW;
    SG(it)=CFunc.SG;
    SR(it)=CFunc.SR;
    SPi(it)=CFunc.SPi;
    Srk(it)=CFunc.Srk;
    SD(it)=CFunc.SD;
    Invrk(it)=CFunc.Invrk;
    InvR(it)=CFunc.InvR;  
    InvPI(it)=CFunc.InvPI;
    InvG(it)=CFunc.InvG;
end

% Final yf
if ~isempty(c0)
    z = [y(loc_x_obs,T-1); w(loc_z_obs,T)];
else
    z = y(loc_x_obs,T-1);
end
if constant == 1
    z = [1; z];
end
yf(:,T)=beta_new'*z;

Returns.y=y;
Returns.yf=yf;
Returns.w=w;
Returns.beta=beta_new;
Returns.beta0=beta0;
Returns.R0=R0;
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

%% Parameters of the consumption function

Returns.c_consfunc=c_consfunc;
Returns.G1=G1;
Returns.G2=G2;
Returns.G3=G3;
Returns.G4=G4;
Returns.SW=SW;
Returns.SG=SG;
Returns.SR=SR;
Returns.SPi=SPi;
Returns.Srk=Srk;
Returns.SD=SD;

%% Parameters of the investment function
Returns.Invrk=Invrk;
Returns.InvR=InvR;
Returns.InvPI=InvPI;
Returns.InvG=InvG;
end

