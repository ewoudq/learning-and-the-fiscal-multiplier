function output = NeoClassicalCFunc(input,shocks,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% INPUTS
%   input: contemporaneous endogenous variables
%   TT: transition (state) matrix
%   XX: state vector
%   beliefs: regression coefficients in forecasting model.
%        rows = observed states
%        columns = all endogenous variables
%   parameters: vector of structural parameters
%

% ** Initialization **
m_obs=length(loc_x_obs);
k_obs=length(loc_z_obs);
nobs=m_obs+k_obs;
m_states=length(states)-length(shocks);

% ** State transition **
TT=[beliefs(loc_x_obs,:);
    zeros(k_obs,m_obs) phi_mat(loc_z_obs,loc_z_obs)];

% ** Parameter values **
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

% auxiliary parameters
phi1=phi*(1-sigma)-1;
phi2=(1-phi)*(1-sigma)*n_ss/(1-n_ss);
ty=gy+(1/beta-1)*by;

% steady state


% ** Input **
k=input(1);
r=input(2);
y=input(3);
c=input(4);
d=input(5);
i=input(6);
n=input(7);
q=input(8);
rk=input(9);
t=input(10);
w=input(11);

% ** Predetermined variables **
endo_states=states(1:m_states);
exo_states=states(1+m_states:end);
k_lagged=endo_states(1);
r_lagged=endo_states(2);
g_lagged=exo_states(1);
z_lagged=exo_states(2);


% ** Shock variables **
eg=shocks(1);
ez=shocks(2);

% ** Exogenous states **
% TFP
z=rhoz*z_lagged+ez;
% government spending
g=rhog*g_lagged+eg;
exo_vars=[g;z];

% ** Observed states **
XX=[input(loc_x_obs)';
    phi_mat(loc_z_obs,loc_z_obs)*exo_vars(loc_z_obs)];


% ** Equations **

% consumption
sub=2;
switch sub
    case 1
        output(7)=phi1*beliefs(4,:)*XX-phi2*beliefs(7,:)*XX...
            +r-phi1*c+phi2*n;
    
        %case 2 %is exactly the same as without debt
    %    wy = ((1-alpha)*(eps-1))/(eps*n_ss);
    %    xi=cy+wy*(1-n_ss);
    %    output(7)=-xi*(sigma*(1-beta)-beta*phi1)*c/(sigma*(1-beta))...
    %        +ky*beta^(-1)*k_lagged+wy*w+ky*rk_ss*rk-gy*g-xi*beta*(phi2*n+r)/(sigma*(1-beta))...
    %        +(wy+xi*(1-phi)*(1-sigma)/sigma)*beta*beliefs(12,:)*inv(eye(nobs)-beta*TT)*XX...
    %        -gy*beta*rhog*g/(1-beta*rhog)...
    %        +dy*d+dy*beta*beliefs(5,:)*inv(eye(nobs)-beta*TT)*XX... % new line for imperfect competition: dividents
    %        +xi*beta/(sigma*(1-beta))*(beliefs(9,:)*inv(eye(nobs)-beta*TT)*XX-beta*beliefs(2,:)*inv(eye(nobs)-beta*TT)*XX)...
    %        +ky*rk_ss*beta*beliefs(10,:)*inv(eye(nobs)-beta*TT)*XX; %adjusted line for adj costs
    
    case 2 %NEW: the same as with lump-sum tax financing
        betatilde=1/(rk_ss+1-delta);
        wy = ((1-alpha)*(eps-1))/(eps*n_ss);
        
        % Old
        %xi=cy+wy*(1-n_ss);
        %output(7)=-xi*(sigma*(1-betatilde)-betatilde*phi1)*c/(sigma*(1-betatilde))...
        %    +ky*betatilde^(-1)*k_lagged+wy*w+ky*rk_ss*rk-gy*g-xi*betatilde*(phi2*n+r)/(sigma*(1-betatilde))...
        %    +(wy+xi*(1-phi)*(1-sigma)/sigma)*betatilde*beliefs(11,:)*inv(eye(nobs)-betatilde*TT)*XX...
        %    -gy*betatilde*rhog*g/(1-betatilde*rhog)...
        %    +dy*d+dy*betatilde*beliefs(5,:)*inv(eye(nobs)-betatilde*TT)*XX... % new line for imperfect competition: dividents
        %    +xi*betatilde/(sigma*(1-betatilde))*(...
        %    -betatilde*beliefs(2,:)*inv(eye(nobs)-betatilde*TT)*XX)...
        %    +ky*rk_ss*betatilde*beliefs(9,:)*inv(eye(nobs)-betatilde*TT)*XX; %adjusted line for adj costs
        
        % Simplified        
        G1=cy/(phi*(1-beta));
        G2=wy-beta*cy*(1-sigma)*(1-phi)/(phi*sigma*(1-beta));
        G3=cy*betatilde/(phi*sigma*(1-betatilde));
        G4=(wy+cy*(1-phi)*(1-sigma)/(phi*sigma));
        
        SW=G4*beliefs(11,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
        SG=gy*betatilde*rhog*g/(1-betatilde*rhog);
        SR=G3*beliefs(2,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
        SPi=NaN;
        Srk=ky*rk_ss*beliefs(9,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
        SD=dy*beliefs(5,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
        
        output.c=(ky*betatilde^(-1)*k_lagged+G2*w+ky*rk_ss*rk...
            +dy*d-gy*g-G3*r+SW-SR+Srk+SD-SG)/G1;
        
        % investment
        %output(9)=-q-r-beliefs(2,:)*beta*inv(eye(nobs)-beta*TT)*XX...            
        %    +beta*rk_ss*beliefs(9,:)*inv(eye(nobs)-beta*TT)*XX;
        InvR=-beliefs(2,:)*beta*inv(eye(nobs)-beta*TT)*XX;
        Invrk=beta*rk_ss*beliefs(9,:)*inv(eye(nobs)-beta*TT)*XX;
end

output.G1=G1;
output.G2=G2;
output.G3=G3;
output.G4=G4;
output.SW=SW;
output.SG=SG;
output.SR=SR;
output.SPi=SPi;
output.Srk=Srk;
output.SD=SD;

output.InvR=InvR;
output.Invrk=Invrk;
output.InvPI=NaN;
%fminsearch
%output=sum(abs(output));
end

