function output = DistortionaryCFunc(input,shocks,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters)
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
theta =parameters(18);
rhopi =parameters(19);
rhor =parameters(20);
by=parameters(21);
tauk=parameters(22);
tauw=parameters(23);
tauc=parameters(24);
financing=parameters(25);
% auxiliary parameters
phi1=phi*(1-sigma)-1;
phi2=(1-phi)*(1-sigma)*n_ss/(1-n_ss);
ty=(1/beta-1)*by+tauc*(1-gy-iy)+tauk*rk_ss*ky+tauw*wny-gy;
ty=(1-1/beta)*by+tauc*(1-gy-iy)+tauk*rk_ss*ky+tauw*wny-gy; % EQ correction oct 2018

% steady state


% ** Input **
%['k ';'r ';'y ';'c ';'d ';'i ';'mc';'n ';'pi';'rk';'t ';'w ';'g ';'ur';'z '] 
k=input(1);
r=input(2);
y=input(3);
c=input(4);
d=input(5);
i=input(6);
mc=input(7);
n=input(8);
pi=input(9);
q=input(10);
rk=input(11);
tr=input(12);
w=input(13);

% ** Predetermined variables **
endo_states=states(1:m_states);
exo_states=states(1+m_states:end);
k_lagged=endo_states(1);
r_lagged=endo_states(2);
g_lagged=exo_states(1);
ur_lagged=exo_states(2);
z_lagged=exo_states(3);


% ** Shock variables **
eg=shocks(1);
er=shocks(2);
ez=shocks(3);

% ** Exogenous states **
% TFP
z=rhoz*z_lagged+ez;
% government spending
g=rhog*g_lagged+eg;
% monetary shock
ur=rhor*ur_lagged+er;
exo_vars=[g;ur;z];

% ** Observed states **
XX=[input(loc_x_obs)';
    phi_mat(loc_z_obs,loc_z_obs)*exo_vars(loc_z_obs)];


% ** Equations **

% consumption
%sub=2;
%switch sub
    %case 1
        %output.c=phi1*beliefs(4,:)*XX-phi2*beliefs(8,:)*XX...
            %+r-beliefs(9,:)*XX-phi1*c+phi2*n;  

    %case 2
        betatilde=1/(rk_ss+1-delta);
        wy = ((1-alpha)*(eps-1))/(eps*n_ss);
        xi=cy+wy*(1-n_ss);
        
        if financing==1 ||financing==0% lump-sum tax financing and capital tax financing
            %New:
            G1=(1+(1-phi)*(1+tauc)/(phi*(1-tauw)))*cy/(1-betatilde);
            G2=wy-(1+(1-phi)*(1+tauc)/(phi*(1-tauw)))*beta*cy*(1-sigma)*(1-phi)/(sigma*(1-betatilde));
            G3=(1+(1-phi)*(1+tauc)/(phi*(1-tauw)))*cy*betatilde/(sigma*(1-betatilde));
            G4=wy+(1+(1-phi)*(1+tauc)/(phi*(1-tauw)))*cy*(1-phi)*(1-sigma)/(sigma);
        
            SW=G4*beliefs(13,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
            SG=gy*betatilde*rhog*g/(1-betatilde*rhog);
            SR=G3*beliefs(2,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
            SPi=G3*beliefs(9,:)*inv(eye(nobs)-betatilde*TT)*XX;
            Srk=ky*rk_ss*beliefs(11,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
            SD=dy*beliefs(5,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
        
            output.c=(ky*betatilde^(-1)*k_lagged+G2*w+ky*rk_ss*rk...
                +dy*d-gy*g-G3*r+SW-SR+SPi+Srk+SD-SG)/G1;
   
        elseif financing==2
            %New:
            G1=(1+(1-phi)*(1+tauc)/(phi*(1-tauw)))*cy/(1-betatilde);
            G2=wy-(1+(1-phi)*(1+tauc)/(phi*(1-tauw)))*beta*cy*(1-sigma)*(1-phi)/(sigma*(1-betatilde));
            G3=(1+(1-phi)*(1+tauc)/(phi*(1-tauw)))*cy*betatilde/(sigma*(1-betatilde));
            G4=wy+(1+(1-phi)*(1+tauc)/(phi*(1-tauw)))*cy*(1-phi)*(1-sigma)/(sigma);
            G5=(1-tauw*n_ss)/((1-tauw)*n_ss)-(1+(1-phi)*(1+tauc)/(phi*(1-tauw)))*...
                betatilde*(1-sigma)*phi*(1-n_ss)/(sigma*(1-betatilde)*(1+tauc)*n_ss);
            G6=(1+(1-phi)*(1+tauc)/(phi*(1-tauw)))*(1-sigma)*(1-n_ss)*phi/(sigma*(1+tauc)*n_ss)...
                +(1-tauw*n_ss)/((1-tauw)*n_ss);
            output.G5=G5;
            output.G6=G6;
            SW=G4*beliefs(13,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
            %SG=gy*betatilde*rhog*g/(1-betatilde*rhog);
            SG=G6*gy*betatilde*rhog*g/(1-betatilde*rhog);
            SR=G3*beliefs(2,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
            SPi=G3*beliefs(9,:)*inv(eye(nobs)-betatilde*TT)*XX;
            Srk=ky*rk_ss*beliefs(11,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
            SD=dy*beliefs(5,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
        
            output.c=(ky*betatilde^(-1)*k_lagged+G2*w+ky*rk_ss*rk...
                +dy*d-G5*gy*g-G3*r+SW-SR+SPi+Srk+SD-SG)/G1;
               
        end
%end

% investment

%switch sub
    %case 1
%output(9)=-q-r+beliefs(9,:)*XX+beta*(rk_ss*(1-tauk)*beliefs(11,:)*XX...
    %+(1-delta)*beliefs(10,:)*XX+sigmai*delta^2*(beliefs(6,:)*XX-k));

    %case 2
        InvR=-beliefs(2,:)*beta*inv(eye(nobs)-beta*TT)*XX;
        InvPI=beliefs(9,:)*inv(eye(nobs)-beta*TT)*XX;
        Invrk=beta*(1-tauk)*rk_ss*beliefs(11,:)*inv(eye(nobs)-beta*TT)*XX;
        
        if financing==1
        %output(9)=-q-r-beliefs(2,:)*beta*inv(eye(nobs)-beta*TT)*XX...
        %    +beliefs(9,:)*inv(eye(nobs)-beta*TT)*XX...
        %    +beta*(1-tauk)*rk_ss*beliefs(11,:)*inv(eye(nobs)-beta*TT)*XX...
        %    -gy/ky*rhog*beta*g/(1-beta*rhog);
        
        InvG=-gy/ky*rhog*beta*g/(1-beta*rhog);
        
        elseif financing==2 || financing==0
        InvG=NaN;    
        %    output(9)=-q-r-beliefs(2,:)*beta*inv(eye(nobs)-beta*TT)*XX...
        %    +beliefs(9,:)*inv(eye(nobs)-beta*TT)*XX...
        %    +beta*(1-tauk)*rk_ss*beliefs(11,:)*inv(eye(nobs)-beta*TT)*XX;
         
        end
%end

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
output.InvPI=InvPI;
output.Invrk=Invrk;
output.InvG=InvG;
end

