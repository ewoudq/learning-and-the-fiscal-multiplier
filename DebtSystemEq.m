function output = DebtSystemEq(input,shocks,states,beliefs,loc_x_obs,loc_z_obs,phi_mat,names_z_obs,parameters)
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
rhob=parameters(26);
% auxiliary parameters
phi1=phi*(1-sigma)-1;
phi2=(1-phi)*(1-sigma)*n_ss/(1-n_ss);
ty=(1/beta-1)*by+tauc*(1-gy-iy)+tauk*rk_ss*ky+tauw*wny-gy;
ty=(1-1/beta)*by+tauc*(1-gy-iy)+tauk*rk_ss*ky+tauw*wny-gy; % EQ correction oct 2018

% steady state


% ** Input **
%['b ';'k ';'r ';'y ';'c ';'d ';'i ';'mc';'n ';'pi';'q ';'rk';'tr';'w ';'g ';'ur';'z ']; 

b=input(1); %different from tr adjustment
k=input(2);
r=input(3);
y=input(4);
c=input(5);
d=input(6);
i=input(7);
mc=input(8);
n=input(9);
pi=input(10);
q=input(11);
rk=input(12);
tr=input(13); %different from tr adjustment
w=input(14);

% ** Predetermined variables **
endo_states=states(1:m_states);
exo_states=states(1+m_states:end);
k_lagged=endo_states(2);
r_lagged=endo_states(3);
b_lagged=endo_states(1); %different from tr adjustment
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

% accounting identity
output(1)=-y+cy*c+gy*g+iy*i;
% production function
output(2)=-y+z+alpha*k_lagged+(1-alpha)*n;
% wage
output(3)=-w+alpha*k_lagged-alpha*n+z+mc;
% rental rate
output(4)=-rk+z+(alpha-1)*k_lagged+(1-alpha)*n+mc;
% capital accumulation
output(5)=-k+(1-delta)*k_lagged+delta*i;
% labour
if financing==2
    output(6)=c+n*n_ss/(1-n_ss)+gy*g/(wny*(1-tauw))-w;
else
    output(6)=c+n*n_ss/(1-n_ss)-w;    
end

% consumption
sub=2;
switch sub
    case 1
        output(7)=phi1*beliefs(5,:)*XX-phi2*beliefs(9,:)*XX...
            +r-beliefs(10,:)*XX-phi1*c+phi2*n;
    
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
    
    case 2
        betatilde=1/(rk_ss+1-delta);
        wy = ((1-alpha)*(eps-1))/(eps*n_ss);
        xi=cy+wy*(1-n_ss);
        
        if financing==1 ||financing==0% lump-sum tax financing and capital tax financing
            %New:
            G1=(1+(1-phi)*(1+tauc)/(phi*(1-tauw)))*cy/(1-betatilde);
            G2=wy-(1+(1-phi)*(1+tauc)/(phi*(1-tauw)))*beta*cy*(1-sigma)*(1-phi)/(sigma*(1-betatilde));
            G3=(1+(1-phi)*(1+tauc)/(phi*(1-tauw)))*cy*betatilde/(sigma*(1-betatilde));
            G4=wy+(1+(1-phi)*(1+tauc)/(phi*(1-tauw)))*cy*(1-phi)*(1-sigma)/(sigma);
        
            SW=G4*beliefs(14,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
            SG=gy*betatilde*rhog*g/(1-betatilde*rhog);
            SR=G3*beliefs(3,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
            SPi=G3*beliefs(10,:)*inv(eye(nobs)-betatilde*TT)*XX;
            Srk=ky*rk_ss*beliefs(12,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
            SD=dy*beliefs(6,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
        
            output(7)=-G1*c+ky*betatilde^(-1)*k_lagged+G2*w+ky*rk_ss*rk...
                +dy*d-gy*g-G3*r+SW-SR+SPi+Srk+SD-SG;
            %Old
            %output(7)=-xi*(sigma*(1-betatilde)-betatilde*phi1)*c/(sigma*(1-betatilde))...
            %    +ky*betatilde^(-1)*k_lagged+wy*w+ky*rk_ss*rk-gy*g-xi*betatilde*(phi2*n+r)/(sigma*(1-betatilde))...
            %    +(wy+xi*(1-phi)*(1-sigma)/sigma)*betatilde*beliefs(13,:)*inv(eye(nobs)-betatilde*TT)*XX...
            %    -gy*betatilde*rhog*g/(1-betatilde*rhog)...
            %    +dy*d+dy*betatilde*beliefs(5,:)*inv(eye(nobs)-betatilde*TT)*XX... % new line for imperfect competition: dividents
            %    +xi*betatilde/(sigma*(1-betatilde))*(beliefs(9,:)*inv(eye(nobs)-betatilde*TT)*XX-betatilde*beliefs(2,:)*inv(eye(nobs)-betatilde*TT)*XX)...
            %    +ky*rk_ss*betatilde*beliefs(11,:)*inv(eye(nobs)-betatilde*TT)*XX; %adjusted line for adj costs
            
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
            SW=G4*beliefs(14,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
            %SG=gy*betatilde*rhog*g/(1-betatilde*rhog);
            SG=G6*gy*betatilde*rhog*g/(1-betatilde*rhog);
            SR=G3*beliefs(3,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
            SPi=G3*beliefs(10,:)*inv(eye(nobs)-betatilde*TT)*XX;
            Srk=ky*rk_ss*beliefs(12,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
            SD=dy*beliefs(6,:)*betatilde*inv(eye(nobs)-betatilde*TT)*XX;
        
            output(7)=-G1*c+ky*betatilde^(-1)*k_lagged+G2*w+ky*rk_ss*rk...
                +dy*d-G5*gy*g-G3*r+SW-SR+SPi+Srk+SD-SG;
            
            %Old
            %output(7)=-xi*(sigma*(1-betatilde)-betatilde*phi1)*c/(sigma*(1-betatilde))...
            %    +ky*betatilde^(-1)*k_lagged+wy*w+ky*rk_ss*rk...
            %    -(1-tauw*n_ss)*gy*g/((1-tauw)*n_ss)... %-gy*g
            %    -xi*betatilde*(phi2*n+r)/(sigma*(1-betatilde))...
            %    +(wy+xi*(1-phi)*(1-sigma)/sigma)*betatilde*beliefs(13,:)*inv(eye(nobs)-betatilde*TT)*XX...
            %    -((1-tauw*n_ss)/((1-tauw)*n_ss)+xi*(1-sigma)*(1-phi)/(sigma*(1-tauw)*wy*n_ss))*gy*betatilde*rhog*g/(1-betatilde*rhog)... %-gy*betatilde*rhog*g/(1-betatilde*rhog)
            %    +dy*d+dy*betatilde*beliefs(5,:)*inv(eye(nobs)-betatilde*TT)*XX... % new line for imperfect competition: dividents
            %    +xi*betatilde/(sigma*(1-betatilde))*(beliefs(9,:)*inv(eye(nobs)-betatilde*TT)*XX-betatilde*beliefs(2,:)*inv(eye(nobs)-betatilde*TT)*XX)...
            %    +ky*rk_ss*betatilde*beliefs(11,:)*inv(eye(nobs)-betatilde*TT)*XX; %adjusted line for adj costs
        end
end

% dividends
output(8)=-dy*d+y-wny*(w+n)-rk_ss*ky*(rk+k_lagged);

% investment

switch sub
    case 1
output(9)=-q-r+beliefs(10,:)*XX+beta*(rk_ss*(1-tauk)*beliefs(12,:)*XX...
    +(1-delta)*beliefs(11,:)*XX+sigmai*delta^2*(beliefs(7,:)*XX-k));
    %case 2
    %    output(9)=-k+k_lagged+(1-beta*(1-delta))/sigmai*beliefs(10,:)*inv(eye(nobs)-beta*TT)*XX...
    %        -r/sigmai-beta/sigmai*beliefs(2,:)*inv(eye(nobs)-beta*TT)*XX+beliefs(9,:)/sigmai*inv(eye(nobs)-beta*TT)*XX;
    %case 2
    %    bk=by*eps*(1/beta-1+delta)/(alpha*(eps-1)*(1-tauk));
    %    gk=gy*eps*(1/beta-1+delta)/(alpha*(eps-1)*(1-tauk));
    %    output(9)=-k+k_lagged+(1/sigmai)*(...
    %        +beta*rk_ss*beliefs(11,:)*inv(eye(nobs)-beta*TT)*XX...
    %        -gk*beta*rhog*g/(1-beta*rhog)...            
    %        +beta*rk_ss*tauk*k...
    %        +beta^2*tauk*rk_ss*beliefs(1,:)*inv(eye(nobs)-beta*TT)*XX...
    %        +(1+bk)*(beliefs(9,:)*inv(eye(nobs)-beta*TT)*XX...
    %        -beliefs(2,:)*beta*inv(eye(nobs)-beta*TT)*XX-r));
    case 2
        if financing==1
        output(9)=-q-r-beliefs(3,:)*beta*inv(eye(nobs)-beta*TT)*XX...
            +beliefs(10,:)*inv(eye(nobs)-beta*TT)*XX...
            +beta*(1-tauk)*rk_ss*beliefs(12,:)*inv(eye(nobs)-beta*TT)*XX...
            -gy/ky*rhog*beta*g/(1-beta*rhog);
        elseif financing==2 || financing==0
            output(9)=-q-r-beliefs(3,:)*beta*inv(eye(nobs)-beta*TT)*XX...
            +beliefs(10,:)*inv(eye(nobs)-beta*TT)*XX...
            +beta*(1-tauk)*rk_ss*beliefs(12,:)*inv(eye(nobs)-beta*TT)*XX;
        end
end

% taylor rule
output(10)=-r+rhopi*pi+ur;

% Phillips curve

switch sub
    case 1
        output(11)=-pi+(1-theta)*(1-beta*theta)*mc/theta+beta*beliefs(10,:)*XX;
    case 2 % expectations about wages, rk, and technology
        output(11)=-pi+(1-theta)*(1-beta*theta)*mc/theta+beta*(1-theta)*(1-beta*theta)*(...
            (1-alpha)*beliefs(14,:)*inv(eye(nobs)-beta*theta*TT)*XX...
            +alpha*beliefs(12,:)*inv(eye(nobs)-beta*theta*TT)*XX...
            -rhoz*z/(1-beta*theta*rhoz))...
            +beta*(1-theta)*beliefs(10,:)*inv(eye(nobs)-beta*theta*TT)*XX;
    case 3 % expectations about marginal costs
        output(11)=-pi+(1-theta)*(1-beta*theta)*mc/theta+beta*(1-theta)*(...
            (1-beta*theta)*beliefs(8,:)*inv(eye(nobs)-beta*theta*TT)*XX...
            +beliefs(10,:)*inv(eye(nobs)-beta*theta*TT)*XX);
end

%output(12)=-tk*tauk*rk_ss*ky+g*gy+by*(r_lagged-pi)/beta-tauk*rk_ss*ky*(k_lagged+rk);
%output(12)=-tr*ty-by*(r_lagged-pi)/beta+tauk*rk_ss*ky*(k_lagged+rk)...
%    +tauc*(1-gy-iy)*c+tauw*wny*(w+n)-g*gy;
if financing==0
    %output(12)=-tr*ty-by*(r_lagged-pi)/beta+tauk*rk_ss*ky*(k_lagged+rk)...
    %+tauc*(1-gy-iy)*c+tauw*wny*(w+n)-g*gy;
    output(12)=-by*(r_lagged-pi)/beta+tauk*rk_ss*ky*(k_lagged+rk)...
    +tauc*(1-gy-iy)*c+tauw*wny*(w+n)-g*gy-ty*tr+by*b-by*b_lagged/beta; %different from tr adjustment
else
    %output(12)=-tr*ty-by*(r_lagged-pi)/beta+tauk*rk_ss*ky*(k_lagged+rk)...
    %+tauc*(1-gy-iy)*c+tauw*wny*(w+n);
    output(12)=-by*(r_lagged-pi)/beta+tauk*rk_ss*ky*(k_lagged+rk)...
    +tauc*(1-gy-iy)*c+tauw*wny*(w+n)-ty*tr+by*b-by*b_lagged/beta; %different from tr adjustment
end
% Tobin's Q
output(13)=-q+delta*sigmai*(i-k_lagged);

% Debt stabilisation
output(14)=-tr-rhob*b_lagged;

%fminsearch
%output=sum(abs(output));
end

