%% Restricted Perceptions Equilibrium
% © 2015 Ewoud Quaghebeur

%%
function [a_RPE,m_RPE,c_RPE,Omega,Pi] = RPE(alpha_mat,beta_mat,delta_mat,kappa_mat,phi_mat,m_states,n_endog,k_exog,loc_z,loc_v,itermax,damp,var_eps,draw)
%RPE Finds the Restricted Perceptions Equilibrium using the SUR method of
%Guse (JEDC, 2008).
%   Detailed explanation goes here. 

nz=m_states+n_endog; % sum of state and other endogenous variables
nv=k_exog; % number of shock variables
sig=10e-6; % significance level of the iterative scheme

%% Write the dynamical system in the appropriate form
%
% Evans and Honkapohja (2001, p. 237) assume the following representation
% of the model:
%
% $$Y_{t} = \alpha + \beta*E_{t}Y_{t+1} + \delta*Y_{t-1} + \kappa*z_{t}$$
%
% $$z_{t} = \phi*z_{t-1} + \epsilon_{t}$$
%
% Since there is no intercept, alpha drops out. But for generality I inlude
% it here nevertheless.
%
% Here, I write the system as
%
% $$z_t=A+BE_t^*z_t+DE_t^*z_{t+1}+Fz_{t-1}+Gv_t $$
% $$v_t=\rho v_{t-1}+\epsilon_t$$
%

A=alpha_mat;
B=zeros(m_states+n_endog);
D=beta_mat;
F=delta_mat;
G=kappa_mat;%eye(nz); %#EQ?
rho=phi_mat;

%% Unprojected Actual Law of Motion
% 
% The (unprojected) Actual Law of Motion is
%
% $z_t=T'_a+T'_mz_{t-1}+h(T'_cv_{t-1})+G\epsilon_t$
% 
a=zeros(nz,1);
m=zeros(nz); 
c=zeros(nz,nv); 
Ta=A+B*a+D*(eye(size(m))+m)*a;
Tm=B*m+D*m^2+F;
Tc=B*c+D*(m*c+c*rho)+G*rho;

%% Restrictions
% 
% The restrictions on the variables used by the agents in their restricted
% forecasting model are passed through the function via _loc_z_ and loc_v_.
% The variables used in the forecasting model have the following
% dimensions:
%
% * $v_{i,t-1}$ is the $(k_i\times 1)$ vector of shocks used to estimate
% variable i.% 
% * $w_{i,t-1}$ is the $(\tilde{k}_i\times 1)$ vector of endogenous
% variables used to estimate variable i.
%
% I use the following integers:
%
% * $k=\sum_{i=1}^n{k_i}$%
% * $\tilde{k}=\sum_{i=1}^n{\tilde{k}_i}$
% * $l=k+\tilde{k}$
%
k=sum(loc_z(:)~=0); % k_i's
ktilde=sum(loc_v(:)~=0); % \tilde{k}_i's
ks=sum(loc_z~=0,2); % k
ktildes=sum(loc_v~=0,2); % \tilde{k}
l=sum(loc_z(:)~=0)+sum(loc_v(:)~=0);

%% Matrix related to error term variance

TG=[G;eye(nv)];
GvarG=TG*var_eps*TG';

%% p-map
% Following Guse (2008) we calculate the "projected T-map( or "p-map")
% which maps the restricted PLM to the projected ALM."
%
%[pA1,pB1,Omega1,Pi1] = pmap2005(a,m,c,rho,GvarG,loc_z,loc_v,A,B,D,F,G,nz,nv);


%% Iterative scheme to find the RPE coefficients

% preparations

Tmvec=Tm(:);
Tcvec=Tc(:);
Tx=[Tmvec;Tcvec];

%test
Tmvec=Tm';
Tcvec=Tc';
Tmvec=Tmvec(:);
Tcvec=Tcvec(:);
Tx=[Tmvec;Tcvec];

vecTc=[Tmvec(logical(loc_z));Tcvec(logical(loc_v))];%Tx(loc_k~=0);
a=nan(nz,itermax);
c=zeros(nz,nv,itermax);
b=nan(l,itermax);
mat=nan(nz+nz*size(c,2)+nz*size(m,2),itermax);
a(:,1)=zeros(nz,1);
c(:,:,1)=Tc;
b(:,1)=vecTc;
m=zeros(nz,nz,itermax);
m(:,:,1)=Tm;
mat(:,1)=[a(:,1);reshape(c(:,:,1),size(a,1)*size(c,2),1);reshape(m(:,:,1),size(a,1)*size(m,2),1)];

% iterative scheme
for i=2:itermax
    % p-mapping
    [pA,pB,Omega,Pi] = pmap2005(a(:,i-1),m(:,:,i-1),c(:,:,i-1),rho,GvarG,loc_z,loc_v,A,B,D,F,G,nz,nv);

    % Updating
    a(:,i)=a(:,i-1)+damp*(pA-a(:,i-1));
    b(:,i)=b(:,i-1)+damp*(pB-b(:,i-1));
    
    % recover c and b from Bcal      
    mvec=b(1:k,i);
    for j=1:nz
       m(j,logical(loc_z(j,:)),i)=mvec(1:ks(j))';
       mvec(1:ks(j))=[];
    end    
    cvec=b(k+1:end,i);
    for j=1:nz
       c(j,logical(loc_v(j,:)),i)=cvec(1:ktildes(j))';
       cvec(1:ktildes(j))=[];
    end
    
    % Matrix:
    mat(:,i)=[a(:,i);reshape(c(:,:,i),size(a,1)*size(c,2),1);reshape(m(:,:,i),size(a,1)*size(m,2),1)];    

    conv=sum(abs(mat(:,i)-mat(:,i-1)));
    %fprintf('Iteration %i, convergence %.3e\n',i-1,conv);
    if conv <sig
        fprintf('RPE found at iteration %i (convergence: %.3e)\n',i-1,conv);
        a_RPE=a(:,i);
        c_RPE=c(:,:,i);
        m_RPE=m(:,:,i);
        %beta_RPE=[a_RPE m_RPE c_RPE];
        break        
    end
    if i==itermax
        disp('No convergence');
    end
end
Tm_RPE=B*m_RPE+D*m_RPE^2+F;
Tc_RPE=B*c_RPE+D*(m_RPE*c_RPE+c_RPE*rho)+G*rho;

Tm_RPEvec=Tm_RPE';
Tc_RPEvec=Tc_RPE';
Tm_RPEvec=Tm_RPEvec(:);
Tc_RPEvec=Tc_RPEvec(:);
Tx_RPE=[Tm_RPEvec;Tc_RPEvec];

a(:,i+1:end)=[];
b(:,i+1:end)=[];
c(:,:,i+1:end)=[];
m(:,:,i+1:end)=[];

%% Figure 
if draw==1
    figure;plot(b');title('RPE coefficients');
end
end

