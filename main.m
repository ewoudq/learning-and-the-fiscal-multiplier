%% Learning and the size of the Government Spending Multiplier
% © 2017 Ewoud Quaghebeur

%% Parameters
alpha = 1/3;                % Output elasticity with respect to capital
beta = 1.04^(-0.25);        % Households subjective discount factor
gain=.02;                   % Gain parameter
delta = 0.025;              % Rate of physical capital depreciation
eps=6;                      % Elasticity of substitution between intermediate goods
theta=.75;                  % Degree of nominal price rigidity
rhog = .9;                  % Government expenditure AR(1) coefficient
rhopi=1.5;                  % Taylor rule inflation rate coefficient
rhor = 0.5;                 % Interest rate AR(1) coefficient
rhoz = .9;                  % Technology shock AR(1) coefficient
sigma = 2;%2;                  % Coefficient of risk aversion
sigmag=.05;                 % Standard deviation of government spending shock.  Units: Percent.
sigmar=.05;                 % Standard deviation of interest rate shock. Units: Percent.
sigmaz=.05;                 % Standard deviation of technology shock. Units: Percent.
sigmai =17;                 % Capital adjustment cost parameter
gy =.2;                     % Steady state government expenditure to output ratio

rk_ss=(beta^(-1) - 1 + delta);% steady state rental rate of captal

mc_ss=(eps-1)/eps;          % steady state marginal cost
ky=mc_ss*alpha/rk_ss;
iy=delta*ky; % I/Y ratio ! different than without imperfect competition
cy=1-gy-iy;                 % C/Y ratio
n_ss=1/3;                   % Steady state labour
wy_ss = ((1-alpha)*(eps-1))/(eps*n_ss);  % steady-state wage rate
by = 0.74;                  % Steady state debt to gdp ratio
phi = cy/(wy_ss*(1-n_ss)+cy); % Preference parameter

wny=(1-alpha)*mc_ss;
dy=1-wny-rk_ss*ky;

ty=gy+(1/beta-1)*by;

parameters=[alpha,beta,delta,eps,rhog,rhoz,sigma,gy,...
    rk_ss,iy,cy,ky,wny,dy,n_ss,phi,sigmai,theta,rhopi,rhor,by];
%% Learning parameters
constant=0;
gain=0.02;

%% Programme settings
irf_shock={'eg'}; %{'eg','er','ez'};
draw=0;
T=40;
S=10000;
itermax=10000;
damp=.9;

%% Neoclassical model
NeoClassicalResults=NeoClassical(parameters,irf_shock,sigmag,sigmaz,T,S,draw,itermax,damp,constant,gain);

%% Plots
Plots(NeoClassicalResults.IHL_eg.y,NeoClassicalResults.IHL_eg.w,...
    NeoClassicalResults.IHL_eg.beta0,NeoClassicalResults.IHL_eg.a,...
    NeoClassicalResults.IHL_eg.b,NeoClassicalResults.IHL_eg.c,...
    NeoClassicalResults.IHL_eg.names_x,...
    NeoClassicalResults.IHL_eg.names_y,...
    NeoClassicalResults.IHL_eg.names_z,...
    NeoClassicalResults.IHL_eg.names_x_obs,...
    NeoClassicalResults.IHL_eg.names_z_obs,constant);

%% Expectations
selected_vars={'r','rk','w','d'};
NeoClassicalExp=PlotExp(NeoClassicalResults,selected_vars);

% coordinates forecast errors
order_alpha=NeoClassicalExp.names;
data_IHL=NeoClassicalExp.errors;
NeoClassicalExp.coordinates.errors= Textcoordinates(data_IHL,order_alpha);
%% Coordinates
%RE
data=NeoClassicalResults.REE_IRF_eg.y';
kn_ratio_re=NeoClassicalResults.REE_IRF_eg.y(:,1)-NeoClassicalResults.REE_IRF_eg.y(:,7);
data = [data; kn_ratio_re'];
order_alpha=NeoClassicalResults.REE_IRF_eg.names;
order_alpha{end+1}='kn';
NeoClassicalResults.coordinates.RE = Textcoordinates(data,order_alpha);

% IHL
data=NeoClassicalResults.IHL_eg.y;
kn_ratio_ihl=NeoClassicalResults.IHL_eg.y(1,:)-NeoClassicalResults.IHL_eg.y(7,:);
data = [data; kn_ratio_ihl];
order_alpha=NeoClassicalResults.IHL_eg.names_xy;
order_alpha{end+1}='kn';
NeoClassicalResults.coordinates.IHL= Textcoordinates(data,order_alpha);

% Expectations
order_alpha=NeoClassicalExp.names;
data_RL=NeoClassicalExp.yf_RL;
data_IHL=NeoClassicalExp.yf_IHL;
data_RE=NeoClassicalExp.yf_RE;
NeoClassicalExp.coordinates.IHL_exp= Textcoordinates(data_IHL,order_alpha);
NeoClassicalExp.coordinates.RE_exp= Textcoordinates(data_RE,order_alpha);
NeoClassicalExp.coordinates.RL_exp= Textcoordinates(data_RL,order_alpha);

%% New Keynesian model
NewKeynesianResults=NewKeynesian(parameters,irf_shock,sigmag,sigmaz,sigmar,T,S,draw,itermax,damp,constant,gain);

%% Expectations
selected_vars={'pi','r','rk','w','d'};
NewKeynesianExp=PlotExp(NewKeynesianResults,selected_vars);

% coordinates forecast errors
order_alpha=NewKeynesianExp.names;
data_IHL=NewKeynesianExp.errors;
NewKeynesianExp.coordinates.errors= Textcoordinates(data_IHL,order_alpha);
%% Coordinates
%RE
data=NewKeynesianResults.REE_IRF_eg.y';
order_alpha=NewKeynesianResults.REE_IRF_eg.names;
NewKeynesianResults.coordinates.RE = Textcoordinates(data,order_alpha);

% IHL
data=NewKeynesianResults.IHL_eg.y;
order_alpha=NewKeynesianResults.IHL_eg.names_xy;
NewKeynesianResults.coordinates.IHL= Textcoordinates(data,order_alpha);

% Expectations
order_alpha=NewKeynesianExp.names;
data_RL=NewKeynesianExp.yf_RL;
data_IHL=NewKeynesianExp.yf_IHL;
data_RE=NewKeynesianExp.yf_RE;
NewKeynesianExp.coordinates.IHL_exp= Textcoordinates(data_IHL,order_alpha);
NewKeynesianExp.coordinates.RE_exp= Textcoordinates(data_RE,order_alpha);
NewKeynesianExp.coordinates.RL_exp= Textcoordinates(data_RL,order_alpha);

%% New Keynesian model with distortionary taxation (*debt adjustment*)

rhob=0.1;

% Paper
tauk=0.3924853571;
tauw=0.3931357143;%0.39;
tauc=0.013076838;

% Leeper et al 2010
%tauk=0.184;
%tauw=0.223;%0.39;
%tauc=0.028;

% adjusted because of capital tax
rk_ss2=(beta^(-1) - 1 + delta)/(1-tauk);% steady state rental rate of captal

ky2=mc_ss*alpha/rk_ss2;
iy2=delta*ky2; % I/Y ratio ! different than without imperfect competition
cy2=1-gy-iy2;                 % C/Y ratio

% adjusted because of labour tax

phi2 = cy2*(1-tauw)^(-1)*(1+tauc)/(wy_ss*(1-n_ss)+cy2*(1+tauc)/(1-tauw)); % Preference parameter

dy2=1-wny-rk_ss2*ky2;

ty2=(1/beta-1)*by+tauc*(1-gy-iy2)+tauk*rk_ss2*ky2+tauw*wny-gy;
ty2=(1-1/beta)*by+tauc*(1-gy-iy2)+tauk*rk_ss2*ky2+tauw*wny-gy; % EQ correction oct 2018

%% Strategy 1
financing=0; % 0=lump-sum financing
parameters_dist=[alpha,beta,delta,eps,rhog,rhoz,sigma,gy,...
    rk_ss2,iy2,cy2,ky2,wny,dy2,n_ss,phi2,sigmai,theta,rhopi,rhor,by,tauk,tauw,tauc,...
    financing,rhob];
Debt0Results=DebtAdjustment(parameters_dist,irf_shock,sigmag,sigmaz,sigmar,T,S,draw,itermax,damp,constant,gain);
Debt0Exp=PlotExp(Debt0Results,selected_vars);

% coordinates forecast errors
Debt0Exp.coordinates.errors=...
    Textcoordinates(Debt0Exp.errors,Debt0Exp.names);

%% New Keynesian model with distortionary taxation (*lump-sum tax adjustment*)

%% Strategy 1
financing=0; % 0=lump-sum financing
parameters_dist=[alpha,beta,delta,eps,rhog,rhoz,sigma,gy,...
    rk_ss2,iy2,cy2,ky2,wny,dy2,n_ss,phi2,sigmai,theta,rhopi,rhor,by,tauk,tauw,tauc,financing];
Distortionary0Results=DistortionaryTaxes(parameters_dist,irf_shock,sigmag,sigmaz,sigmar,T,S,draw,itermax,damp,constant,gain);
Distortionary0Exp=PlotExp(Distortionary0Results,selected_vars);

% coordinates forecast errors
Distortionary0Exp.coordinates.errors=...
    Textcoordinates(Distortionary0Exp.errors,Distortionary0Exp.names);
%% Coordinates
%RE
data=Distortionary0Results.REE_IRF_eg.y';
order_alpha=Distortionary0Results.REE_IRF_eg.names;
Distortionary0Results.coordinates.RE = Textcoordinates(data,order_alpha);

% IHL
data=Distortionary0Results.IHL_eg.y;
order_alpha=Distortionary0Results.IHL_eg.names_xy;
Distortionary0Results.coordinates.IHL= Textcoordinates(data,order_alpha);

% Expectations
order_alpha=Distortionary0Exp.names;
data_RL=Distortionary0Exp.yf_RL;
data_IHL=Distortionary0Exp.yf_IHL;
data_RE=Distortionary0Exp.yf_RE;
Distortionary0Exp.coordinates.IHL_exp= Textcoordinates(data_IHL,order_alpha);
Distortionary0Exp.coordinates.RE_exp= Textcoordinates(data_RE,order_alpha);
Distortionary0Exp.coordinates.RL_exp= Textcoordinates(data_RL,order_alpha);

%% Strategy 2
financing=1; % 1=capital tax financing
parameters_dist(end)=financing;
Distortionary1Results=DistortionaryTaxes(parameters_dist,irf_shock,sigmag,sigmaz,sigmar,T,S,draw,itermax,damp,constant,gain);
Distortionary1Exp=PlotExp(Distortionary1Results,selected_vars);

% path of tauk
loc_g=find(strcmp(Distortionary1Results.REE_IRF_eg.names,'g'));
Distortionary1Results.tauk=gy*Distortionary1Results.REE_IRF_eg.y(:,loc_g)/(ky2*tauk*rk_ss2);

% coordinates forecast errors
Distortionary1Exp.coordinates.errors=...
    Textcoordinates(Distortionary1Exp.errors,Distortionary1Exp.names);
%% Coordinates
%RE
data=Distortionary1Results.REE_IRF_eg.y';
order_alpha=Distortionary1Results.REE_IRF_eg.names;
Distortionary1Results.coordinates.RE = Textcoordinates(data,order_alpha);

% tauk
name={'tauk'};
Distortionary1Results.coordinates.tauk=...
    Textcoordinates(Distortionary1Results.tauk',name);

% IHL
data=Distortionary1Results.IHL_eg.y;
order_alpha=Distortionary1Results.IHL_eg.names_xy;
Distortionary1Results.coordinates.IHL= Textcoordinates(data,order_alpha);

% Expectations
order_alpha=Distortionary1Exp.names;
data_RL=Distortionary1Exp.yf_RL;
data_IHL=Distortionary1Exp.yf_IHL;
data_RE=Distortionary1Exp.yf_RE;
Distortionary1Exp.coordinates.IHL_exp= Textcoordinates(data_IHL,order_alpha);
Distortionary1Exp.coordinates.RE_exp= Textcoordinates(data_RE,order_alpha);
Distortionary1Exp.coordinates.RL_exp= Textcoordinates(data_RL,order_alpha);


%% Strategy 3
financing=2; % 2=labour tax financing
parameters_dist(end)=financing;
Distortionary2Results=DistortionaryTaxes(parameters_dist,irf_shock,sigmag,sigmaz,sigmar,T,S,draw,itermax,damp,constant,gain);
Distortionary2Exp=PlotExp(Distortionary2Results,selected_vars);

% path of tauw
loc_g=find(strcmp(Distortionary2Results.REE_IRF_eg.names,'g'));
Distortionary2Results.tauw=gy*Distortionary2Results.REE_IRF_eg.y(:,loc_g)/(wny*tauw);

% coordinates forecast errors
Distortionary2Exp.coordinates.errors=...
    Textcoordinates(Distortionary2Exp.errors,Distortionary2Exp.names);
%% Coordinates
%RE
data=Distortionary2Results.REE_IRF_eg.y';
order_alpha=Distortionary2Results.REE_IRF_eg.names;
Distortionary2Results.coordinates.RE = Textcoordinates(data,order_alpha);

% tauw
name={'tauw'};
Distortionary2Results.coordinates.tauw=...
    Textcoordinates(Distortionary2Results.tauw',name);

% IHL
data=Distortionary2Results.IHL_eg.y;
order_alpha=Distortionary2Results.IHL_eg.names_xy;
Distortionary2Results.coordinates.IHL= Textcoordinates(data,order_alpha);

% Expectations
order_alpha=Distortionary2Exp.names;
data_RL=Distortionary2Exp.yf_RL;
data_IHL=Distortionary2Exp.yf_IHL;
data_RE=Distortionary2Exp.yf_RE;
Distortionary2Exp.coordinates.IHL_exp= Textcoordinates(data_IHL,order_alpha);
Distortionary2Exp.coordinates.RE_exp= Textcoordinates(data_RE,order_alpha);
Distortionary2Exp.coordinates.RL_exp= Textcoordinates(data_RL,order_alpha);


%% Combine coordinates

% prefix
name='y';
Series1.cor=NewKeynesianResults.coordinates;
Series2.cor=Distortionary1Results.coordinates;
Series3.cor=Distortionary2Results.coordinates;

%\addplot[color=red,solid,line width=0.75pt]

%% 6. Impulse responses to an increase in government spending for different degrees of non-separability
% This part of the code calculates the impulse responses for diffent
% degrees of the $$\gamma$$ parameter and plots the responses of output and
% private consumption (see figure 3 in the manuscript).
para_of_interest='rhog';
switch para_of_interest
    case 'sigma'
        interval=1:.5:4; % sigma
    case 'sigmai'
        interval=1:3:20; % sigmai
    case 'rhog'
        interval=.7:.05:.9; % rhog
    case 'rhopi'
        interval=1:.125:2; % rhopi
end
for j=1:size(NewKeynesianResults.REE_IRF_eg.names,1) % locates output variable in cell array of endogenous variables
    if strcmp(NewKeynesianResults.REE_IRF_eg.names(j),'y')            
        loc_output=j;
    end
end
for j=1:size(NewKeynesianResults.REE_IRF_eg.names,1) % locates consumption variable in cell array of endogenous variables
    if strcmp(NewKeynesianResults.REE_IRF_eg.names(j),'c')            
        loc_consumption=j;
    end
end
Output_RE=zeros(T,size(interval,2)); % impulse responses of output under rational expectations
Consumption_RE=zeros(T,size(interval,2)); % impulse responses of consumption under rational expectations
Output_IHL=zeros(T,size(interval,2)); % impulse responses of output under Infinite Horizon learning
Consumption_IHL=zeros(T,size(interval,2)); % impulse responses of consumption under Infinite Horizon learning
Output_RL=zeros(T,size(interval,2)); % impulse responses of output under Infinite Horizon learning
Consumption_RL=zeros(T,size(interval,2)); % impulse responses of consumption under Euler eq. learning
k=1;
switch para_of_interest
    case 'sigma'
        sigma_old=sigma; %initial value of sigma parameter is saved
    case 'sigmai'
        sigmai_old=sigmai; 
    case 'rhog'
        rhog_old=rhog;
     case 'rhopi'
        rhopi_old=rhopi; 
end

parameters_alt=parameters;
for i=interval
    switch para_of_interest
    case 'sigma'
        parameters_alt(7)=i; % sigma
    case 'sigmai'
        parameters_alt(17)=i; %sigmai
    case 'rhog'
        parameters_alt(5)=i;
    case 'rhopi'
        parameters_alt(19)=i;
    end
    
    NKResults_sigma=NewKeynesian(parameters_alt,irf_shock,sigmag,sigmaz,sigmar,T,S,draw,itermax,damp,constant,gain);
    
    % REE    
    Output_RE(:,k)=NKResults_sigma.REE_IRF_eg.y(:,loc_output);
    Consumption_RE(:,k)=NKResults_sigma.REE_IRF_eg.y(:,loc_consumption);
    
    % Infinite Horizon Learning
    Output_IHL(:,k)=NKResults_sigma.IHL_eg.y(loc_output,:)';
    Consumption_IHL(:,k)=NKResults_sigma.IHL_eg.y(loc_consumption,:)';
    
    % Euler Equation Learning
    Output_RL(:,k)=NKResults_sigma.RL_eg.y(loc_output,:)';
    Consumption_RL(:,k)=NKResults_sigma.RL_eg.y(loc_consumption,:)';    
    k=k+1;
end
switch para_of_interest
    case 'sigma'
        parameter ='$\sigma$';
    case 'sigmai'
        parameter ='$\varsigma_I$';
    case 'rhog'
        parameter ='$\rho_G$';
    case 'rhopi'
        parameter ='$\rho_\Pi$';
end

ticks = size(interval,2);% number of ticks on y axis
lim = size(interval,2);                % total number of parameter values considered
h = figure('Position',[100 100 1100 550],'Name','Impulse responses for different degrees of non-separability (IH Learning)');
subplot(1,2,1); hold on;grid on;surf(Output_IHL','FaceAlpha',.6,'FaceColor','blue'); surf(Output_RE','FaceAlpha',.6,'FaceColor','green'); hold off; set(gca,'YLim',[1 lim]); set(gca,'YTick',[1:ticks]); set(gca,'YTickLabel',interval);ylabel(parameter,'Interpreter','LaTex');view(20,30); title('Output');
subplot(1,2,2); hold on;grid on;surf(Consumption_IHL','FaceAlpha',.6,'FaceColor','blue');surf(Consumption_RE','FaceAlpha',.6,'FaceColor','green'); hold off; set(gca,'YLim',[1 lim]); set(gca,'YTick',[1:ticks]); set(gca,'YTickLabel',interval); ylabel(parameter,'Interpreter','LaTex'); view(20,30); title('Consumption');
print(h,'-dpng', 'Figure 3');
%%
h = figure('Name','Impulse responses for different degrees of non-separability (RL Learning)');
subplot(1,2,1); hold on;grid on;surf(Output_RL','FaceAlpha',.6,'FaceColor','blue'); surf(Output_RE','FaceAlpha',.6,'FaceColor','green'); hold off; set(gca,'YLim',[1 lim]); set(gca,'YTick',[1:ticks]); set(gca,'YTickLabel',interval);ylabel(parameter,'Interpreter','LaTex');view(20,30); title('Output');
subplot(1,2,2); hold on;grid on;surf(Consumption_RL','FaceAlpha',.6,'FaceColor','blue');surf(Consumption_RE','FaceAlpha',.6,'FaceColor','green'); hold off; set(gca,'YLim',[1 lim]); set(gca,'YTick',[1:ticks]); set(gca,'YTickLabel',interval); ylabel(parameter,'Interpreter','LaTex'); view(20,30); title('Consumption');
print(h,'-dpng', 'Figure 4');
h = figure('Name','Impulse responses for different degrees of non-separability (EE (blue) vs IH Learning)');
subplot(1,2,1); hold on;grid on;surf(Output_RL','FaceAlpha',.6,'FaceColor','blue'); surf(Output_IHL','FaceAlpha',.6,'FaceColor','green'); hold off; set(gca,'YLim',[1 lim]); set(gca,'YTick',[1:ticks]); set(gca,'YTickLabel',interval);ylabel(parameter,'Interpreter','LaTex');view(20,30); title('Output');
subplot(1,2,2); hold on;grid on;surf(Consumption_RL','FaceAlpha',.6,'FaceColor','blue');surf(Consumption_IHL','FaceAlpha',.6,'FaceColor','green'); hold off; set(gca,'YLim',[1 lim]); set(gca,'YTick',[1:ticks]); set(gca,'YTickLabel',interval); ylabel(parameter,'Interpreter','LaTex'); view(20,30); title('Consumption');
print(h,'-dpng', 'Figure 5');

%% Coordinates
data2.y=Output_RE';
data2.c=Consumption_RE';
order_alpha={'y','c'};
Figure5.RE = Textcoordinates3D(data2,order_alpha);
data3.y=Output_IHL';
data3.c=Consumption_IHL';
Figure5.IHL = Textcoordinates3D(data3,order_alpha);
%% 7. Impact multipliers for different degrees of price rigidity
% Replicates figure 4 of the manuscript: the impact multipliers of output,
% consumption, and investment for different degrees of price rigidity.

thetas=[0 .6 .75 .85]; % range of theta values considered
OutputMultiplier=zeros(size(thetas,2),2);
ConsumptionMultiplier=zeros(size(thetas,2),2);
InvestmentMultiplier=zeros(size(thetas,2),2);

%Neoclassical

%REE
OutputMultiplier(1,1)= 100*NeoClassicalResults.REE_IRF_eg.y(1,strcmp(NeoClassicalResults.REE_IRF_eg.names,'y'));
ConsumptionMultiplier(1,1)= 100*NeoClassicalResults.REE_IRF_eg.y(1,strcmp(NeoClassicalResults.REE_IRF_eg.names,'c'))*cy;
InvestmentMultiplier(1,1)= 100*NeoClassicalResults.REE_IRF_eg.y(1,strcmp(NeoClassicalResults.REE_IRF_eg.names,'i'))*iy;

% Infinite Horizon Learning
OutputMultiplier(1,2)= 100*NeoClassicalResults.IHL_eg.y(strcmp(NeoClassicalResults.IHL_eg.names_xy,'y'),1)';
ConsumptionMultiplier(1,2)= 100*NeoClassicalResults.IHL_eg.y(strcmp(NeoClassicalResults.IHL_eg.names_xy,'c'),1)'*cy;
InvestmentMultiplier(1,2)= 100*NeoClassicalResults.IHL_eg.y(strcmp(NeoClassicalResults.IHL_eg.names_xy,'i'),1)'*iy;

% Euler Equation Learning
OutputMultiplier(1,3)= 100*NeoClassicalResults.RL_eg.y(strcmp(NeoClassicalResults.RL_eg.names_xy,'y'),1)';
ConsumptionMultiplier(1,3)= 100*NeoClassicalResults.RL_eg.y(strcmp(NeoClassicalResults.RL_eg.names_xy,'c'),1)'*cy;
InvestmentMultiplier(1,3)= 100*NeoClassicalResults.RL_eg.y(strcmp(NeoClassicalResults.RL_eg.names_xy,'i'),1)'*iy;

parameters_alt2=parameters;
for i=2:size(thetas,2)
    parameters_alt2(18)=thetas(i);
    
    NKResults_theta=NewKeynesian(parameters_alt2,irf_shock,sigmag,sigmaz,sigmar,T,S,draw,itermax,damp,constant,gain);
    
    % REE
    OutputMultiplier(i,1)=100*NKResults_theta.REE_IRF_eg.y(1,strcmp(NKResults_theta.REE_IRF_eg.names,'y'));
    ConsumptionMultiplier(i,1)=100*NKResults_theta.REE_IRF_eg.y(1,strcmp(NKResults_theta.REE_IRF_eg.names,'c'))*cy;
    InvestmentMultiplier(i,1)=100*NKResults_theta.REE_IRF_eg.y(1,strcmp(NKResults_theta.REE_IRF_eg.names,'i'))*iy;
    
    % Infinite Horizon Learning
    OutputMultiplier(i,2)=100*NKResults_theta.IHL_eg.y(strcmp(NKResults_theta.IHL_eg.names_xy,'y'),1)';
    ConsumptionMultiplier(i,2)=100*NKResults_theta.IHL_eg.y(strcmp(NKResults_theta.IHL_eg.names_xy,'c'),1)'*cy;
    InvestmentMultiplier(i,2)=100*NKResults_theta.IHL_eg.y(strcmp(NKResults_theta.IHL_eg.names_xy,'i'),1)'*iy;
    
    % Euler Equation Learning
    OutputMultiplier(i,3)= 100*NKResults_theta.RL_eg.y(strcmp(NKResults_theta.RL_eg.names_xy,'y'),1)';
    ConsumptionMultiplier(i,3)= 100*NKResults_theta.RL_eg.y(strcmp(NKResults_theta.RL_eg.names_xy,'c'),1)'*cy;
    InvestmentMultiplier(i,3)= 100*NKResults_theta.RL_eg.y(strcmp(NKResults_theta.RL_eg.names_xy,'i'),1)'*iy;
end

h = figure('Name','Impact multipliers for different degrees of price rigidity');
subplot(2,2,1);bar(OutputMultiplier);set(gca,'XTickLabel',thetas)
title('Output');
x=[1:size(thetas,2)];
for i=1:numel(OutputMultiplier(1,:))
    text(x,OutputMultiplier(:,i),num2str(OutputMultiplier(:,i),'%0.2f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
subplot(2,2,2);bar(ConsumptionMultiplier);set(gca,'XTickLabel',thetas)
title('Consumption');
for i=1:numel(ConsumptionMultiplier(1,:))
    text(x,ConsumptionMultiplier(:,i),num2str(ConsumptionMultiplier(:,i),'%0.2f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
%colormap([.83 .83 .83; .04 .12 .38]);
subplot(2,2,3);bar(InvestmentMultiplier);set(gca,'XTickLabel',thetas)
title('Investment');
for i=1:numel(InvestmentMultiplier(1,:))
    text(x,InvestmentMultiplier(:,i),num2str(InvestmentMultiplier(:,i),'%0.2f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
hLegend = legend('Rational Expectations','IH Learning', 'EE Learning');
%set(hLegend,'Position', newPosition,'Units', newUnits);
print(h,'-dpng', 'Figure 6'); % figure is saved to file

%% X. Present-value fiscal multipliers

PVMultipliers=nan(T,45);

% Output
%  - Neoclassical
PVMultipliers(:,1) = PVMultiplier(NeoClassicalResults.REE_IRF_eg.y(:,strcmp(NeoClassicalResults.REE_IRF_eg.names,'y'))...
    ,NeoClassicalResults.REE_IRF_eg.y(:,strcmp(NeoClassicalResults.REE_IRF_eg.names,'g')),1/beta,gy,T);
PVMultipliers(:,2) = PVMultiplier(NeoClassicalResults.IHL_eg.y(strcmp(NeoClassicalResults.IHL_eg.names_xy,'y'),:)',...
    NeoClassicalResults.IHL_eg.w(strcmp(NeoClassicalResults.IHL_eg.names_z,'g'),:)',1/beta,gy,T);
PVMultipliers(:,3) = PVMultiplier(NeoClassicalResults.RL_eg.y(strcmp(NeoClassicalResults.RL_eg.names_xy,'y'),:)',...
    NeoClassicalResults.RL_eg.w(strcmp(NeoClassicalResults.RL_eg.names_z,'g'),:)',1/beta,gy,T);
%  - New Keynesian
PVMultipliers(:,4) = PVMultiplier(NewKeynesianResults.REE_IRF_eg.y(:,strcmp(NewKeynesianResults.REE_IRF_eg.names,'y')),...
    NewKeynesianResults.REE_IRF_eg.y(:,strcmp(NewKeynesianResults.REE_IRF_eg.names,'g')),1/beta,gy,T);
PVMultipliers(:,5) = PVMultiplier(NewKeynesianResults.IHL_eg.y(strcmp(NewKeynesianResults.IHL_eg.names_xy,'y'),:)',...
    NewKeynesianResults.IHL_eg.w(strcmp(NewKeynesianResults.IHL_eg.names_z,'g'),:)',1/beta,gy,T);
PVMultipliers(:,6) = PVMultiplier(NewKeynesianResults.RL_eg.y(strcmp(NewKeynesianResults.RL_eg.names_xy,'y'),:)',...
    NewKeynesianResults.RL_eg.w(strcmp(NewKeynesianResults.RL_eg.names_z,'g'),:)',1/beta,gy,T);
%  - Distortionary0
PVMultipliers(:,7) = PVMultiplier(Distortionary0Results.REE_IRF_eg.y(:,strcmp(Distortionary0Results.REE_IRF_eg.names,'y')),...
    Distortionary0Results.REE_IRF_eg.y(:,strcmp(Distortionary0Results.REE_IRF_eg.names,'g')),1/beta,gy,T);
PVMultipliers(:,8) = PVMultiplier(Distortionary0Results.IHL_eg.y(strcmp(Distortionary0Results.IHL_eg.names_xy,'y'),:)',...
    Distortionary0Results.IHL_eg.w(strcmp(Distortionary0Results.IHL_eg.names_z,'g'),:)',1/beta,gy,T);
PVMultipliers(:,9) = PVMultiplier(Distortionary0Results.RL_eg.y(strcmp(Distortionary0Results.RL_eg.names_xy,'y'),:)',...
    Distortionary0Results.RL_eg.w(strcmp(Distortionary0Results.RL_eg.names_z,'g'),:)',1/beta,gy,T);
%  - Distortionary1
PVMultipliers(:,10) = PVMultiplier(Distortionary1Results.REE_IRF_eg.y(:,strcmp(Distortionary1Results.REE_IRF_eg.names,'y')),...
    Distortionary1Results.REE_IRF_eg.y(:,strcmp(Distortionary1Results.REE_IRF_eg.names,'g')),1/beta,gy,T);
PVMultipliers(:,11) = PVMultiplier(Distortionary1Results.IHL_eg.y(strcmp(Distortionary1Results.IHL_eg.names_xy,'y'),:)',...
    Distortionary1Results.IHL_eg.w(strcmp(Distortionary1Results.IHL_eg.names_z,'g'),:)',1/beta,gy,T);
PVMultipliers(:,12) = PVMultiplier(Distortionary1Results.RL_eg.y(strcmp(Distortionary1Results.RL_eg.names_xy,'y'),:)',...
    Distortionary1Results.RL_eg.w(strcmp(Distortionary1Results.RL_eg.names_z,'g'),:)',1/beta,gy,T);
%  - Distortionary2
PVMultipliers(:,13) = PVMultiplier(Distortionary2Results.REE_IRF_eg.y(:,strcmp(Distortionary2Results.REE_IRF_eg.names,'y')),...
    Distortionary2Results.REE_IRF_eg.y(:,strcmp(Distortionary2Results.REE_IRF_eg.names,'g')),1/beta,gy,T);
PVMultipliers(:,14) = PVMultiplier(Distortionary2Results.IHL_eg.y(strcmp(Distortionary2Results.IHL_eg.names_xy,'y'),:)',...
    Distortionary2Results.IHL_eg.w(strcmp(Distortionary2Results.IHL_eg.names_z,'g'),:)',1/beta,gy,T);
PVMultipliers(:,15) = PVMultiplier(Distortionary2Results.RL_eg.y(strcmp(Distortionary2Results.RL_eg.names_xy,'y'),:)',...
    Distortionary2Results.RL_eg.w(strcmp(Distortionary2Results.RL_eg.names_z,'g'),:)',1/beta,gy,T);
% Consumption
%  - Neoclassical
PVMultipliers(:,16) = PVMultiplier(NeoClassicalResults.REE_IRF_eg.y(:,strcmp(NeoClassicalResults.REE_IRF_eg.names,'c')),...
    NeoClassicalResults.REE_IRF_eg.y(:,strcmp(NeoClassicalResults.REE_IRF_eg.names,'g')),1/beta,gy/cy,T);
PVMultipliers(:,17) = PVMultiplier(NeoClassicalResults.IHL_eg.y(strcmp(NeoClassicalResults.IHL_eg.names_xy,'c'),:)',...
    NeoClassicalResults.IHL_eg.w(strcmp(NeoClassicalResults.IHL_eg.names_z,'g'),:)',1/beta,gy/cy,T);
PVMultipliers(:,18) = PVMultiplier(NeoClassicalResults.RL_eg.y(strcmp(NeoClassicalResults.RL_eg.names_xy,'c'),:)',...
    NeoClassicalResults.RL_eg.w(strcmp(NeoClassicalResults.RL_eg.names_z,'g'),:)',1/beta,gy/cy,T);
%  - New Keynesian
PVMultipliers(:,19) = PVMultiplier(NewKeynesianResults.REE_IRF_eg.y(:,strcmp(NewKeynesianResults.REE_IRF_eg.names,'c')),...
    NewKeynesianResults.REE_IRF_eg.y(:,strcmp(NewKeynesianResults.REE_IRF_eg.names,'g')),1/beta,gy/cy,T);
PVMultipliers(:,20) = PVMultiplier(NewKeynesianResults.IHL_eg.y(strcmp(NewKeynesianResults.IHL_eg.names_xy,'c'),:)',...
    NewKeynesianResults.IHL_eg.w(strcmp(NewKeynesianResults.IHL_eg.names_z,'g'),:)',1/beta,gy/cy,T);
PVMultipliers(:,21) = PVMultiplier(NewKeynesianResults.RL_eg.y(strcmp(NewKeynesianResults.RL_eg.names_xy,'c'),:)',...
    NewKeynesianResults.RL_eg.w(strcmp(NewKeynesianResults.RL_eg.names_z,'g'),:)',1/beta,gy/cy,T);
%  - Distortionary0
PVMultipliers(:,22) = PVMultiplier(Distortionary0Results.REE_IRF_eg.y(:,strcmp(Distortionary0Results.REE_IRF_eg.names,'c')),...
    Distortionary0Results.REE_IRF_eg.y(:,strcmp(Distortionary0Results.REE_IRF_eg.names,'g')),1/beta,gy/cy,T);
PVMultipliers(:,23) = PVMultiplier(Distortionary0Results.IHL_eg.y(strcmp(Distortionary0Results.IHL_eg.names_xy,'c'),:)',...
    Distortionary0Results.IHL_eg.w(strcmp(Distortionary0Results.IHL_eg.names_z,'g'),:)',1/beta,gy/cy,T);
PVMultipliers(:,24) = PVMultiplier(Distortionary0Results.RL_eg.y(strcmp(Distortionary0Results.RL_eg.names_xy,'c'),:)',...
    Distortionary0Results.RL_eg.w(strcmp(Distortionary0Results.RL_eg.names_z,'g'),:)',1/beta,gy/cy,T);
%  - Distortionary1
PVMultipliers(:,25) = PVMultiplier(Distortionary1Results.REE_IRF_eg.y(:,strcmp(Distortionary1Results.REE_IRF_eg.names,'c')),...
    Distortionary1Results.REE_IRF_eg.y(:,strcmp(Distortionary1Results.REE_IRF_eg.names,'g')),1/beta,gy/cy,T);
PVMultipliers(:,26) = PVMultiplier(Distortionary1Results.IHL_eg.y(strcmp(Distortionary1Results.IHL_eg.names_xy,'c'),:)',...
    Distortionary1Results.IHL_eg.w(strcmp(Distortionary1Results.IHL_eg.names_z,'g'),:)',1/beta,gy/cy,T);
PVMultipliers(:,27) = PVMultiplier(Distortionary1Results.RL_eg.y(strcmp(Distortionary1Results.RL_eg.names_xy,'c'),:)',...
    Distortionary1Results.RL_eg.w(strcmp(Distortionary1Results.RL_eg.names_z,'g'),:)',1/beta,gy/cy,T);
%  - Distortionary2
PVMultipliers(:,28) = PVMultiplier(Distortionary2Results.REE_IRF_eg.y(:,strcmp(Distortionary2Results.REE_IRF_eg.names,'c')),...
    Distortionary2Results.REE_IRF_eg.y(:,strcmp(Distortionary2Results.REE_IRF_eg.names,'g')),1/beta,gy/cy,T);
PVMultipliers(:,29) = PVMultiplier(Distortionary2Results.IHL_eg.y(strcmp(Distortionary2Results.IHL_eg.names_xy,'c'),:)',...
    Distortionary2Results.IHL_eg.w(strcmp(Distortionary2Results.IHL_eg.names_z,'g'),:)',1/beta,gy/cy,T);
PVMultipliers(:,30) = PVMultiplier(Distortionary2Results.RL_eg.y(strcmp(Distortionary2Results.RL_eg.names_xy,'c'),:)',...
    Distortionary2Results.RL_eg.w(strcmp(Distortionary2Results.RL_eg.names_z,'g'),:)',1/beta,gy/cy,T);
% Investment
%  - Neoclassical
PVMultipliers(:,31) = PVMultiplier(NeoClassicalResults.REE_IRF_eg.y(:,strcmp(NeoClassicalResults.REE_IRF_eg.names,'i')),...
    NeoClassicalResults.REE_IRF_eg.y(:,strcmp(NeoClassicalResults.REE_IRF_eg.names,'g')),1/beta,gy/iy,T);
PVMultipliers(:,32) = PVMultiplier(NeoClassicalResults.IHL_eg.y(strcmp(NeoClassicalResults.IHL_eg.names_xy,'i'),:)',...
    NeoClassicalResults.IHL_eg.w(strcmp(NeoClassicalResults.IHL_eg.names_z,'g'),:)',1/beta,gy/iy,T);
PVMultipliers(:,33) = PVMultiplier(NeoClassicalResults.RL_eg.y(strcmp(NeoClassicalResults.RL_eg.names_xy,'i'),:)',...
    NeoClassicalResults.RL_eg.w(strcmp(NeoClassicalResults.RL_eg.names_z,'g'),:)',1/beta,gy/iy,T);
%  - New Keynesian
PVMultipliers(:,34) = PVMultiplier(NewKeynesianResults.REE_IRF_eg.y(:,strcmp(NewKeynesianResults.REE_IRF_eg.names,'i')),...
    NewKeynesianResults.REE_IRF_eg.y(:,strcmp(NewKeynesianResults.REE_IRF_eg.names,'g')),1/beta,gy/iy,T);
PVMultipliers(:,35) = PVMultiplier(NewKeynesianResults.IHL_eg.y(strcmp(NewKeynesianResults.IHL_eg.names_xy,'i'),:)',...
    NewKeynesianResults.IHL_eg.w(strcmp(NewKeynesianResults.IHL_eg.names_z,'g'),:)',1/beta,gy/iy,T);
PVMultipliers(:,36) = PVMultiplier(NewKeynesianResults.RL_eg.y(strcmp(NewKeynesianResults.RL_eg.names_xy,'i'),:)',...
    NewKeynesianResults.RL_eg.w(strcmp(NewKeynesianResults.RL_eg.names_z,'g'),:)',1/beta,gy/iy,T);
%  - Distortionary0
PVMultipliers(:,37) = PVMultiplier(Distortionary0Results.REE_IRF_eg.y(:,strcmp(Distortionary0Results.REE_IRF_eg.names,'i')),...
    Distortionary0Results.REE_IRF_eg.y(:,strcmp(Distortionary0Results.REE_IRF_eg.names,'g')),1/beta,gy/iy,T);
PVMultipliers(:,38) = PVMultiplier(Distortionary0Results.IHL_eg.y(strcmp(Distortionary0Results.IHL_eg.names_xy,'i'),:)',...
    Distortionary0Results.IHL_eg.w(strcmp(Distortionary0Results.IHL_eg.names_z,'g'),:)',1/beta,gy/iy,T);
PVMultipliers(:,39) = PVMultiplier(Distortionary0Results.RL_eg.y(strcmp(Distortionary0Results.RL_eg.names_xy,'i'),:)',...
    Distortionary0Results.RL_eg.w(strcmp(Distortionary0Results.RL_eg.names_z,'g'),:)',1/beta,gy/iy,T);
%  - Distortionary1
PVMultipliers(:,40) = PVMultiplier(Distortionary1Results.REE_IRF_eg.y(:,strcmp(Distortionary1Results.REE_IRF_eg.names,'i')),...
    Distortionary1Results.REE_IRF_eg.y(:,strcmp(Distortionary1Results.REE_IRF_eg.names,'g')),1/beta,gy/iy,T);
PVMultipliers(:,41) = PVMultiplier(Distortionary1Results.IHL_eg.y(strcmp(Distortionary1Results.IHL_eg.names_xy,'i'),:)',...
    Distortionary1Results.IHL_eg.w(strcmp(Distortionary1Results.IHL_eg.names_z,'g'),:)',1/beta,gy/iy,T);
PVMultipliers(:,42) = PVMultiplier(Distortionary1Results.RL_eg.y(strcmp(Distortionary1Results.RL_eg.names_xy,'i'),:)',...
    Distortionary1Results.RL_eg.w(strcmp(Distortionary1Results.RL_eg.names_z,'g'),:)',1/beta,gy/iy,T);
%  - Distortionary2
PVMultipliers(:,43) = PVMultiplier(Distortionary2Results.REE_IRF_eg.y(:,strcmp(Distortionary2Results.REE_IRF_eg.names,'i')),...
    Distortionary2Results.REE_IRF_eg.y(:,strcmp(Distortionary2Results.REE_IRF_eg.names,'g')),1/beta,gy/iy,T);
PVMultipliers(:,44) = PVMultiplier(Distortionary2Results.IHL_eg.y(strcmp(Distortionary2Results.IHL_eg.names_xy,'i'),:)',...
    Distortionary2Results.IHL_eg.w(strcmp(Distortionary2Results.IHL_eg.names_z,'g'),:)',1/beta,gy/iy,T);
PVMultipliers(:,45) = PVMultiplier(Distortionary2Results.RL_eg.y(strcmp(Distortionary2Results.RL_eg.names_xy,'i'),:)',...
    Distortionary2Results.RL_eg.w(strcmp(Distortionary2Results.RL_eg.names_z,'g'),:)',1/beta,gy/iy,T);

%Table4=[PVMultipliers(1,4) PVMultipliers(4,4) PVMultipliers(16,4) PVMultipliers(24,4)...
%    PVMultipliers(1,5) PVMultipliers(4,5) PVMultipliers(16,5) PVMultipliers(24,5)...
%    PVMultipliers(1,6) PVMultipliers(4,6) PVMultipliers(16,6) PVMultipliers(24,6);
%    PVMultipliers(1,10) PVMultipliers(4,10) PVMultipliers(16,10) PVMultipliers(24,10)...
%    PVMultipliers(1,11) PVMultipliers(4,11) PVMultipliers(16,11) PVMultipliers(24,11)...
%    PVMultipliers(1,12) PVMultipliers(4,12) PVMultipliers(16,12) PVMultipliers(24,12);
%    PVMultipliers(1,16) PVMultipliers(4,16) PVMultipliers(16,16) PVMultipliers(24,16)...
%    PVMultipliers(1,17) PVMultipliers(4,17) PVMultipliers(16,17) PVMultipliers(24,17)...
%    PVMultipliers(1,18) PVMultipliers(4,18) PVMultipliers(16,18) PVMultipliers(24,18)];