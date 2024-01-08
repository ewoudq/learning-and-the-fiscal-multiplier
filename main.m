%% Learning and the size of the Government Spending Multiplier
% © 2017 Ewoud Quaghebeur

%% Parameters
alpha = 1/3;                % Output elasticity with respect to capital
beta = 1.04^(-0.25);        % Households subjective discount factor
gain = .02;                 % Gain parameter
delta = 0.025;              % Rate of physical capital depreciation
eps=6;                      % Elasticity of substitution between intermediate goods
theta=.75;                  % Degree of nominal price rigidity
rhog = .9;                  % Government expenditure AR(1) coefficient
rhopi = 1.5;                % Taylor rule inflation rate coefficient
rhor = 0.5;                 % Interest rate AR(1) coefficient
rhoz = .9;                  % Technology shock AR(1) coefficient
sigma = 2;%2;               % Coefficient of risk aversion
sigmag = .05;               % Standard deviation of government spending shock.  Units: Percent.
sigmar = .05;               % Standard deviation of interest rate shock. Units: Percent.
sigmaz = .05;               % Standard deviation of technology shock. Units: Percent.
sigmai = 17;                % Capital adjustment cost parameter
gy = .2;                    % Steady state government expenditure to output ratio

rk_ss = (beta^(-1) - 1 + delta);% steady state rental rate of captal

mc_ss = (eps-1)/eps;          % steady state marginal cost
ky = mc_ss*alpha/rk_ss;
iy = delta*ky; % I/Y ratio ! different than without imperfect competition
cy = 1 - gy - iy;             % C/Y ratio
n_ss = 1/3;                   % Steady state labor
wy_ss = ((1-alpha)*(eps-1))/(eps*n_ss);  % steady-state wage rate
by = 0.74;                  % Steady state debt to gdp ratio
phi = cy/(wy_ss*(1-n_ss)+cy); % Preference parameter

wny = (1-alpha)*mc_ss;
dy = 1 - wny - rk_ss*ky;

ty = gy + (1/beta-1)*by;

parameters = [alpha,beta,delta,eps,rhog,rhoz,sigma,gy,...
    rk_ss,iy,cy,ky,wny,dy,n_ss,phi,sigmai,theta,rhopi,rhor,by];

%% Learning parameters
constant = 0;
gain = 0.02;

%% Programme settings
draw = 0;       % Switch to draw additional figures (0 = no figures)
T = 40;         % Number of periods for impulse response plots
S = 10000;      % Number of simulation periods
itermax = 10000;% RPE: maximum number of iterations
damp = .9;      % RPE: damping factor

%% Neoclassical model
NeoClassicalResults = NeoClassical(parameters,{'eg'},sigmag,sigmaz,T,S,draw,itermax,damp,constant,gain);

% Figure 1. Impulse responses to an increase in government spending of 1%
% of GDP in the neoclassical specification of the model.
PlotIRFs(NeoClassicalResults,'Neoclassical');

%% New Keynesian model
NewKeynesianResults = NewKeynesian(parameters,{'eg'},sigmag,sigmaz,sigmar,T,S,draw,itermax,damp,constant,gain);

% Figure 2. Impulse responses to an increase in government spending of 1%
% of GDP in the new Keynesian model.
PlotIRFs(NewKeynesianResults,'New Keynesian');

%% Figure 3. Impulse responses to an increase in government spending for different degrees of non-separability
% Calculate the impulse responses for diffent degrees of the $$\sigma$$
% parameter and plot the responses of output and private consumption.

interval = 1:.5:4; % sigma

Output_RE = zeros(T,size(interval,2)); % impulse responses of output under rational expectations
Consumption_RE = zeros(T,size(interval,2)); % impulse responses of consumption under rational expectations
Output_IHL = zeros(T,size(interval,2)); % impulse responses of output under Infinite Horizon learning
Consumption_IHL = zeros(T,size(interval,2)); % impulse responses of consumption under Infinite Horizon learning
k = 1;
sigma_old = sigma; %initial value of sigma parameter is saved

parameters_alt = parameters;

for i = interval
    parameters_alt(7) = i; % sigma
    NKResults_sigma = NewKeynesian(parameters_alt,{'eg'},sigmag,sigmaz,sigmar,T,S,draw,itermax,damp,constant,gain);

    % REE
    Output_RE(:,k) = 100*NKResults_sigma.REE_IRF_eg.y(:,strcmp(NewKeynesianResults.REE_IRF_eg.names,'y'));
    Consumption_RE(:,k) = 100*NKResults_sigma.REE_IRF_eg.y(:,strcmp(NewKeynesianResults.REE_IRF_eg.names,'c'));
    
    % Infinite Horizon Learning
    Output_IHL(:,k) = 100*NKResults_sigma.IHL_eg.y(strcmp(NewKeynesianResults.REE_IRF_eg.names,'y'),:)';
    Consumption_IHL(:,k) = 100*NKResults_sigma.IHL_eg.y(strcmp(NewKeynesianResults.REE_IRF_eg.names,'c'),:)'; 
    k = k + 1;
end

parameter = '$\sigma$';
ticks = size(interval,2);% number of ticks on y axis
lim = size(interval,2);                % total number of parameter values considered
figure('Name','Impulse responses for different degrees of non-separability');
subplot(1,2,1); hold on;grid on;surf(Output_IHL','FaceAlpha',.6,'FaceColor','blue'); surf(Output_RE','FaceAlpha',.6,'FaceColor','green'); hold off; set(gca,'YLim',[1 lim]); set(gca,'YTick',[1:ticks]); set(gca,'YTickLabel',interval);ylabel(parameter,'Interpreter','LaTex');view(20,30); title('Output');
subplot(1,2,2); hold on;grid on;surf(Consumption_IHL','FaceAlpha',.6,'FaceColor','blue');surf(Consumption_RE','FaceAlpha',.6,'FaceColor','green'); hold off; set(gca,'YLim',[1 lim]); set(gca,'YTick',[1:ticks]); set(gca,'YTickLabel',interval); ylabel(parameter,'Interpreter','LaTex'); view(20,30); title('Consumption');
legend('Adaptive learning','Rational Expectations','Location','northwest');

%% Table 3. Present-value multipliers in the new Keynesian model

% Output
RE_Multipliers = PVMultiplier(NewKeynesianResults.REE_IRF_eg.y(:,strcmp(NewKeynesianResults.REE_IRF_eg.names,'y')),...
    NewKeynesianResults.REE_IRF_eg.y(:,strcmp(NewKeynesianResults.REE_IRF_eg.names,'g')),1/beta,gy,T);
Table3_PVMultipliers = RE_Multipliers([1 4 4*4 4*6])';
IHL_Multipliers = PVMultiplier(NewKeynesianResults.IHL_eg.y(strcmp(NewKeynesianResults.IHL_eg.names_xy,'y'),:)',...
    NewKeynesianResults.IHL_eg.w(strcmp(NewKeynesianResults.IHL_eg.names_z,'g'),:)',1/beta,gy,T);
Table3_PVMultipliers = [Table3_PVMultipliers IHL_Multipliers([1 4 4*4 4*6])'];

% Consumption
RE_Multipliers = PVMultiplier(NewKeynesianResults.REE_IRF_eg.y(:,strcmp(NewKeynesianResults.REE_IRF_eg.names,'c')),...
    NewKeynesianResults.REE_IRF_eg.y(:,strcmp(NewKeynesianResults.REE_IRF_eg.names,'g')),1/beta,gy/cy,T);

IHL_Multipliers = PVMultiplier(NewKeynesianResults.IHL_eg.y(strcmp(NewKeynesianResults.IHL_eg.names_xy,'c'),:)',...
    NewKeynesianResults.IHL_eg.w(strcmp(NewKeynesianResults.IHL_eg.names_z,'g'),:)',1/beta,gy/cy,T);
Table3_PVMultipliers = [Table3_PVMultipliers;RE_Multipliers([1 4 4*4 4*6])' IHL_Multipliers([1 4 4*4 4*6])'];

% Investment
RE_Multipliers = PVMultiplier(NewKeynesianResults.REE_IRF_eg.y(:,strcmp(NewKeynesianResults.REE_IRF_eg.names,'i')),...
    NewKeynesianResults.REE_IRF_eg.y(:,strcmp(NewKeynesianResults.REE_IRF_eg.names,'g')),1/beta,gy/iy,T);
IHL_Multipliers = PVMultiplier(NewKeynesianResults.IHL_eg.y(strcmp(NewKeynesianResults.IHL_eg.names_xy,'i'),:)',...
    NewKeynesianResults.IHL_eg.w(strcmp(NewKeynesianResults.IHL_eg.names_z,'g'),:)',1/beta,gy/iy,T);
Table3_PVMultipliers = [Table3_PVMultipliers; RE_Multipliers([1 4 4*4 4*6])' IHL_Multipliers([1 4 4*4 4*6])'];

clear RE_Multipliers IHL_Multipliers;
Table3_tmp = Table3_PVMultipliers;
Table3_PVMultipliers = array2table(Table3_PVMultipliers,'RowNames',{'Output','Consumption','Investment'}, ...
    'VariableNames',{'RE Impact','RE 1 yr','RE 4 yrs','RE 6 yrs',...
    'AL Impact','AL 1 yr','AL 4 yrs','AL 6 yrs'});

%% Figure 4. Impact multipliers for different degrees of price rigidity
% Impact multipliers of output, consumption, and investment for different
% degrees of price rigidity.

thetas = [0 .6 .75 .85]; % range of theta values considered
Figure4OutputMultiplier = zeros(size(thetas,2),2);
Figure4ConsumptionMultiplier = zeros(size(thetas,2),2);
Figure4InvestmentMultiplier = zeros(size(thetas,2),2);

%%% Neoclassical (no price rigidity)

%REE
Figure4OutputMultiplier(1,1) = 100*NeoClassicalResults.REE_IRF_eg.y(1,strcmp(NeoClassicalResults.REE_IRF_eg.names,'y'));
Figure4ConsumptionMultiplier(1,1) = 100*NeoClassicalResults.REE_IRF_eg.y(1,strcmp(NeoClassicalResults.REE_IRF_eg.names,'c'))*cy;
Figure4InvestmentMultiplier(1,1) = 100*NeoClassicalResults.REE_IRF_eg.y(1,strcmp(NeoClassicalResults.REE_IRF_eg.names,'i'))*iy;

% Infinite Horizon Learning
Figure4OutputMultiplier(1,2) = 100*NeoClassicalResults.IHL_eg.y(strcmp(NeoClassicalResults.IHL_eg.names_xy,'y'),1)';
Figure4ConsumptionMultiplier(1,2) = 100*NeoClassicalResults.IHL_eg.y(strcmp(NeoClassicalResults.IHL_eg.names_xy,'c'),1)'*cy;
Figure4InvestmentMultiplier(1,2) = 100*NeoClassicalResults.IHL_eg.y(strcmp(NeoClassicalResults.IHL_eg.names_xy,'i'),1)'*iy;

%%% New Keynesian
parameters_alt2 = parameters;
for i = 2:size(thetas,2)
    parameters_alt2(18) = thetas(i);
    
    NKResults_theta = NewKeynesian(parameters_alt2,{'eg'},sigmag,sigmaz,sigmar,T,S,draw,itermax,damp,constant,gain);
    
    % REE
    Figure4OutputMultiplier(i,1) = 100*NKResults_theta.REE_IRF_eg.y(1,strcmp(NKResults_theta.REE_IRF_eg.names,'y'));
    Figure4ConsumptionMultiplier(i,1) = 100*NKResults_theta.REE_IRF_eg.y(1,strcmp(NKResults_theta.REE_IRF_eg.names,'c'))*cy;
    Figure4InvestmentMultiplier(i,1) = 100*NKResults_theta.REE_IRF_eg.y(1,strcmp(NKResults_theta.REE_IRF_eg.names,'i'))*iy;
    
    % Infinite Horizon Learning
    Figure4OutputMultiplier(i,2) = 100*NKResults_theta.IHL_eg.y(strcmp(NKResults_theta.IHL_eg.names_xy,'y'),1)';
    Figure4ConsumptionMultiplier(i,2) = 100*NKResults_theta.IHL_eg.y(strcmp(NKResults_theta.IHL_eg.names_xy,'c'),1)'*cy;
    Figure4InvestmentMultiplier(i,2) = 100*NKResults_theta.IHL_eg.y(strcmp(NKResults_theta.IHL_eg.names_xy,'i'),1)'*iy; 
end

figure('Name','Impact multipliers for different degrees of price rigidity');
subplot(2,2,1);bar(Figure4OutputMultiplier);set(gca,'XTickLabel',thetas)
title('Output');
x = 1:size(thetas,2);
for i = 1:numel(Figure4OutputMultiplier(1,:))
    text(x,Figure4OutputMultiplier(:,i),num2str(Figure4OutputMultiplier(:,i),'%0.2f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
subplot(2,2,2);bar(Figure4ConsumptionMultiplier);set(gca,'XTickLabel',thetas)
title('Consumption');
for i = 1:numel(Figure4ConsumptionMultiplier(1,:))
    text(x,Figure4ConsumptionMultiplier(:,i),num2str(Figure4ConsumptionMultiplier(:,i),'%0.2f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
subplot(2,2,3);bar(Figure4InvestmentMultiplier);set(gca,'XTickLabel',thetas)
title('Investment');
for i = 1:numel(Figure4InvestmentMultiplier(1,:))
    text(x,Figure4InvestmentMultiplier(:,i),num2str(Figure4InvestmentMultiplier(:,i),'%0.2f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
hLegend = legend('Rational Expectations','Adaptive Learning','Location','southeast');

%%% Create tables
RowNames = {['theta=' num2str(thetas(1))],['theta=' num2str(thetas(2))],...
    ['theta=' num2str(thetas(3))],['theta=' num2str(thetas(4))]};
VariableNames = {'RatExp','Learning'};
Figure4OutputMultiplier = array2table(Figure4OutputMultiplier,...
    "RowNames",RowNames,"VariableNames",VariableNames);
Figure4ConsumptionMultiplier = array2table(Figure4ConsumptionMultiplier,...
    "RowNames",RowNames,"VariableNames",VariableNames);
Figure4InvestmentMultiplier = array2table(Figure4InvestmentMultiplier,...
    "RowNames",RowNames,"VariableNames",VariableNames);

%% Alternative specification of fiscal policy (distortionary taxes) 

% Distortionary tax rates
tauk = 0.3924853571;
tauw = 0.3931357143;
tauc = 0.013076838;

% Steady state rental rate of captal, adjusted because of capital tax
rk_ss2 = (beta^(-1) - 1 + delta)/(1-tauk);% 

ky2 = mc_ss*alpha/rk_ss2;
iy2 = delta*ky2; % I/Y ratio
cy2 = 1-gy-iy2;  % C/Y ratio

% Preference parameter, adjusted because of labor tax
phi2 = cy2*(1-tauw)^(-1)*(1+tauc)/(wy_ss*(1-n_ss)+cy2*(1+tauc)/(1-tauw));

dy2 = 1-wny-rk_ss2*ky2;
ty2 = (1-1/beta)*by+tauc*(1-gy-iy2)+tauk*rk_ss2*ky2+tauw*wny-gy;

%%% Capital tax financing
financing = 1;
parameters_dist = [alpha,beta,delta,eps,rhog,rhoz,sigma,gy,...
    rk_ss2,iy2,cy2,ky2,wny,dy2,n_ss,phi2,sigmai,theta,rhopi,rhor,by,tauk,tauw,tauc,financing];
Distortionary1Results = DistortionaryTaxes(parameters_dist,{'eg'},sigmag,sigmaz,sigmar,T,S,draw,itermax,damp,constant,gain);

% path of tauk
loc_g = find(strcmp(Distortionary1Results.REE_IRF_eg.names,'g'));
Distortionary1Results.tauk = gy*Distortionary1Results.REE_IRF_eg.y(:,loc_g)/(ky2*tauk*rk_ss2);

%%% Labor tax financing
financing = 2; % 2=labor tax financing
parameters_dist(end) = financing;
Distortionary2Results = DistortionaryTaxes(parameters_dist,{'eg'},sigmag,sigmaz,sigmar,T,S,draw,itermax,damp,constant,gain);

% path of tauw
loc_g = find(strcmp(Distortionary2Results.REE_IRF_eg.names,'g'));
Distortionary2Results.tauw = gy*Distortionary2Results.REE_IRF_eg.y(:,loc_g)/(wny*tauw);

%% Table 4. Present-value multipliers for different specifications of fiscal policy
% Present-value multipliers for different specifications of fiscal policy
% in the new Keynesian model with rational expectations and with adaptive
% learning

%%% Lump-sum financing (baseline model)
Table4_PVMultipliers = Table3_tmp;

%%% Capital tax financing
% - Output
RE_Multipliers = PVMultiplier(Distortionary1Results.REE_IRF_eg.y(:,strcmp(Distortionary1Results.REE_IRF_eg.names,'y')),...
    Distortionary1Results.REE_IRF_eg.y(:,strcmp(Distortionary1Results.REE_IRF_eg.names,'g')),1/beta,gy,T);
IHL_Multipliers = PVMultiplier(Distortionary1Results.IHL_eg.y(strcmp(Distortionary1Results.IHL_eg.names_xy,'y'),:)',...
    Distortionary1Results.IHL_eg.w(strcmp(Distortionary1Results.IHL_eg.names_z,'g'),:)',1/beta,gy,T);
Table4_PVMultipliers = [Table4_PVMultipliers; RE_Multipliers([1 4 4*4 4*6])' IHL_Multipliers([1 4 4*4 4*6])'];
% - Consumption
RE_Multipliers = PVMultiplier(Distortionary1Results.REE_IRF_eg.y(:,strcmp(Distortionary1Results.REE_IRF_eg.names,'c')),...
    Distortionary1Results.REE_IRF_eg.y(:,strcmp(Distortionary1Results.REE_IRF_eg.names,'g')),1/beta,gy/cy,T);
IHL_Multipliers = PVMultiplier(Distortionary1Results.IHL_eg.y(strcmp(Distortionary1Results.IHL_eg.names_xy,'c'),:)',...
    Distortionary1Results.IHL_eg.w(strcmp(Distortionary1Results.IHL_eg.names_z,'g'),:)',1/beta,gy/cy,T);
Table4_PVMultipliers = [Table4_PVMultipliers; RE_Multipliers([1 4 4*4 4*6])' IHL_Multipliers([1 4 4*4 4*6])'];
% - Investment
RE_Multipliers = PVMultiplier(Distortionary1Results.REE_IRF_eg.y(:,strcmp(Distortionary1Results.REE_IRF_eg.names,'i')),...
    Distortionary1Results.REE_IRF_eg.y(:,strcmp(Distortionary1Results.REE_IRF_eg.names,'g')),1/beta,gy/iy,T);
IHL_Multipliers = PVMultiplier(Distortionary1Results.IHL_eg.y(strcmp(Distortionary1Results.IHL_eg.names_xy,'i'),:)',...
    Distortionary1Results.IHL_eg.w(strcmp(Distortionary1Results.IHL_eg.names_z,'g'),:)',1/beta,gy/iy,T);
Table4_PVMultipliers = [Table4_PVMultipliers; RE_Multipliers([1 4 4*4 4*6])' IHL_Multipliers([1 4 4*4 4*6])'];

%%% Labor tax financing
%  - Output
RE_Multipliers = PVMultiplier(Distortionary2Results.REE_IRF_eg.y(:,strcmp(Distortionary2Results.REE_IRF_eg.names,'y')),...
    Distortionary2Results.REE_IRF_eg.y(:,strcmp(Distortionary2Results.REE_IRF_eg.names,'g')),1/beta,gy,T);
IHL_Multipliers = PVMultiplier(Distortionary2Results.IHL_eg.y(strcmp(Distortionary2Results.IHL_eg.names_xy,'y'),:)',...
    Distortionary2Results.IHL_eg.w(strcmp(Distortionary2Results.IHL_eg.names_z,'g'),:)',1/beta,gy,T);
Table4_PVMultipliers = [Table4_PVMultipliers; RE_Multipliers([1 4 4*4 4*6])' IHL_Multipliers([1 4 4*4 4*6])'];
%  - Consumption
RE_Multipliers = PVMultiplier(Distortionary2Results.REE_IRF_eg.y(:,strcmp(Distortionary2Results.REE_IRF_eg.names,'c')),...
    Distortionary2Results.REE_IRF_eg.y(:,strcmp(Distortionary2Results.REE_IRF_eg.names,'g')),1/beta,gy/cy,T);
IHL_Multipliers = PVMultiplier(Distortionary2Results.IHL_eg.y(strcmp(Distortionary2Results.IHL_eg.names_xy,'c'),:)',...
    Distortionary2Results.IHL_eg.w(strcmp(Distortionary2Results.IHL_eg.names_z,'g'),:)',1/beta,gy/cy,T);
Table4_PVMultipliers = [Table4_PVMultipliers; RE_Multipliers([1 4 4*4 4*6])' IHL_Multipliers([1 4 4*4 4*6])'];
%  - Investment
RE_Multipliers = PVMultiplier(Distortionary2Results.REE_IRF_eg.y(:,strcmp(Distortionary2Results.REE_IRF_eg.names,'i')),...
    Distortionary2Results.REE_IRF_eg.y(:,strcmp(Distortionary2Results.REE_IRF_eg.names,'g')),1/beta,gy/iy,T);
IHL_Multipliers = PVMultiplier(Distortionary2Results.IHL_eg.y(strcmp(Distortionary2Results.IHL_eg.names_xy,'i'),:)',...
    Distortionary2Results.IHL_eg.w(strcmp(Distortionary2Results.IHL_eg.names_z,'g'),:)',1/beta,gy/iy,T);
Table4_PVMultipliers = [Table4_PVMultipliers; RE_Multipliers([1 4 4*4 4*6])' IHL_Multipliers([1 4 4*4 4*6])'];

clear RE_Multipliers IHL_Multipliers;
RowNames = {'Lump-sum Y','Lump-sum C','Lump-sum I','Capital tax Y',...
    'Capital tax C','Capital tax I','Labor tax Y','Labor tax C','Labor tax I'};
Table4_PVMultipliers = array2table(Table4_PVMultipliers,'RowNames',RowNames, ...
    'VariableNames',{'RE Impact','RE 1 yr','RE 4 yrs','RE 6 yrs',...
    'AL Impact','AL 1 yr','AL 4 yrs','AL 6 yrs'});

%% Figure 5. Impulse responses for different fiscal policy specifications
% Impulse responses to an increase in government spending of 1% of GDP of
% the new Keynesian model for different fiscal policy specifications.

variablesMAT = {
    'c'     'Consumption ($C$)'
    'g'     'Government spending ($G$)'
    'i'     'Investment ($I$)'
    'k'     'Capital stock ($K$)'
    'mc'    'Real marginal cost ($MC$)'
    'n'     'Hours worked ($N$)'
    'pi'    'Inflation ($\Pi$)'
    'q'     'Tobin''s Q ($Q$)'
    'rk'    'Rental rate of capital ($r^k$)'
    'r'     'Nominal interest rate ($R$)'
    't'     'Lump-sum tax ($T$)'
    'w'     'Real wages ($W$)'
    'y'     'Output ($Y$)'
    'tauk'  'Capital income tax ($\tau^k$)'
    'tauw'  'Labour income tax ($\tau^w$)'
    };

figure('Name','Impulse responses for different fiscal policy specifications');
for iPlot = 1:size(variablesMAT,1)
    subplot(5,3,iPlot);
    hold on;
    plot(1:T,zeros(T,1),'-','Color',[.5 .5 .5]);
    % Rational expectations
    if strcmp(variablesMAT{iPlot,1},'t')        
        plot(1:T,NewKeynesianResults.REE_IRF_eg.y(:, ...
            strcmp(NewKeynesianResults.REE_IRF_eg.names,variablesMAT{iPlot,1})), ...
            'b')
    elseif strcmp(variablesMAT{iPlot,1},'tauk')
        plot(1:T,Distortionary1Results.tauk,'r');
    elseif strcmp(variablesMAT{iPlot,1},'tauw')
        plot(1:T,Distortionary2Results.tauw,'g');
    else
        a = plot(1:T,NewKeynesianResults.REE_IRF_eg.y(:, ...
            strcmp(NewKeynesianResults.REE_IRF_eg.names,variablesMAT{iPlot,1})), ...
            'b');
        b = plot(1:T,Distortionary1Results.REE_IRF_eg.y(:, ...
            strcmp(Distortionary1Results.REE_IRF_eg.names,variablesMAT{iPlot,1})), ...
            'r');
        c = plot(1:T,Distortionary2Results.REE_IRF_eg.y(:, ...
            strcmp(Distortionary2Results.REE_IRF_eg.names,variablesMAT{iPlot,1})), ...
            'g');
    end

    % Infinite Horizon Learning
    if strcmp(variablesMAT{iPlot,1},'t')
        Data = [NewKeynesianResults.IHL_eg.y; NewKeynesianResults.IHL_eg.w];
        VarNames = [NewKeynesianResults.IHL_eg.names_xy NewKeynesianResults.IHL_eg.names_z];
        plot(1:T,...
            Data(strcmp(VarNames,variablesMAT{iPlot,1}),:), ...
            '--b');
    elseif strcmp(variablesMAT{iPlot,1},'tauk')
        plot(1:T,Distortionary1Results.tauk,'r');
    elseif strcmp(variablesMAT{iPlot,1},'tauw')
        plot(1:T,Distortionary2Results.tauw,'g');
    else
        Data = [NewKeynesianResults.IHL_eg.y; NewKeynesianResults.IHL_eg.w];
        VarNames = [NewKeynesianResults.IHL_eg.names_xy NewKeynesianResults.IHL_eg.names_z];
        plot(1:T,...
            Data(strcmp(VarNames,variablesMAT{iPlot,1}),:), ...
            '--b');
        Data = [Distortionary1Results.IHL_eg.y; Distortionary1Results.IHL_eg.w];
        VarNames = [Distortionary1Results.IHL_eg.names_xy Distortionary1Results.IHL_eg.names_z];
        plot(1:T,...
            Data(strcmp(VarNames,variablesMAT{iPlot,1}),:), ...
            '--r');
        Data = [Distortionary2Results.IHL_eg.y; Distortionary2Results.IHL_eg.w];
        VarNames = [Distortionary2Results.IHL_eg.names_xy Distortionary2Results.IHL_eg.names_z];
        plot(1:T, ...
            Data(strcmp(VarNames,variablesMAT{iPlot,1}),:),'--g');
    end

    hold off;
    title(variablesMAT(iPlot,2),'interpreter','latex');    
end
legend([a b c],'Baseline','Capital tax financing','Labour tax financing',...
    'Orientation','horizontal','Location','best');