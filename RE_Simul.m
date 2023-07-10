%% Rational Expectations Equilibrium Simulation
% © 2014 Ewoud Quaghebeur

%%
function [Results] = RE_Simul(Omega,Phi,T,epsilon,names,draw)

%%
% Generates simulated data for the Rational Expectations Equilibrium (REE) 
%
% $$X_t = \Omega X_{t-1} + \Phi\epsilon_t$$
%
% * |Omega| and |Phi| are transition matrices
% * |T| is the lenght of the simulation period
% * |mu| and |sigma| are the mean and the st.dev. of the normal 
%     distribution for the error terms
% * |names| contains the labels of the variables 
% * |epsilon| is a k times T matrix with innovations in the shock processes.
% * |draw| determines whether the result should be plotted (1) or not 0.
%

k=size(Omega,1); % number of disturbances/shocks
y(:,1) = Omega*zeros(k,1)+Phi*epsilon(:,1); %endogenous variables in period 1
yf(:,1)= Omega*y(:,1); % expectations

for t=2:T % realisation of endogenous variables for periods 2...T
    y(:,t) = Omega*y(:,t-1)+Phi*epsilon(:,t);
    yf(:,t) = Omega*y(:,t);
end

if draw==1 % figures of simulated variables.
    i=1;
    figure('Name','REE simulation')
    for j=1:k           
        subplot(3,3,i); plot(1:T,y(j,:));
        title(names(j));
        i=i+1;
        if mod(j,9)==0
            figure('Name','REE simulation') 
            i=1;
        end         
    end
end
Results.y = y'; % function returns simulated data
Results.yf = yf'; % function returns expactations
Results.names = names; % function returns names of variables
end

