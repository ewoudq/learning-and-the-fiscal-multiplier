function Plots(y,w,beta0,a_new,b_new,c_new,names_x,names_y,names_z,names_x_obs,names_z_obs,constant)
%Plots Summary of this function goes here
%   Detailed explanation goes here

names_xy = [names_x names_y];
T=size(y,2);
nm=size(names_xy,2);
m=size(names_x,2);
k=size(names_z,2);
m_obs=size(names_x_obs,2); %new
k_obs=size(names_z_obs,2); %new

b0=beta0(1+constant:m_obs+constant,:)';
c0=beta0(m_obs+constant+1:end,:)';


i=1;% Endogenous variables
figure('Name','Restricted Learning simulation: endogenous variables');
for j=1:nm
    subplot(3,3,i);plot(1:T,y(j,:));
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
    subplot(3,3,i);plot(1:T,w(j,:));
    title(names_z(j));
    i=i+1;
    if mod(j,9)==0
        figure('Name','Restricted Learning simulation: stochastic processes');
        i=1;
    end
end
% Beliefs
figure('Name','Restricted Learning simulation: Beliefs a');%   beliefs a
plot(1:T,a_new,'-',1:T,a_new(:,1)*ones(1,T),'--');
title('Beliefs a')
% beliefs b
for xx=1:m_obs
    titl = sprintf('Restricted Learning simulation: Beliefs b: state %s.', names_x_obs{xx});
    figure('Name',titl);
    i=1;
    series=b_new((xx-1)*(nm)+1:(nm)*xx,:);
    for j=1:size(series,1)
        str = sprintf('Beliefs b for %s.', names_xy{j});
        %subplot(3,3,i);plot(1:T,series(j,:),'-',1:T,series(j,1)*ones(1,T),'--');
        subplot(3,3,i);plot(1:T,series(j,:),'-',1:T,b0(j,xx)*ones(1,T),'--');
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
    series=c_new((zz-1)*(nm)+1:(nm)*zz,:);
    for j=1:size(series,1)
        str = sprintf('Beliefs c for %s.', names_xy{j});
        %subplot(3,3,i);plot(1:T,series(j,:),'-',1:T,series(j,1)*ones(1,T),'--');
        subplot(3,3,i);plot(1:T,series(j,:),'-',1:T,c0(j,zz)*ones(1,T),'--');
        title(str)
        i=i+1;
        if mod(j,9)==0
            figure('Name',titl);
            i=1;
        end
    end
end
end