function Output=PlotExp(Results,selected_vars)
%PlotExp Summary of this function goes here
%   Detailed explanation goes here

%selected_vars={'c','i','n','pi','q','rk','d','w','r'};

RL_names=Results.RL_eg.names_xy;
IHL_names=Results.IHL_eg.names_xy;

figure('Name','Expectations');

iz=1;
for ix=1:size(selected_vars,2)
    iy=find(strcmp(RL_names,selected_vars(ix)));
    iq=find(strcmp(IHL_names,selected_vars(ix)));
    yf_RE(ix,:)=Results.RL_eg.y(iy,2:end);
    yf_RL(ix,:)=Results.RL_eg.yf(iy,1:end-1);
    yf_IHL(ix,:)=Results.IHL_eg.yf(iy,2:end);
    error(ix,:)=yf_RE(ix,:)-yf_IHL(ix,:);
    subplot(3,3,iz);
    hold on;
    plot(yf_RE(ix,:));
    plot(yf_RL(ix,:),'--');
    plot(yf_IHL(ix,:),'--');
    hold off;
    title(selected_vars(ix));
    iz=iz+1;
    if iz==10        
        iz=1;
        if size(selected_vars,2)>9
            legend('RE','EE learning','IH learning');
            figure('Name','Expectations');            
        end
    end
end
legend('RE','EE learning','IH learning');

% errors
figure('Name','Forecast errors');
iz=1;
for ix=1:size(selected_vars,2)    
    subplot(3,3,iz);
    plot(error(ix,:));
    title(selected_vars(ix));
    iz=iz+1;
    if iz==10        
        iz=1;
        if size(selected_vars,2)>9
            legend('IHL error');
            figure('Name','Forecast errors');            
        end
    end
end

Output.names=selected_vars;
Output.yf_RL=yf_RL;
Output.yf_IHL=yf_IHL;
Output.yf_RE=yf_RE;
Output.errors=error;
end

