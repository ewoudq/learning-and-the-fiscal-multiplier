%% 
% _Order and reorder vectors_
%
% This cell reorders the variables alphabetically. The reordered variables 
% of |VARNAMES| are stored in the cell array |names|.
%
% * |order_alpha|: variable names  of |VARNAMES| in alphabetical order
% * |reorder_alpha|: vector to reorder variables in |VARNAMES| from x-y-z 
% to alphabetical order.

names = cellstr(VARNAMES);
[order_alpha,reorder_alpha] = sort(names);
names_x=sort(names_x);
names_y=sort(names_y);
names_z=sort(names_z);
names_disturbances=sort(names_disturbances);
names_x_obs=sort(names_x_obs);
names_z_obs=sort(names_z_obs);
irf_shock=sort(irf_shock);

%%
% _Location of the observed variables and shocks of interest_
%
% Determines the location of the observed states, the observed shocks, and
% the distrurbance(s) for the impulse response analysis in |names_x|, 
% |names_z|, % and |names_disturbances|, respectively. 

loc_x_obs=zeros(1,size(names_x_obs,2)); %location of the observed states
for i=1:size(names_x_obs,2)
    for j=1:size(names_x,2)
        if strcmp(names_x_obs(i),names_x(j))            
            loc_x_obs(i)=j;
        end
    end
end
loc_z_obs=zeros(1,size(names_z_obs,2)); % location of the observed shocks
for i=1:size(names_z_obs,2)
    for j=1:size(names_z,2)
        if strcmp(names_z_obs(i),names_z(j))            
            loc_z_obs(i)=j;
        end
    end
end
loc_shock=zeros(1,size(irf_shock,2)); % location of the disturbance(s) for the impulse response analysis
for j=1:size(irf_shock,2)
    for i=1:size(names_disturbances,2)
        if strcmp(irf_shock(j),names_disturbances(i))            
            loc_shock(j)=i;
        end
    end
end