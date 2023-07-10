%% Map from the Perceived Law of Motion (PLM) to the Actual Law of Motion (ALM)
% Ewoud Quaghebeur, Ghent University, April 2014.

%%
function [a_ALM,b_ALM,c_ALM] = T_map(a_PLM,b_PLM,c_PLM,alpha_mat,beta_mat,delta_mat,kappa_mat,phi_mat)

%%
% This function returns the coefficients of the Actual Law of Motion (ALM)
%
% $$y_{t}=a+by_{t-1}+cw_t$$
%
% of the macroeconomic model 
%
% $$y_t=\alpha+\beta E^*_ty_{t+1}+\delta y_{t-1}+\kappa w_t$$
%
% $$w_t=\varphi w_{t-1}+e_t$$
%
% when the Perceived Law of Motion is 
%
% $$E_ty_{t+1}=a+by_t+c\varphi w_t$$
%
% The function inputs are
%
% * The coefficients "a_PLM", "b_PLM", and "c_PLM" of the Perceived Law of Motion
% * The parameters of the model "alpha_mat", "beta_mat" , "delta_mat",
% "kappa_mat", and "phi_mat"
%
% The function returns the coefficients of the Actual law of motion:
% "a_ALM", "b_ALM", and "c_ALM".
% See Evans and Honkapohja, 2001, p. 237-238.
%
% #EQ 2015061This is the 'unprojected' ALM. See Evans and Honkapohja, 2001, p. 321-322
% and Guse, 2007, JEDC, p. 1522

pre_mat=eye(size(alpha_mat,1))-beta_mat*b_PLM;
a_ALM=pre_mat\(alpha_mat+beta_mat*a_PLM);
b_ALM=pre_mat\delta_mat;
c_ALM=pre_mat\(kappa_mat+beta_mat*c_PLM*phi_mat);
end

