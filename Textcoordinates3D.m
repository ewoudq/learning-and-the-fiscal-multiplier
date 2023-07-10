function [result] = Textcoordinates3D(data,order_alpha)
%TEXTCOORDINATES transforms data matrix to strings with coordinates for
%Tikz file
%   INPUTS
%   - 'data' is a k*j matrix with k parameter values and j simulation periods
%   - 'order_alpha' is the 1*i cell array with names for the i variables
%

precision = 3; % number of significant digits
result = cell(size(order_alpha,2),2);
eval(['nrpar = size(data.', order_alpha{1}, ',1)']) % # of parameters
eval(['nrper = size(data.', order_alpha{1}, ',2)'])% # of periods

for i=1:size(order_alpha,2) % # of variables
    result{i,1} = order_alpha{i};
    result{i,2}='coordinates{';
    for k=1:nrpar % # of parameters (y)
        for j=1:nrper % # of periods (x)
            eval(['value = data.',order_alpha{i},'(k,j)'])
            result{i,2}=[char(result{i,2}) '(' int2str(j) ',' int2str(k) ',' num2str(value, precision) ') '];
        end
    end 
    result{i,2}=[char(result{i,2}), '};'];
end
end

