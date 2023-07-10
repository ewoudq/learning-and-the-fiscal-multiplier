function [result] = Textcoordinates(data,order_alpha)
%TEXTCOORDINATES transforms data matrix to strings with coordinates for
%Tikz file
%   INPUTS
%   - 'data' is a i*j matrix with i variables and j simulation periods
%   - 'order_alpha' is the 1*i cell array with names for the i variables
%

precision = 3; % number of significant digits
result = cell(size(data,1),2);

for i=1:size(data,1) % # of variables
    result{i,1} = order_alpha{i};
    result{i,2}='coordinates{';
    for j=1:size(data,2) % # of periods        
        result{i,2}=[char(result{i,2}) '(' int2str(j) ',' num2str(data(i,j), precision) ') '];
    end 
    result{i,2}=[char(result{i,2}), '};'];
end
end

