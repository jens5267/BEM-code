% Test script to fill in missing values from an xlsx file

clear all
close all
clc

%% Read data
% Read all data
data = txtread('NACA_2408_75k.txt','13:'); 

% Find the missing data
[row,col] = find(isnan(data)); 

% create loop to fill in each gap separately
for i = 1:length(row)
    % get row and col for each loop
    row_i = row(i);
    col_i = col(i);
    
    % find first value above this missing value in the same column
    % above
    k_up = find(~isnan(data(1:row_i,col_i)));
    index_up = row_i-k_up(end);
    
    % below
    index_down = find(~isnan(data(row_i+1:end, col_i)),1);
    
    % Fill in the weighted average in the data gap 
    data(row_i,col_i) = (data(k_up(end),col_i)*index_up...
        +data(row_i+index_down,col_i)*index_down)/(index_up+index_down)
    
end