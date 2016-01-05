function [ M2 ] = softmaxtup( theta , data )
%SOFTMAXTUP Summary of this function goes here
%   Detailed explanation goes here

M1 = theta * data;
M1 = bsxfun(@minus, M1, max(M1, [], 1)); %%Preventing overflows 
M2 = exp(M1);
M2 = bsxfun(@rdivide, M2, sum(M2));

end

