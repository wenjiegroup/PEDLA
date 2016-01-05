function [ x_new ] = format_x_locally( x ,neibour )
%FORMAT_X_LOCALLY Summary of this function goes here
%   Detailed explanation goes here

x_new=[];
[m,n]=size(x);
x=[zeros(neibour,n);x;zeros(neibour,n)];
for i=0:2*neibour
    x_new=[x_new x(1+i:m+i,:)];
end


end

