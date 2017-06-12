clear all;
% Forecasting private sector jobs using an autoregressive forecasting model
% Author: Jesse Weisberg
% 
% Data from: http://research.stlouisfed.org/fred2/series/USPRIV
%

[X,Y1,Z1] = xlsread('Jobs(2)');
M = 6; %number of months ahead for forecast
N = length(X);
n = 12; %number of months of data used in predicting
m = N-M-n+1;
H = zeros(m,n);
Z = zeros(m,1);
P = zeros(m,2);
for I = 1:m
    for J = 1:n
        %H(I,J)= X(n+I-J)/X(I+M+n-1); %scale actual results to 1 (error
        %minimized is the sume of sqaures of percent errors
        H(I,J)= X(n+I-J);
    end
    %Z(I) = X(I+M+n-1)/X(I+M+n-1); %scale actual results to 1
    Z(I) = X(I+M+n-1)
end
W = (H'*H)^-1*H'*Z; %weights calculated for for model
Y = H*W;
SumSquares = (Z-Y)'*(Z-Y);
Error = sqrt(SumSquares/m);
for I = 1:m
 P(I,1) = Y(I); %estimated
 P(I,2) = Z(I); %measured
end
plot(P)