function main()
clc;clear all;close all
load('EEGData')
Kmax = 15; % The length of actions
%% a six-dimensional phase space
subplot(121)
Y = x(:,1:5); 
F = YaghoobiFractalDimension(Y,Kmax);
disp(['Fractal Dimension: ',num2str(F)])
%% a three-dimensional phase space
subplot(122)
Y = x(:,7:9); 
F = YaghoobiFractalDimension(Y,Kmax);
disp(['Fractal Dimension: ',num2str(F)])
set(gcf,'pos',[100 100 900 400])
end

function YFD = YaghoobiFractalDimension(Y,Kmax)
%% Estimate Yaghoobi fractal dimension
% Y: a high dimension trajectory
% Kmax: The length of actions
% YFD: Fractal Dimension of the input signal

%% Copyright RYK, 2021
% Reza Yaghoobi Karimui
% Freelance Researchers Association
% E-mail: reza_yaghoby@yahoo.ca
% Website: https://apexpg.jimdofree.com
% Reference: A new approach to measure the fractal dimension of a trajectory in the high-dimensional phase space
%
dim = size(Y,2);
N = size(Y,1);
for k = 1:Kmax
    clear xx L
    for i = 1:k
        a = Y(i:k:(i+fix((N-i)/k)*k),:);
        d = diff(a);
        r = sqrt(sum(d.^2,2));
        alpha = [N-1]/(k*floor([N-i]/k));
        L(i) = alpha*sum(r)/(k^(dim-1));
    end
    LM(k) = mean(L);
end
LM(LM==0)=1e-200;
lLM = log10(LM);
lk = log10(1:Kmax);
lk(LM==0) = [];
lLM(LM==0) =[];
%% LSE
[a,b] = LSE(lk,lLM);
%% Plot
plot(lk,a*lk+b,':k');hold on
plot(lk,lLM,'.','MarkerSize',15);
xlabel('Log_1_0(k)')
ylabel('Log_1_0(L_k)')
xlim([lk(1) lk(end)])
YFD = abs(a);
title(['FD = ',num2str(YFD)])
end
%% fitting
function [a,b] = LSE(x,y)
N = size(x,2);
a = [N*x*y'-sum(x)*sum(y)]/[N*sum(x.^2)-sum(x)^2];
b = [sum(x.^2)*sum(y)-x*y'*sum(x)]/[N*sum(x.^2)-sum(x)^2];
end

