function [X,Y] = loaddata(ind)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
if ind == 1
    load('glass_uni.mat');
elseif ind ==2
    load('yeast_uni.mat');
elseif ind ==3
    load('wine_uni.mat');
elseif ind ==4
    load('german_uni.mat');
elseif ind ==5
    load('dermatology_uni.mat');
elseif ind ==6
    load('JAFFE_1024.mat');
elseif ind ==7
    load('coil-20_1000.mat');
elseif ind ==8
    load('MSRA25_1024.mat');
elseif ind ==9
    load('minist_5000_50.mat');
elseif ind ==10
    load('PIE_vec.mat');
elseif ind ==11
    load('X_highdimen.mat');    
end
X = double(X);
Y = Y;
end


