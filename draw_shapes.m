%% 
% Title: Draw estimated shapes by varying scores on principal components
% Author: Krishna Manaswi Digumarti
% Version: 3.0
% Date: Sep 2020
% Description: This is a less interactive version of the pcs analysis code.
% We start with the shape constructed out of the mean values of scores
% on principal components. Then, the score on each component is varyied by
% one and two standard deviation on each side, while keeping the score on
% other components at the mean. The resultant shapes are drawn next to each
% other. 
%
% Each row: mean - 2*std, mean - std, mean, mean + std, mean + 2*std.
% Other coponents fixed at their means
%
% I would appreciate it if you cite the following paper for which this code
% was originally developed 
% Digumarti KM, Trimmer B, Conn AT, Rossiter J. 
% "Quantifying Dynamic Shapes in Soft Morphologies."
% Soft Robotics. 6(6), pp.733-744. 2019

%% Tabula rasa
clear all
close all
clc

%% load coefficients
load('coefficients.mat')
data = coeffs_Mat(:,8:end);
s1 = size(data,1);

nHarmonics = size(coeffs_Mat,2)/4 - 1; % deduce num harmonics frm saved data
nSynthesis = 100;

%% reduce data dimensions
% make data zero mean
data_zeroMean = data - mean(data);

% reduce data to a specified number of principal components
numEigs = 3; % number of eigenvectors to reduce the coefficeints to
covMat=cov(data_zeroMean);
covMat_norm = 2*(covMat - min(min(covMat))) ./ (max(max(covMat)) - min(min(covMat)))-1;
[eV, eD] = eigs(covMat,numEigs);

% represent the data in terms of the principal components
newData = data_zeroMean*eV;
new_covMat = cov(newData);
new_covMat2 = eV' * covMat * eV;
mean_ = mean(newData);

stdDev_ = diag(new_covMat2)';
stdDev_ = sqrt(stdDev_);

%% draw the mean shape
% Note: this code is written for 3 principal components. Add more weights
% and subplots if you want to look at more components. The code is
% deliberately written in a repeating manner to make it easier to read.

mn = mean(data);
ext = mean(coeffs_Mat(:,5:7));

f = figure(4);

% varying PC 1
subplot(3,5,1)
w1 = -2; w2 = 0; w3 = 0;
drawShape(mn,mean_,stdDev_,ext,eV,w1,w2,w3,nHarmonics,nSynthesis);

subplot(3,5,2)
w1 = -1; w2 = 0; w3 = 0;
drawShape(mn,mean_,stdDev_,ext,eV,w1,w2,w3,nHarmonics,nSynthesis);

subplot(3,5,3)
w1 = 0; w2 = 0; w3 = 0; 
drawShape(mn,mean_,stdDev_,ext,eV,w1,w2,w3,nHarmonics,nSynthesis);

subplot(3,5,4)
w1 = 1; w2 = 0; w3 = 0; 
drawShape(mn,mean_,stdDev_,ext,eV,w1,w2,w3,nHarmonics,nSynthesis);

subplot(3,5,5)
w1 = 2; w2 = 0; w3 = 0;
drawShape(mn,mean_,stdDev_,ext,eV,w1,w2,w3,nHarmonics,nSynthesis);

% varying PC 2
subplot(3,5,6)
w2 = -2; w1 = 0; w3 = 0; 
drawShape(mn,mean_,stdDev_,ext,eV,w1,w2,w3,nHarmonics,nSynthesis);

subplot(3,5,7)
w2 = -1; w1 = 0; w3 = 0; 
drawShape(mn,mean_,stdDev_,ext,eV,w1,w2,w3,nHarmonics,nSynthesis);

subplot(3,5,8)
w2 = 0; w1 = 0; w3 = 0; 
drawShape(mn,mean_,stdDev_,ext,eV,w1,w2,w3,nHarmonics,nSynthesis);

subplot(3,5,9)
w2 = 1; w1 = 0; w3 = 0; 
drawShape(mn,mean_,stdDev_,ext,eV,w1,w2,w3,nHarmonics,nSynthesis);

subplot(3,5,10)
w2 = 2; w1 = 0; w3 = 0; 
drawShape(mn,mean_,stdDev_,ext,eV,w1,w2,w3,nHarmonics,nSynthesis);

% varying PC 3
subplot(3,5,11)
w3 = -2; w1 = 0; w2 = 0;
drawShape(mn,mean_,stdDev_,ext,eV,w1,w2,w3,nHarmonics,nSynthesis);

subplot(3,5,12)
w3 = -1; w1 = 0; w2 = 0; 
drawShape(mn,mean_,stdDev_,ext,eV,w1,w2,w3,nHarmonics,nSynthesis);

subplot(3,5,13)
w3 = 0; w1 = 0; w2 = 0; 
drawShape(mn,mean_,stdDev_,ext,eV,w1,w2,w3,nHarmonics,nSynthesis);

subplot(3,5,14)
w3 = 1; w1 = 0; w2 = 0; 
drawShape(mn,mean_,stdDev_,ext,eV,w1,w2,w3,nHarmonics,nSynthesis);

subplot(3,5,15)
w3 = 2; w1 = 0; w2 = 0; w4 = 0; w5 = 0;
drawShape(mn,mean_,stdDev_,ext,eV,w1,w2,w3,nHarmonics,nSynthesis);

function drawShape(mn,mean_,stdDev_,ext,eV,w1,w2,w3,nHarmonics,nSynthesis)
    weights = [w1,w2,w3];
    add2 = [0;0;0;eV*(mean_+weights.*stdDev_)']';
    p1 = [ext, mn] + add2; 
%     disp(mean_);
%     disp(weights);
%     disp(weights.*stdDev_);
    p = reshape(p1,[4,nHarmonics])';

    a = p(:,1);
    b = p(:,2);
    c = p(:,3);
    d = p(:,4);

    crd = zeros(nSynthesis+1,2);
    for t = 1 : nSynthesis                                                                                                                                                        
            x_ = 0.0;
            y_ = 0.0;

            for i = 1 : nHarmonics
                x_ = x_ + (a(i) * cos(2 * i * pi * t / nSynthesis) + b(i) * sin(2 * i * pi * t / nSynthesis));
                y_ = y_ + (c(i) * cos(2 * i * pi * t / nSynthesis) + d(i) * sin(2 * i * pi * t / nSynthesis));
            end

            crd(t,1) = x_;
            crd(t,2) = y_;

    end
    crd(end,:) = crd(1,:); % close the contour
    
    % draw the shape
    plot(crd(:,1), crd(:,2),'color',[0.2316735 ,  0.3181058 ,  0.54483444] ,'LineWidth',3)
    set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[]);
    set(gca,'Visible','off')
    axis equal
end