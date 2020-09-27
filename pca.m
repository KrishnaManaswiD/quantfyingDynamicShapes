%% 
% Title: Dimensionality reduction and principal component analysis
% Author: Krishna Manaswi Digumarti
% Version: 3.0
% Date: Sep 2020
% Description: First, the code takes the set of harmonic coefficients and
% reduces them down to a chosen (3 here) number of principal components.
% Then it provides an interactive way to see what the shape looks like when
% scores on the principal components are changed.
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

%% plot data points along each principal componet
figure(2)
for k = 1:size(newData,2)
    subplot(size(newData,2),1,k)
    plot(1:s1, newData(1:s1,k), 'o',...
        'MarkerFaceColor',[0.2316735 ,  0.3181058 ,  0.54483444],...
        'MarkerEdgeColor', [0.2316735 ,  0.3181058 ,  0.54483444]);
    %set(gca,'Xticklabel',[]) % rid of the numbers but leave the ticks.
    hold on
    xlabel('Frame number')
    ylabel(['PC ', num2str(k)])
    if k == 1
        title('Coefficients data represented in terms of principal components')
    end
end

%% draw the mean shape
% Note: This code is written for 3 principal components. If you are
% analysing more components, add weights, sliders and update functions.

mn = mean(data);
ext = mean(coeffs_Mat(:,5:7));

f = figure(3);
ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);
myHandles = guihandles(f);
weights = zeros(1,3);
myHandles.w1 = weights(1);
myHandles.w2 = weights(2);
myHandles.w3 = weights(3);

myHandles.mn = mn;
myHandles.mean_ = mean_;
myHandles.stdDev_ = stdDev_;
myHandles.ext = ext;
myHandles.eV = eV;
myHandles.nHarmonics = nHarmonics;
myHandles.nSynthesis = nSynthesis;
guidata(f,myHandles);

% create sliders
slider1 = uicontrol('Parent',f,'Style','slider','Position',[81,110,419,20],...
              'value',weights(1), 'min',-2, 'max',2, 'CreateFcn', @updatePlot, 'Callback', @updateW1);
slider2 = uicontrol('Parent',f,'Style','slider','Position',[81,85,419,20],...
              'value',weights(2), 'min',-2, 'max',2, 'CreateFcn', @updatePlot, 'Callback', @updateW2);
slider3 = uicontrol('Parent',f,'Style','slider','Position',[81,60,419,20],...
              'value',weights(3), 'min',-2, 'max',2, 'CreateFcn', @updatePlot, 'Callback', @updateW3);          

function updateW1(source,~)
    myHandles = guidata(gcbo);
    myHandles.w1 = source.Value;
    guidata(gcbo,myHandles)
    updatePlot()
end

function updateW2(source,~)
    myHandles = guidata(gcbo);
    myHandles.w2 = source.Value;
    guidata(gcbo,myHandles);
    updatePlot()
end

function updateW3(source,~)
    myHandles = guidata(gcbo);
    myHandles.w3 = source.Value;
    guidata(gcbo,myHandles)
    updatePlot()
end

function updatePlot(~,~)
    myHandles = guidata(gcbo);
    mn = myHandles.mn;
    ext = myHandles.ext;
    mean_ = myHandles.mean_;
    stdDev_ = myHandles.stdDev_;
    nHarmonics_ = myHandles.nHarmonics;
    nSynthesis_ = myHandles.nSynthesis;
    eV = myHandles.eV;
    weights = [myHandles.w1,myHandles.w2,myHandles.w3];
    add2 = [0;0;0;eV*(mean_+weights.*stdDev_)']';
    p1 = [ext, mn] + add2; 
%     disp(mean_);
%     disp(weights);
%     disp(weights.*stdDev_);
    p = reshape(p1,[4,nHarmonics_])';

    a = p(:,1);
    b = p(:,2);
    c = p(:,3);
    d = p(:,4);
    
    coordinates = zeros(nSynthesis_,2);
    for t = 1 : nSynthesis_                                                                                                                                                        
            x_ = 0.0;
            y_ = 0.0;

            for i = 1 : nHarmonics_
                x_ = x_ + (a(i) * cos(2 * i * pi * t / nSynthesis_) + b(i) * sin(2 * i * pi * t / nSynthesis_));
                y_ = y_ + (c(i) * cos(2 * i * pi * t / nSynthesis_) + d(i) * sin(2 * i * pi * t / nSynthesis_));
            end

            coordinates(t,1) = x_;
            coordinates(t,2) = y_;
    end

    contour = [coordinates; coordinates(1,:)]; % close the contour
    plot(contour(:,1), contour(:,2),'b','LineWidth',3)
    axis equal
end