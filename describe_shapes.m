%% 
% Title: Describe shapes in terms of elliptic Fourier descriptors
% Author: Krishna Manaswi Digumarti
% Version: 3.0
% Date: Sep 2020
% Description: This code take a set of images and for each image and
% computes the elliptic Fourier descriptors of shape.
%
% I would appreciate it if you cite the following paper for which this code
% was originally developed 
% Digumarti KM, Trimmer B, Conn AT, Rossiter J. 
% "Quantifying Dynamic Shapes in Soft Morphologies."
% Soft Robotics. 6(6), pp.733-744. 2019
%
% Based on: 
% F. P. Kuhl and C. R. Giardina, “Elliptic fourier features of a closed contour,”
% Computer graphics and image processing, vol. 18, no. 3, pp. 236–258, 1982.
% 
% Freeman chain encoding is based on:
% H. Freeman, “Computer processing of line-drawing images,” 
% ACM Computing Surveys (CSUR), vol. 6, no. 1, pp. 57–97, 1974.
%
% Parts of the code have been adapted from:
% Auralius Manurung, "Elliptic fourier for shape analysis"
% https://uk.mathworks.com/matlabcentral/fileexchange/32800-elliptic-fourier-for-shape-analysis

%% tabula rasa
clear all
close all
clc

%% some settable properties
nHarmonics = 6; % number of harmonics to use for estimation
nSynthesis = 100; % num of pts to reconstruct contour from estimate

shouldNormalize = 0; % 1 = yes / 0 = no: use 0 for visualization, 1 for pca
shouldVisualize = 1; % 1 = yes / 0 = no: to visualised estimated contour
% Note: only normalised coefficients are used in pca. 
% To speed up processing, avoid visualisation when computing coefficients 
% for further analysis. Set shouldNormalize to 1 and shouldVisualize to 0

color = [1 0 0];
lineWidth = 2;

% specify folder to read image frames from
folderForFrames = 'frames';
recordingName = 'movie';
subFolder = strcat(folderForFrames, '/', recordingName);

% figure out how many frames there are in total
a = dir([subFolder, '/*.png']);
numberOfFrames = numel(a);

%% compute chain code and harmonic coeffieints for each frame
coeffs_Mat = []; % a variable that stores the Fourier coefficients

for frameNum = 1:1:numberOfFrames
    % read image and binarize it
    img = imread(strcat(subFolder, '/frame', num2str(frameNum), '.png'));
    img_gray = rgb2gray(img);           % convert image to grayscale
    img_binary = imbinarize(img_gray);  % make it a binary image
    img_bw = double(img_binary);        % convert from logical to double
    img_bw(img_bw == 1) = 255;          % make black = 0, white = 255
    
    if shouldVisualize
        figure(1)
        imshow(img_bw)
        hold on
    end
    
    % get left most point on boundary - look for first non zero element
    [startRow,startCol] = find(img_bw,1,'first'); 
    % cols go left to right, rows go top to bottom.

    % determine the direction in which to start moving along boundary
    if img_bw(startRow+1,startCol+1) == 255
        firstStep = 'SE';
    elseif img_bw(startRow,startCol+1) == 255
        firstStep = 'E';
    elseif img_bw(startRow-1,startCol+1) == 255
        firstStep = 'NE';
    elseif img_bw(startRow-1,startCol) == 255
        firstStep = 'N';
    elseif img_bw(startRow-1,startCol-1) == 255
        firstStep = 'NW';
    elseif img_bw(startRow,startCol-1) == 255
        firstStep = 'W';
    elseif img_bw(startRow+1,startCol-1) == 255
        firstStep = 'SW';
    elseif img_bw(startRow+1,startCol) == 255
        firstStep = 'S';
    else
        assert(0);
    end

    % trae the boundary using the starting position and direction
    chain = bwtraceboundary(img_bw,[startRow,startCol],firstStep);
    chain = rot90(chain,2);
    chainCode=[];

    % assign the Freeman chain code
    for i=1:1:size(chain,1)-1
        dx = chain(i+1,1)-chain(i,1);
        dy = chain(i+1,2)-chain(i,2);

        if dx == 1 && dy == 0
            chainCode(i,1) = 0;
        elseif dx == 1 && dy == 1
            chainCode(i,1) = 1;
        elseif dx == 0 && dy == 1
            chainCode(i,1) = 2;
        elseif dx == -1 && dy == 1
            chainCode(i,1) = 3;
        elseif dx == -1 && dy == 0
            chainCode(i,1) = 4;
        elseif dx == -1 && dy == -1
            chainCode(i,1) = 5;
        elseif dx == 0 && dy == -1
            chainCode(i,1) = 6;
        elseif dx == 1 && dy == -1
            chainCode(i,1) = 7;
        end 
    end
    
    % get harmonic coefficients
    coefficients = fourier_approx(transpose(chainCode), nHarmonics, shouldNormalize);

    A0 = coefficients(1,1);
    C0 = coefficients(1,3);
    a = coefficients(2:end,1);
    b = coefficients(2:end,2);
    c = coefficients(2:end,3);
    d = coefficients(2:end,4);
    
    % collect coefficients from each frame in a matrix
    coeffs = reshape(coefficients', [1, size(coefficients,1)*size(coefficients,2)]);
    coeffs_Mat = [coeffs_Mat; coeffs];
    
    % Optional - visualize the estimated contour
    if shouldVisualize
        % synthesize the estimated contour
        coordinates = zeros(nSynthesis,2);
        for j = 1 : nSynthesis
            x_ = 0.0;
            y_ = 0.0;

            for i = 1 : nHarmonics
                x_ = x_ + (a(i) * cos(2 * i * pi * j / nSynthesis) + b(i) * sin(2 * i * pi * j / nSynthesis));
                y_ = y_ + (c(i) * cos(2 * i * pi * j / nSynthesis) + d(i) * sin(2 * i * pi * j / nSynthesis));
            end

            coordinates(j,1) = A0 + x_;
            coordinates(j,2) = C0 + y_;
        end

        % correct location
        coordinates = coordinates + [startCol*ones(size(x_,1),1) startRow*ones(size(x_,1),1)];

        % Make it closed contour
        contour = [coordinates; coordinates(1,1) coordinates(1,2)];

        % draw the synthesized contour
        figure(1)
        plot(contour(:,1), contour(:,2), 'color', color, 'linewidth', lineWidth);
        drawnow
    end
end

% save coefficients - Note: only normalized coefficients are used in pca
if shouldNormalize
    save('coefficients.mat','coeffs_Mat')
end