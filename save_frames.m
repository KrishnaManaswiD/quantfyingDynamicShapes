%% 
% Title: Extract frames from a video and save them as png
% Author: Krishna Manaswi Digumarti
% Version: 3.0
% Date: Sep 2020
%
% I would appreciate it if you cite the following paper for which this code
% was originally developed 
% Digumarti KM, Trimmer B, Conn AT, Rossiter J. 
% "Quantifying Dynamic Shapes in Soft Morphologies."
% Soft Robotics. 6(6), pp.733-744. 2019

%% tabula rasa
clear all
close all
clc

%% construct subfolder to save frames from movie
folderForFrames = 'frames';
recordingName = 'movie';           

% create sub folder
subFolder = strcat(folderForFrames, '/', recordingName);

if ~exist(subFolder, 'dir')
   mkdir(subFolder); 
end

%% save frames from a video
v = VideoReader(strcat(recordingName, '.mp4'));
frameNum = 0;

while hasFrame(v)  
    img_frame = readFrame(v);
    frameNum=frameNum+1;
    imwrite(img_frame,[subFolder, '\', 'frame', int2str(frameNum), '.png']);
end
