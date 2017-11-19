% This code analyze real canvas images in our paper. The data "RyderPart" is
% available upon request. The other data is from the RKD database available
% at http://english.rkd.nl/Services/image-sharing
%
% Images have to be stored in a mat file with a variable name as X. Please
% see the code fastAnalysis.m for more information about the input.
%
% By Haizhao Yang

% Clean data
%fastAnalysis('F205',1,0.8,1e-4,1,[8,5],15,90,2,4,4);

% Noisy data
%fastAnalysis('F659',0.8,0.625,1e-4,1,[8,5],15,90,2,4,4);
%fastAnalysis('L11',0.8,0.625,1e-4,1,[8,5],15,90,2,4,4);
%fastAnalysis('L17',0.75,0.625,1e-4,1,[8,15],15,90,2,4,4);
%fastAnalysis('L30',0.8,0.625,1e-4,1,[8,5],15,90,2,4,4);

% Twill
%fastAnalysis('RyderPart',0.8,0.625,1e-7,1,[10,20],6,180,1,8,2);


