function setpath()
%   MAKE adds paths of the DeComp to Matlab.

%  Copyright (c) 2017 Haizhao Yang, National University of Singapore
%  This file is distributed under the terms of the MIT License.

global SetPath
global CSPT
global MSPT

type = computer;

if strcmp(type,'MAC2'),
    CSPT = ':';
    SetPath = [pwd, CSPT];
    MSPT = ';';
elseif isunix,
    % Mac OS X returns isunix=1
    CSPT = '/';
    SetPath = [pwd, CSPT];
    MSPT = ':';
elseif strcmp(type(1:2),'PC');
    CSPT = '\';
    SetPath = [pwd, CSPT];
    MSPT = ';';
end

disp('Begin to set MATLAB path...')

file_path = mfilename('fullpath');
tmp = strfind(file_path,'setpath');
file_path = file_path(1:(tmp(end)-1));

% Foulder for all soource files recursively
addpath(genpath([file_path 'Applications']));
addpath(genpath([file_path 'results']));
addpath(genpath([file_path 'Source']));
addpath(genpath([file_path 'External']));

disp('Begin to compile MEX files...');
rootDir = pwd;
cd(['Source' CSPT 'SynLab' CSPT 'src' CSPT 'SS_CT_2D' CSPT 'src' CSPT]);
mex SS_polar.c;
mex SS_polar_old.c;
mex SS_polar_v2.c;
mex SS_polar_v1.c;
cd(rootDir);
cd(['Applications' CSPT 'MaterialsSci' CSPT 'SynCrystal' CSPT 'SSTmethod' CSPT 'src' CSPT]);
mex LocWeight.c;
mex LocSmooth.c;
cd(rootDir);
cd(['Applications' CSPT 'MaterialsSci' CSPT 'SynCrystal' CSPT 'VarSSTmethod' CSPT 'src' CSPT 'srcSST' CSPT]);
mex LocWeight.c;
mex LocWeight_v2.c;
cd(rootDir);
cd(['Applications' CSPT 'Art&Hist' CSPT 'CanvasThreadCount' CSPT 'src' CSPT]);
mex LocWavVec_v2.c;
cd(rootDir);


disp('Path set!');

clear tempPath front back
clear SetPath MATLABVERSION CSPT
clear type MSPT
end