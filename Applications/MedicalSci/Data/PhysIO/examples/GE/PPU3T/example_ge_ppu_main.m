% Script that executes pulse oximetry (PPU) 3T GE logfiles. Just press play (F5)
%
% 
% Note: 
% - This is the input script to the PhysIO toolbox. Only this file has to be adapted for your study.
% - For documentation of any of the defined substructures here, please
% see also tapas_physio_new.m or the Manual_PhysIO-file.
%
% Copyright (C) 2013, Institute for Biomedical Engineering, ETH/Uni Zurich.
%
% This file is part of the PhysIO toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$
%

%% 0. Put code directory into path; for some options, SPM should also be in the path
pathRETROICORcode = fullfile(fileparts(mfilename('fullpath')), ...
    '../../../../public/code');

addpath(genpath(pathRETROICORcode));

physio      = tapas_physio_new();
log_files   = physio.log_files;
preproc     = physio.preproc;
sqpar       = physio.scan_timing.sqpar;
sync        = physio.scan_timing.sync;
model       = physio.model;
verbose     = physio.verbose;

%% 1. Define Input Files

log_files.vendor            = 'GE';

% Simple case
log_files.cardiac           = 'ECGData_epiRT_1025201115_44_19';      
log_files.respiration       = 'RespData_epiRT_1025201115_44_19'; 

% Difficult case
% log_files.cardiac           = 'ECGData_epiRT_phys_0921201215_38_08';      
% log_files.respiration       = 'RespData_epiRT_phys_0921201215_38_08'; 


log_files.sampling_interval = 25e-3;


%% 2. Define Nominal Sequence Parameter (Scan Timing)

% 2.1. Counting scans and dummy volumes from end of run, i.e. logfile
sqpar.Nslices           = 35;
sqpar.NslicesPerBeat    = 35;
sqpar.TR                = 1.925;
sqpar.Ndummies          = 5;
sqpar.Nscans            = 434;
sqpar.onset_slice       = [1 20 30];

sqpar.time_slice_to_slice = sqpar.TR/sqpar.Nslices; % equidistant slice spacing
sqpar.Nprep = 0; % start counting from beginning, not end of file


%% 3. Order of RETROICOR-expansions for cardiac, respiratory and
%% interaction terms. Option to orthogonalise regressors

model.type = 'none'; % 'RETROICOR';
% model.type = 'RETROICOR+HRV+RVT'; % 'RETROICOR';
model.order = struct('c',3,'r',4,'cr',1, 'orthogonalise', 'none');
model.input_other_multiple_regressors = ''; % either txt-file or mat-file with variable R
model.output_multiple_regressors = 'multiple_regressors.txt';

%% 4. Define Gradient Thresholds to Infer Gradient Timing (Philips only)
sync.method = 'nominal';

%% 5. Define which Cardiac Data Shall be Used

%% 4.1. Using plethysmograph curve with peak preprocolding
preproc.cardiac.modality = 'OXY'; % 'ECG' or 'OXY' (for pulse oximetry)
preproc.cardiac.initial_cpulse_select.min = 0.4;
preproc.cardiac.initial_cpulse_select.method = 'auto_matched';


%% 6. Output Figures to be generated

verbose.level = 2;
% 0 = none; 
% 1 = main plots (default); 
% 2 = debugging plots: for missed slice/volume events, missed heartbeats, 1D time series of created regressors
% 3 = all plots, incl. cardiac/respiratory phase estimation,
%     slice-to-volume assignment
verbose.fig_output_file = 'PhysIO_output.ps';

%% 7. Run the main script with defined parameters

physio.log_files            = log_files;
physio.preproc              = preproc;
physio.scan_timing.sqpar    = sqpar;
physio.scan_timing.sync     = sync;
physio.model                = model;
physio.verbose              = verbose;

[physio_out, R, ons_secs] = tapas_physio_main_create_regressors(physio);
