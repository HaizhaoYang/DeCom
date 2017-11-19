% Performs PhysIO-Regressor generation from Siemens ECG 3T logfile
%
% The physiological ECG logile used is courtesy of
%   Miriam Sebold, Charite Berlin and
%   Quentin Huys, TNU Zurich
%
% IN
%       '2001_PAV.ecg'          ECG logfile from Siemens 3T scanner
% OUT
%       'main_regressors.txt'   file with physiologial regressors
%
%   See also tapas_physio_new tapas_physio_read_physlogfiles_siemens
%
%
% Author: Lars Kasper, adapting similar code for file read-in 
%                      by Miriam Sebold, Charite Berlin (2014)
%
% Created: 2014-08-24
% Copyright (C) 2014 TNU, Institute for Biomedical Engineering, University of Zurich and ETH Zurich.
%
% This file is part of the TAPAS PhysIO Toolbox, which is released under the terms of the GNU General Public
% Licence (GPL), version 3. You can redistribute it and/or modify it under the terms of the GPL
% (either version 3 or, at your option, any later version). For further details, see the file
% COPYING or <http://www.gnu.org/licenses/>.
%
% $Id$

physio      = tapas_physio_new();   % create structure, numbering according to *PhysIO_PhysNoiseBackground.pptx
log_files   = physio.log_files;     % 1a) Read logfiles
preproc     = physio.preproc;
sqpar       = physio.scan_timing.sqpar;
sync        = physio.scan_timing.sync;
model       = physio.model;         % 3)/4) Model physiological time series
verbose     = physio.verbose;       % Auxiliary: Output


%% 1. Define Input Files

log_files.vendor            = 'Siemens';
log_files.cardiac           = 'siemens_PAV.ecg';
log_files.respiration       = '';
log_files.sampling_interval = 1/400; % in seconds;
log_files.relative_start_acquisition = 0 ; % in seconds

%% 2. Define Nominal Sequence Parameter (Scan Timing)

sqpar.Nslices           = 20;
sqpar.NslicesPerBeat    = 20;   % typically equivalent to Nslices; exception: heartbeat-triggered sequence
sqpar.TR                = 2.41; % in seconds
sqpar.Ndummies          = 5;
sqpar.Nscans            = 400;  % ??? guessed
sqpar.onset_slice       = 11;

% Set to >=0 to count scans and dummy
% volumes from beginning of run, i.e. logfile,
% includes counting of preparation gradientssqpar.Nprep             = [];
sqpar.time_slice_to_slice  = sqpar.TR / sqpar.Nslices;



%% 3. Order of RETROICOR-expansions for cardiac, respiratory and
%% interaction terms. Option to orthogonalise regressors

model.type = 'RETROICOR';
model.order = struct('c',3,'r',4,'cr',1, 'orthogonalise', 'none');
model.input_other_multiple_regressors = ''; % either .txt-file or .mat-file (saves variable R)
model.output_multiple_regressors = 'multiple_regressors.txt';

%% 4. Define Gradient Thresholds to Infer Gradient Timing (Philips only)
%
% method to determine slice onset times
% 'nominal' - to derive slice acquisition timing from sqpar directly
% 'gradient' or 'gradient_log' - derive from logged gradient time courses
%                                in SCANPHYSLOG-files (Philips only)
sync.method = 'nominal'; %'gradient_log'; 'nominal'


%% 5. Define which Cardiac Data Shall be Used

preproc.cardiac.modality = 'ECG';
%preproc.cardiac.initial_cpulse_select.method = 'load_from_logfile'; % using
%Siemens cardiac pulse trigger on signals...try out for comparison
preproc.cardiac.initial_cpulse_select.method = 'auto_matched'; % auto detection using cross-correlation to a self-calibrated template
preproc.cardiac.posthoc_cpulse_select.method = 'off';



%% 6. Output Figures to be generated

verbose.level           = 2; % 0 = none; 1 = main plots (default);  2 = debugging plots, for setting up new study; 3 = all plots
verbose.fig_output_file = 'PhysIO.ps'; % Physio.tiff, .ps, .fig possible


%% 7. Run the main script with defined parameters

physio.log_files            = log_files;
physio.preproc              = preproc;
physio.scan_timing.sqpar    = sqpar;
physio.scan_timing.sync     = sync;
physio.model                = model;
physio.verbose              = verbose;

[physio_out, R, ons_secs] = tapas_physio_main_create_regressors(physio);