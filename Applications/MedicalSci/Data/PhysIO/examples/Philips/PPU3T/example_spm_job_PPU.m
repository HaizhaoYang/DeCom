%-----------------------------------------------------------------------
% Job saved on 06-Jan-2015 10:12:37 by cfg_util (rev $Rev: 6134 $)
% spm SPM - SPM12 (6225)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.tools.physio.save_dir = {''};
matlabbatch{1}.spm.tools.physio.log_files.vendor = 'Philips';
matlabbatch{1}.spm.tools.physio.log_files.cardiac = {'/Users/kasperla/Documents/code/matlab/smoothing_trunk/PhysIOToolbox/examples/Philips/PPU/SCANPHYSLOG.log'};
matlabbatch{1}.spm.tools.physio.log_files.respiration = {'/Users/kasperla/Documents/code/matlab/smoothing_trunk/PhysIOToolbox/examples/Philips/PPU/SCANPHYSLOG.log'};
matlabbatch{1}.spm.tools.physio.log_files.scan_timing = {''};
matlabbatch{1}.spm.tools.physio.log_files.sampling_interval = [];
matlabbatch{1}.spm.tools.physio.log_files.relative_start_acquisition = 0;
matlabbatch{1}.spm.tools.physio.sqpar.Nslices = 32;
matlabbatch{1}.spm.tools.physio.sqpar.NslicesPerBeat = [];
matlabbatch{1}.spm.tools.physio.sqpar.TR = 2;
matlabbatch{1}.spm.tools.physio.sqpar.Ndummies = 0;
matlabbatch{1}.spm.tools.physio.sqpar.Nscans = 180;
matlabbatch{1}.spm.tools.physio.sqpar.onset_slice = 17;
matlabbatch{1}.spm.tools.physio.sqpar.time_slice_to_slice = [];
matlabbatch{1}.spm.tools.physio.sqpar.Nprep = [];
matlabbatch{1}.spm.tools.physio.model.type = 'RETROICOR_HRV_RVT';
matlabbatch{1}.spm.tools.physio.model.order.c = 3;
matlabbatch{1}.spm.tools.physio.model.order.r = 4;
matlabbatch{1}.spm.tools.physio.model.order.cr = 1;
matlabbatch{1}.spm.tools.physio.model.order.orthogonalise = 'none';
matlabbatch{1}.spm.tools.physio.model.input_other_multiple_regressors = {''};
matlabbatch{1}.spm.tools.physio.model.output_multiple_regressors = 'multiple_regressors.txt';
matlabbatch{1}.spm.tools.physio.thresh.scan_timing.gradient_log.grad_direction = 'y';
matlabbatch{1}.spm.tools.physio.thresh.scan_timing.gradient_log.zero = 700;
matlabbatch{1}.spm.tools.physio.thresh.scan_timing.gradient_log.slice = 1800;
matlabbatch{1}.spm.tools.physio.thresh.scan_timing.gradient_log.vol = [];
matlabbatch{1}.spm.tools.physio.thresh.scan_timing.gradient_log.vol_spacing = [];
matlabbatch{1}.spm.tools.physio.thresh.cardiac.modality = 'PPU';
matlabbatch{1}.spm.tools.physio.thresh.cardiac.initial_cpulse_select.auto_template.min = 1;
matlabbatch{1}.spm.tools.physio.thresh.cardiac.initial_cpulse_select.auto_template.file = 'initial_cpulse_kRpeakfile.mat';
matlabbatch{1}.spm.tools.physio.thresh.cardiac.posthoc_cpulse_select.off = struct([]);
matlabbatch{1}.spm.tools.physio.verbose.level = 2;
matlabbatch{1}.spm.tools.physio.verbose.fig_output_file = 'PhysIO_output_level2.jpg';
matlabbatch{1}.spm.tools.physio.verbose.use_tabs = false;
