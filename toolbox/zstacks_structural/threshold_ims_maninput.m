% threshold_ims_maninput: script run by hand while batch_threshold_ims is
%   stopped in the debugger (it uses that function's workspace: ipars,
%   Data, Data_smooth and its local functions). Edit the thresholds below,
%   run section 1 until the mask looks right, then 'dbcont' in the command
%   window. Section 2 is optional and only displays results.

%% 1) Edit manually input_int_ths by looking at the intensity in plotted planes

ipars.input_int_ths = [75 75 75 90 95];

% and re-run
[Data_mask, ipars.int_ths_per_z] = ...
    generate_binary_mask(ipars.input_int_ths, ...
    ipars.z_slices_per_int_ths, Data_smooth);
displayresults(Data_smooth, Data_mask, ipars)

%% 2) display results/save video

pi = [];
pi.range = [0 200];
slice3Dmatrix(double(Data).*Data_mask, pi);

pi.range = [0 500];
% pi.range = [0 10^5];
pi.Y1 = Data_mask;
pi.vgate = 1;
pi.vgate = 'test';
slice3Dmatrix(double(Data), pi);
