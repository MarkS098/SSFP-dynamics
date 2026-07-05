% Reads Bruker's procs file and extracts processing parameters
% relevant for SSFP reconstruction. Syntax:
%
%              parameters=read_ssfp_procs(directory)
%
% Parameters:
%
%    directory - path to the directory containing
%                the procs file
%
% Outputs the following fields of the parameters structure:
%
%    .SI           - processed spectrum size
%
%    .nc_proc      - processing scaling factor
%                    (NC_proc)

function parameters=read_ssfp_procs(directory)

% Read acquisition parameters file
file_name=[directory filesep 'pdata',filesep,'1',filesep,'procs'];
disp(['Bruker processing parameter file:      ' file_name]);
procs=readlines(file_name);


% Reading the SI from the proc file
SI = procs(startsWith(procs,'##$SI='));
parameters.SI = str2double(regexp(SI,'\d+','match'));
if isempty(parameters.SI)
    error('SI not found.');
else
    disp(['SI:            ' int2str(parameters.SI)]);
end

% Reading the scaling factor from the proc file
nc_proc = procs(startsWith(procs,'##$NC_proc='));
parameters.nc_proc = str2double(regexp(nc_proc,'\d+','match'));
if isempty(nc_proc)
    error('scale factor not found.');
else
    disp(['scale factor:            ' int2str(parameters.nc_proc)]);
end




