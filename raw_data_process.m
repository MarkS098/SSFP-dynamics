clc; close all; clearvars;

% Filepath
data_dir = '/home/mark/NMR Data/XmX processing/XmX raw data/LF_Acetylacetone_24_11_25_DMSO_NaOH_ao'; % main experiment directory
data_mat_name = 'ACAC_peak_A_ao_FA_10';

% Processing parameters
plot_spect = true;
save_data_mat = true;
file_type = 'abs'; % options: abs,fid, 1r, 1i 
dir_num = [10:110]; % experiment folder names, first file is reference
TR_vals = zeros(1,numel(dir_num)); % pre-allocating TR array
peaks = zeros(1,numel(dir_num)); % pre-allocating peak intensity array
noise_var = zeros(1,numel(dir_num));

n_skip = 0; % number of points to remove from the start of the FID
bounds = [28,34]; % maxima calculation boundaries
noise_bound = [10,20]; % in case a peak is mixed with the noise we take the mean value of the noise in its expected region
min_height = 1e5; % minimum intensity for a point in the spectrum to be considered a peak

for j = 1:numel(dir_num)

% Reading the data
    directory = [data_dir, filesep, num2str(dir_num(j))];
    
    % Read acquisition and processing parameters
    parameters = read_ssfp_acqus(directory);
    proc_parameters = read_ssfp_procs(directory);
    
    switch file_type
        case '1r'
            file_ID = fopen([directory,filesep,'pdata',filesep,'1',filesep,'1r']);
            data = fread(file_ID,'int32');
            fclose(file_ID);
            
        case '1i'
            file_ID = fopen([directory,filesep,'pdata',filesep,'1',filesep,'1i']);
            data = fread(file_ID,'int32');
            fclose(file_ID);

        case 'fid'
            data = read_bruker_data(directory,parameters.version);
            data = data(n_skip+1:end);
            spect = abs(fftshift(fft(data)));

        case 'abs'
            file_ID = fopen([directory,filesep,'pdata',filesep,'1',filesep,'1i']);
            data_imag = fread(file_ID,'int32');
            fclose(file_ID);

            file_ID = fopen([directory,filesep,'pdata',filesep,'1',filesep,'1r']);
            data_real = fread(file_ID,'int32');
            fclose(file_ID);

            data = data_real + 1i*data_imag;
            data = abs(data);
    end

    % Acquisition parameters
    TD = parameters.n_acqps; n_points = TD;
    
    if parameters.version >= 4
        
        if mod(n_points,128) ~= 0
            n_points = round(n_points/128)*128;
        end
    
    elseif mod(n_points,256) ~= 0
        n_points = round(n_points/256)*256;
    end

    SFO1 = parameters.carrier_frq;
    O1 = parameters.offset_frq;
    SWH = parameters.sp_width;
    TR_vals(j) = parameters.rep_time*1e3;
    FA = parameters.flip_angle;

    if (strcmp(file_type,'1r') || strcmp(file_type,'1i') || strcmp(file_type,'abs'))
        
        % Processing parameters
        SI = proc_parameters.SI;
        NC_proc = proc_parameters.nc_proc;
        n_skip = 0;
        n_points = SI;

        % Getting the spectrum via FFT and scaling
        spect = flip(data)/(2^NC_proc);
    end

    % Calculating frequency axis
    Hz_axis = linspace(-SWH/2,SWH/2,n_points - n_skip) + O1; % We remove the corresponding number of points from the axis for plotting
    ppm_axis = Hz_axis/(SFO1*1e-6); % Converting the frequency axis to ppm

    % Taking a section from the spectrum for maxima calculations
    sub_axis = ppm_axis(ppm_axis > min(bounds) & ppm_axis < max(bounds));
    sub_spect = spect(ppm_axis > min(bounds) & ppm_axis < max(bounds));

    % Noise variance calculation
    noise_var(j) = var(spect(ppm_axis > min(noise_bound) & ppm_axis < max(noise_bound)));

    % Determining the maxima of the spectrum
    [pks,locs] = findpeaks(sub_spect,sub_axis,"MinPeakHeight",min_height);
    
    if ~isempty(pks)
        peaks(j) = max(pks);
    else
        peaks(j) = mean(sub_spect(ppm_axis > min(noise_bound) & ppm_axis < max(noise_bound)));
    end

    if plot_spect == true
        % Plotting the spectrum with removed points for every dataset
        figure()
        hold on
        title(['TR = ',num2str(TR_vals(j)),' ms'])
        plot(ppm_axis,spect)
        xlabel('\delta^{13} (ppm)')
        set(gca,'XDir','reverse')
    end
end


% Plotting
figure()
hold on
plot(TR_vals, peaks,'--bo','LineWidth',2)
ylabel('Intensity')
xlabel('TR (ms)')
set(gca,'FontSize',12)


% Saving the .mat file in the data directory
if save_data_mat == true
    save([data_dir,filesep,data_mat_name],'peaks','TR_vals','FA')
end