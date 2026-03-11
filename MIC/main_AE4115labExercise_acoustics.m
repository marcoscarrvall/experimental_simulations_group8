close all
clear
clc

cd 'C:\Users\SID-DRW\OneDrive\Escritorio\MDO\Assigment\XDSM\experimental_simulations_group8\MIC'
%% Inputs
% path to folder containing the measurement data
fnFolder = '.\DATA';

% structure of filenames to main txt file containing the averaged data - you can add multiple filenames here
fn = {'propOn_dE000_dR000_J16.txt'}; 

% propeller diameter (used to compute advance ratio)
D = 0.4064;

% inputs for phase averaging
phIntp = linspace(0,2*pi,361); % azimuthal grid for phase averaging 
phIntp(end)=[]; % remove the grid point at 2*pi - phase-averaged signal is periodic so should be the same as at 0 deg
dPh = 0; % offset in phase between measurement data and 1P signal data
fS = 51.2e3; % sampling frequency [Hz]

%% Loop over all TDMS files of name "Measurement_i.tdms)" in the specified folder
for i=1:length(fn)
   
    % load data operating file
    AVGpath    = [fnFolder,'\',fn{i}];
    AVGdata{i} = load(AVGpath);
    
    opp{i}.DPN    = AVGdata{i}(:,1);
    opp{i}.vInf   = AVGdata{i}(:,7); % freestream velocity [m/s]
    opp{i}.AoA    = AVGdata{i}(:,13);  % angle of attack [deg]
    opp{i}.AoS    = AVGdata{i}(:,14);  % angle of sideslip [deg]
    opp{i}.RPS_M1 = AVGdata{i}(:,15);  % RPS motor 1 [Hz]
    opp{i}.RPS_M2 = AVGdata{i}(:,22);  % RPS motor 2 [Hz]
    opp{i}.J_M1   = opp{i}.vInf./(opp{i}.RPS_M1*D); % advance ratio motor 1
    opp{i}.J_M2   = opp{i}.vInf./(opp{i}.RPS_M2*D); % advance ratio motor 2
    
    % load microphone data
    for j=1:length(opp{i}.DPN) % loop over all the datapoints for this configuration
        
        % Construct filename (required in case of duplicate files)
        runNo = 1;
        TDMSpath = [fnFolder '\' fn{i}(1:end-4) '_run',num2str(opp{i}.DPN(j)),'_',sprintf('%03.0f',runNo),'.tdms'];
        
        % load data
        rawData = ReadFile_TDMS(TDMSpath);
        disp(['Loaded file ' TDMSpath])
        
        % Extract the microphone pressure and write to 
        MIC{i}.pMic{j} = rawData{1}(:,1); % the data are stored with calibration factor applied [Pa]
        MIC{i}.oneP{j} = rawData{1}(:,2:3); % these are the data from the one-per-revolution trigger (pulse whenever the propellers pass a predefined azimuthal position)
    
        % Perform phase-averaging of signals (optional)
        [MIC{i}.yAvg(:,j),~,~,~,~,~] = phaseAvgData(MIC{i}.pMic{j},MIC{i}.oneP{j}(:,1),fS,opp{i}.RPS_M1(j),1,phIntp,dPh);

    end
    
end % end while loop over files


%% Plot phase-averaged data
figure
for j=1:2
subplot(2,1,j)
plot(phIntp,MIC{1}.yAvg(:,j))
title(['AoA=',sprintf('%.1f',opp{1}.AoA(j)),'deg'])
xlabel('Phase angle [rad]')
ylabel('Acoustic pressure [Pa]')
end


figure,plot((opp{i}.AoA),rms(MIC{i}.yAvg),'*b')
xlabel('AoA [deg]')
ylabel('pRMS tonal content [Pa]')

%% Tonal Noise Analysis via FFT & Non-Dimensional Scaling
% Append this to the end of main_AE4115labExercise_acoustics.m

Nb = 2; % Define the number of blades for your propeller (update if 3 or more)

% Initialize a new figure for the Pi_noise spectra
figure('Name', 'Non-Dimensional Acoustic Spectra');
hold on; grid on;

for i = 1:length(fn)
    for j = 1:length(opp{i}.DPN)
        
        % 1. Extract the phase-averaged pressure for the current run
        % Using column 1 assuming microphone 1 is the primary focus
        % 1. Extract the phase-averaged pressure for the current run
        % Use index 'j' to correctly grab the data for the current operating point
        p_avg = MIC{i}.yAvg(:, j);
        N = length(p_avg);
        
        % 2. Perform the Fast Fourier Transform (FFT)
        P_fft = fft(p_avg);
        
        % 3. Compute single-sided RMS pressure spectrum
        % Divide by N to normalize, take absolute value for magnitude
        P_mag = abs(P_fft / N);
        
        % Extract the single-sided spectrum (up to the Nyquist frequency)
        P1 = P_mag(1:floor(N/2)+1);
        
        % Double the non-DC and non-Nyquist bins to conserve spectral energy
        P1(2:end-1) = 2 * P1(2:end-1);
        
        % Convert amplitude peak to RMS 
        p_rms_f = P1 / sqrt(2);
        
        % 4. Normalize the frequency axis
        rps = opp{i}.RPS_M1(j); % Rotational speed [Hz]
        f_BPF = rps * Nb;       % Blade Passing Frequency [Hz]
        
        % Since p_avg represents exactly one revolution, 
        % the FFT bins are exact harmonics of the rotation rate (1P, 2P, etc.)
        f_axis = rps * (0:floor(N/2))'; 
        f_norm = f_axis / f_BPF; % Normalized to f/f_BPF
        
        % 5. Non-Dimensional Scaling (Pi_noise)
        % IMPORTANT: Replace T_dummy with your actual Thrust data array
        T_dummy = 15; % [N] Placeholder thrust value 
        
        Pi_noise = (p_rms_f * (D^2)) / T_dummy;
        
        % Store the processed data back into the structure
        MIC{i}.f_norm{j} = f_norm;
        MIC{i}.Pi_noise{j} = Pi_noise;
        MIC{i}.p_rms_f_spectrum{j} = p_rms_f;
        
        % Plot the result (using a log scale for the y-axis, standard for acoustics)
        plot(f_norm, 20*log10(Pi_noise), 'DisplayName', ['AoA = ', num2str(opp{i}.AoA(j)), '^\circ']);
    end
end

xlabel('Normalized Frequency ($f / f_{BPF}$)', 'Interpreter', 'latex');
ylabel('$20 \log_{10}(\Pi_{noise})$ [dB]', 'Interpreter', 'latex');
title('Non-Dimensional Tonal Noise Spectra');
xlim([0 100]); % Zoom in on the first 10 BPF harmonics
legend('Location', 'best');