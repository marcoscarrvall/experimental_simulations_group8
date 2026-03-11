close all
clear
clc

%% Setup Paths and Inputs
% Update this to your working directory
cd 'C:\Users\SID-DRW\OneDrive\Escritorio\MDO\Assigment\XDSM\experimental_simulations_group8\MIC'

fnFolder = '.\DATA';
fn = {'propOn_dE000_dR000_J16.txt'}; 
% MATLAB script to reformat acoustic test data
% Input: 34 columns -> Output: 41 columns

filename = [fnFolder '\' 'acoustic_test_1_run31_001.txt'];
output_filename = [fnFolder '\' 'reformatted_test_data.txt'];

raw_data = readmatrix(filename);

% 2. Calculate the mean of all columns to get a single measurement line
% This averages pressures, velocities, and motor performance 
avg_values = mean(raw_data, 1); 

% 3. Initialize the 41-column output row with zeros
summary_row = zeros(1, 41);

% 4. Map the averaged values to the new 41-column structure
summary_row(1)     = 1;                     % 01) Datapoint number
summary_row(2:8)   = avg_values(2:8);       % 02-08) Aero data (deltaP to Reynolds) 
% 09-12) remain 0 (not used)
summary_row(13:14) = avg_values(12:13);     % 13-14) Angles (AoA, AoS) 
summary_row(15:21) = avg_values(14:20);     % 15-21) Motor 1 Primary set 
summary_row(22:28) = avg_values(21:27);     % 22-28) Motor 1 Secondary set 
% 29-35) remain 0 (not used)
summary_row(36:41) = avg_values(29:34);     % 36-41) Timestamp (D, M, Y, H, M, S) 

% 5. Save as a single line in a text file
writematrix(summary_row, output_filename, 'Delimiter', ',');

fprintf('Process complete. Single averaged line saved to %s\n', output_filename);
% Constants
D = 0.4064;         % Propeller diameter [m]
Nb = 2;            % Number of blades
fS = 51.2e3;       % Sampling frequency [Hz]
p_ref = 20e-6;     % Reference pressure [Pa]
N_fft = 2^12;      % Block size for broadband analysis

% Phase averaging settings
phIntp = linspace(0,2*pi,361); 
phIntp(end)=[]; 
dPh = 0; 

%% Processing Loop
for i=1:length(fn)
    % Load operating conditions
    AVGpath    = [fnFolder,'\',fn{i}];
    AVGdata{i} = load(AVGpath);
    
    opp{i}.DPN    = AVGdata{i}(:,1);
    opp{i}.vInf   = AVGdata{i}(:,7);  %
    opp{i}.AoA    = AVGdata{i}(:,13); %
    opp{i}.RPS_M1 = AVGdata{i}(:,15); %
    
    for j=1:length(opp{i}.DPN)
        % Construct filename and load TDMS
        runNo = 1;
        TDMSpath = [fnFolder '\' fn{i}(1:end-4) '_run',num2str(opp{i}.DPN(j)),'_',sprintf('%03.0f',runNo),'.tdms'];
        rawData = ReadFile_TDMS(TDMSpath);
        
        % 1. Extract raw pressure
        p_raw = rawData{1}(:,1); 
        MIC{i}.pMic{j} = p_raw;
        MIC{i}.oneP{j} = rawData{1}(:,2:3);
    
        % 2. Tonal Analysis: Phase-averaging
        [MIC{i}.yAvg(:,j),~,~,~,~,~] = phaseAvgData(p_raw, MIC{i}.oneP{j}(:,1), fS, opp{i}.RPS_M1(j), 1, phIntp, dPh);
        
        % 3. Broadband Analysis: Ensemble-averaged PSD
        % bst = 2 indicates time-series data
        [~, f_broad, ~, df, Gpp, ~] = fcn_spectrumN_V1(N_fft, 1/fS, p_raw, 2);
        
        % Convert to Sound Pressure Spectral Level (SPSL)
        MIC{i}.f_broad{j} = f_broad;
        MIC{i}.SPSL{j} = 20*log10(sqrt(Gpp/p_ref^2));
        
        % 4. Non-Dimensional Tonal Analysis
        p_avg = MIC{i}.yAvg(:, j);
        N_pavg = length(p_avg);
        P_fft = fft(p_avg);
        P_mag = abs(P_fft / N_pavg);
        P1 = P_mag(1:floor(N_pavg/2)+1);
        P1(2:end-1) = 2 * P1(2:end-1);
        p_rms_f = P1 / sqrt(2);
        
        % Scaling logic
        T_dummy = 15; % Replace with actual thrust from performance data
        MIC{i}.Pi_noise{j} = (p_rms_f * (D^2)) / T_dummy;
        MIC{i}.f_norm{j} = (opp{i}.RPS_M1(j) * (0:floor(N_pavg/2))') / (opp{i}.RPS_M1(j) * Nb);
    end
end

%% Integrated Visualization


% Plot 1: Broadband SPSL vs Frequency
figure('Name', 'Broadband Analysis');
hold on; grid on;
for j = 1:length(opp{1}.DPN)
    semilogx(MIC{1}.f_broad{j}, MIC{1}.SPSL{j}, 'DisplayName', ['AoA = ', num2str(opp{1}.AoA(j)), ' deg']);
end
xlabel('Frequency [Hz]'); ylabel('SPSL [dB/Hz]');
title('Broadband Acoustic Spectrum (Ensemble Averaged)');
legend;

% Plot 2: Non-Dimensional Tonal Noise
figure('Name', 'Non-Dimensional Tonal Noise');
hold on; grid on;
for j = 1:length(opp{1}.DPN)
    plot(MIC{1}.f_norm{j}, 20*log10(MIC{1}.Pi_noise{j}), 'DisplayName', ['AoA = ', num2str(opp{1}.AoA(j)), ' deg']);
end
xlabel('f / f_{BPF}'); ylabel('20 log_{10}(\Pi_{noise})');
title('Non-Dimensional Tonal Harmonics');
xlim([0 10]); % Focus on first 10 harmonics
legend;