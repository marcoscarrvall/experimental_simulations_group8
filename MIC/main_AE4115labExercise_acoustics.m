%% Aeroacoustic Data Processing Script
% Characterizes propeller noise using physics-based non-dimensional scaling.
% Characterizes broadband and tonal noise components across varying AoA/J.

close all; clear; clc;

%% 1. Setup Paths and Inputs
% Working directory
baseDir = 'C:\Users\SID-DRW\OneDrive\Escritorio\MDO\Assigment\XDSM\experimental_simulations_group8\MIC';
cd(baseDir);

% Decide comparison 1 for J and 0 for AoA
option = 0;
if option == 0
    fnFolder = '.\DATA\alpha_change';
    fn = {'acoustic_test_a.txt'}; 
else
    fnFolder = '.\DATA\J_change';
    fn = {'acoustic_test_J.txt'}; 
end

% Constants
D     = 0.2032;    % Propeller diameter [m]
Nb    = 6;         % Number of blades
fS    = 51.2e3;    % Sampling frequency [Hz]
p_ref = 20e-6;     % Reference pressure [Pa]
N_fft = 2^12;      % Block size for broadband analysis

% Phase averaging settings
phIntp = linspace(0, 2*pi, 361);
phIntp(end) = [];
dPh = 0;

%% 2. Processing Loop
for i = 1:length(fn)
    % Load operating conditions
    AVGpath = fullfile(fnFolder, fn{i});
    AVGdata{i} = load(AVGpath);
    
    opp{i}.DPN    = AVGdata{i}(:,1);
    opp{i}.vInf   = AVGdata{i}(:,7); 
    opp{i}.AoA    = AVGdata{i}(:,13);
    opp{i}.RPS_M1 = AVGdata{i}(:,15);
    
    % Pre-allocate MIC structure fields to avoid errors during plotting
    numRuns = length(opp{i}.DPN);
    MIC{i}.J        = cell(numRuns, 1);
    MIC{i}.f_norm   = cell(numRuns, 1);
    MIC{i}.Pi_noise = cell(numRuns, 1);
    MIC{i}.f_broad  = cell(numRuns, 1);
    MIC{i}.SPSL     = cell(numRuns, 1);
    MIC{i}.yAvg     = NaN(length(phIntp), numRuns);

    for j = 1:numRuns
        % --- Construct filename and load TDMS ---
        runName = fn{i}(1:end-4);
        TDMSpath = fullfile(fnFolder, sprintf('%s_run%d_001.tdms', runName, opp{i}.DPN(j)));
        
        if ~exist(TDMSpath, 'file'), continue; end
        rawData = ReadFile_TDMS(TDMSpath);
        
        p_raw = rawData{1}(:,1);
        oneP  = rawData{1}(:,2:3);
        
        % --- Calculate Advance Ratio J ---
        % Important: J = V_inf / (n * D)
        MIC{i}.J{j} = opp{i}.vInf(j) / (opp{i}.RPS_M1(j) * D);

        % --- Tonal Analysis: Phase-averaging ---
        try
            [MIC{i}.yAvg(:,j),~,~,~,~,~] = phaseAvgData(p_raw, oneP(:,1), fS, opp{i}.RPS_M1(j), 1, phIntp, dPh);
            
            % Non-Dimensional Scaling
            p_avg  = MIC{i}.yAvg(:, j);
            N_pavg = length(p_avg);
            P_fft  = fft(p_avg);
            P_mag  = abs(P_fft / N_pavg);
            P1     = P_mag(1:floor(N_pavg/2)+1);
            P1(2:end-1) = 2 * P1(2:end-1);
            p_rms_f     = P1 / sqrt(2);
            
            T_dummy = 15; % Placeholder thrust
            MIC{i}.Pi_noise{j} = (p_rms_f * (D^2)) / T_dummy;
            MIC{i}.f_norm{j}   = (opp{i}.RPS_M1(j) * (0:floor(N_pavg/2))') / (opp{i}.RPS_M1(j) * Nb);
        catch
            fprintf('Warning: Run %d failed tonal analysis.\n', opp{i}.DPN(j));
        end

        % --- Broadband Analysis ---
        [~, f_broad, ~, ~, Gpp, ~] = fcn_spectrumN_V1(N_fft, 1/fS, p_raw, 2);
        MIC{i}.f_broad{j} = f_broad;
        MIC{i}.SPSL{j}    = 20*log10(sqrt(Gpp/p_ref^2));
    end
end

%% 3. Integrated Visualization
set(0, 'DefaultLineLineWidth', 1.2);

% Plot 1: Broadband SPSL
figure('Name', 'Broadband Analysis', 'Color', 'w'); hold on; grid on;
for j = 1:numRuns
    if ~isempty(MIC{1}.f_broad{j})
        lbl = sprintf('AoA = %.1f^o, J = %.2f', opp{1}.AoA(j), MIC{1}.J{j});
        semilogx(MIC{1}.f_broad{j}, MIC{1}.SPSL{j}, 'DisplayName', lbl);
    end
end
xlabel('Frequency [Hz]'); ylabel('SPSL [dB/Hz]'); legend('Location', 'eastoutside');

% Plot 2: Non-Dimensional Tonal Noise (Collapse Check)

figure('Name', 'Non-Dimensional Tonal Noise', 'Color', 'w'); hold on; grid on;
for j = 1:numRuns
    if ~isempty(MIC{1}.f_norm{j})
        lbl = sprintf('AoA = %.1f^o, J = %.2f', opp{1}.AoA(j), MIC{1}.J{j});
        plot(MIC{1}.f_norm{j}, 20*log10(MIC{1}.Pi_noise{j}), 'DisplayName', lbl);
    end
end
xlabel('Harmonic Order (f / f_{BPF})'); ylabel('20 log_{10}(\Pi_{noise})'); xlim([0 10]);
legend('Location', 'eastoutside');

% Plot 3: Phase-Averaged Signal
figure('Name', 'Phase-Averaged Signal', 'Color', 'w'); hold on; grid on;
for j = 1:numRuns
    if ~all(isnan(MIC{1}.yAvg(:, j)))
        lbl = sprintf('AoA = %.1f^o, J = %.2f', opp{1}.AoA(j), MIC{1}.J{j});
        plot(phIntp, MIC{1}.yAvg(:, j), 'DisplayName', lbl);
    end
end
xlabel('Phase Angle [rad]'); ylabel('Acoustic Pressure [Pa]');
legend('Location', 'eastoutside');