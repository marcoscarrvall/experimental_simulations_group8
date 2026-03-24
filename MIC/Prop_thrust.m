function [resultsTable] = Prop_thrust()
    %% 1. Load Files
    % Load propeller-off (reference) corrected data
    optsOff = detectImportOptions('fully_corrected_data_propoff.csv');
    dataOff = readtable('fully_corrected_data_propoff.csv', optsOff);
    
    % Load propeller-on experimental corrected data
    optsOn = detectImportOptions('fully_corrected_data.csv');
    dataOn = readtable('fully_corrected_data.csv', optsOn);

    %% 2. Setup Baseline Vectors for Interpolation
    % We extract unique AoA values from the prop-off file to build the reference
    [aoa_ref, idx] = unique(dataOff.AoA); 
    cl_ref_vec = dataOff.CL(idx);
    cd_ref_vec = dataOff.CD(idx);

    %% 3. Constants
    S = 0.165;    % Wing reference area [m^2]
    D = 0.2032;   % Propeller diameter [m]
    rho = 1.225;  % Air density [kg/m^3] (Update if Rho is a column in your CSV)

    %% 4. Extract Prop-On variables
    runs    = dataOn.run;
    alphas  = dataOn.AoA;
    vs      = dataOn.V;
    nrotors = dataOn.rpsM2; 
    cl_on   = dataOn.CL;
    cd_on   = dataOn.CD;
    
    % Calculate dynamic pressure for each run
    qs = 0.5 * rho * vs.^2;

    %% 5. Process Each Run
    numRuns = height(dataOn);
    thrusts = zeros(numRuns, 1);
    Js      = zeros(numRuns, 1);

    for j = 1:numRuns
        % Calculate Advance Ratio J
        if nrotors(j) > 0
            Js(j) = vs(j) / (nrotors(j) * D);
        else
            Js(j) = 0;
        end
        
        alpha_rad = deg2rad(alphas(j));
        
        % --- Calculate X_on (Propeller On Body-Axis Force) ---
        % Transforming the CL/CD from the 'on' file to body-axis axial force
        cx_on = cd_on(j) * cos(alpha_rad) - cl_on(j) * sin(alpha_rad);
        x_on  = cx_on * qs(j) * S;

        % --- Calculate X_off (Propeller Off Body-Axis Force) ---
        % Interpolate the baseline coefficients from the 'propoff' file
        cl_base_interp = interp1(aoa_ref, cl_ref_vec, alphas(j), 'linear', 'extrap');
        cd_base_interp = interp1(aoa_ref, cd_ref_vec, alphas(j), 'linear', 'extrap');
        
        % Transforming the interpolated baseline to body-axis axial force
        cx_off = cd_base_interp * cos(alpha_rad) - cl_base_interp * sin(alpha_rad);
        x_off  = cx_off * qs(j) * S;
        
        % --- Calculate Thrust ---
        % T = -(X_on - X_off). Divided by 2 for thrust per propeller.
        thrusts(j) = abs(-(x_on - x_off) / 2); 
    end

    %% 6. Final Table
    resultsTable = table(runs, alphas, vs, Js, thrusts, ...
        'VariableNames', {'Run', 'Alpha_deg', 'V_ms', 'J', 'Thrust_Prop_N'});
    
    % Optional: display the first few rows to verify
    disp('Propeller Thrust Estimation Table (Using Fully Corrected CSVs):');
    disp(head(resultsTable, 20)); 
end