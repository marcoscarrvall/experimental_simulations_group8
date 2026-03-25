cd('C:\Users\SID-DRW\OneDrive\Escritorio\MDO\Assigment\XDSM\experimental_simulations_group8\MIC');
resultsTable = Prop_thrust();
% 1. Load the data into a table
cd('C:\Users\SID-DRW\OneDrive\Escritorio\MDO\Assigment\XDSM\experimental_simulations_group8');
filename = 'fully_corrected_data.csv';
data = readtable(filename);

% 2. Create the new column
% Let's calculate the Lift-to-Drag ratio: L_D = CL / CD
% (Note: MATLAB tables allow you to access columns using the dot notation)
data.C_T = resultsTable.CT; % Assuming 'CT' is the column name in resultsTable for the thrust coefficient

% 3. Alternatively, if you want a column of constant values or zeros:
% data.NewColumn = zeros(height(data), 1);

% 4. Save the table back to a new CSV (or overwrite the old one)
writetable(data, 'fully_corrected_data.csv');

disp('New column appended and file saved successfully.');