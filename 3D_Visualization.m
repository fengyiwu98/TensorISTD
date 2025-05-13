%% Create figure window
figure('Units', 'pixels', 'Position', [100 100 703 641]); % Set window properties (units/position)

%% Data input and visualization
A = imread('.\all_result\sequence1\STPA_FCTN\0064.bmp'); % Load BMP image data
mesh(A);                     % Create 3D mesh surface plot

%% Axes configuration
set(gca,...
    'linewidth', 1, ...      % Axis line width 1pt
    'Fontname', 'Times New Roman', ... % Western font typeface
    'FontSize', 27);         % Font size 27pt

%% X/Y/Z-axis ticks configuration
xticks(0:50:250);            % X-axis tick interval 50
xtickformat('%d');           % X-axis integer format
yticks(0:50:250);            % Y-axis tick interval 50
ytickformat('%d');           % Y-axis integer format
zticks(0:50:250);            % Z-axis tick interval 50
ztickformat('%d');           % Z-axis integer format

%% Axes limits and view settings
xlim([0 250]);               % X-axis display range [0-250]
ylim([0 250]);               % Y-axis display range [0-250]
zlim([0 250]);               % Z-axis display range [0-250]
xtickangle(0);               % X-tick labels horizontal alignment
ytickangle(0);               % Y-tick labels horizontal alignment
ztickangle(0);               % Z-tick labels horizontal alignment
view(-37.5,30);              % Set 3D view perspective (azimuth -37.5°, elevation 30°)