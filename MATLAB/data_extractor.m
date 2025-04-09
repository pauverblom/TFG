clear;

max_angle_values = [pi];
corona_multiplicity_values = [25];
alpha_values = [0.4];


for m_idx = 1:length(max_angle_values)
    current_max_angle = max_angle_values(m_idx);
    for row = 1:length(corona_multiplicity_values)
        for col = 1:length(alpha_values)
            current_rc = corona_multiplicity_values(row);
            current_alpha = alpha_values(col);
            
            % Load saved figure
            saved_filename = sprintf('corona_%d_angle_%.2f_alpha_%.1f.fig', ...
                current_rc, current_max_angle, current_alpha);
            saved_filepath = fullfile('figures', saved_filename);
            
            % Open the figure (but do not display it)
            fig = openfig(saved_filepath, 'invisible');
            
            % Get all axes in the figure
            axesHandles = findall(fig, 'Type', 'axes');
            
            % Initialize data storage
            filtered_data = {};
            
            % Define the color you are looking for (red)
            target_color = [1, 0, 0];  % RGB for red
            
            % Counter to track which dataset we're on
            scatter_index = 0;
            
            % Loop through each axis to extract plotted data
            for ax = axesHandles'
                % Get all scatter objects in the current axis
                scatters = findall(ax, 'Type', 'Scatter');
            
                % Extract data from scatter plots
                for scatter = scatters'
                    scatter_index = scatter_index + 1;
                    if isequal(scatter.CData, target_color)  % Check if color is red
                        x_data = scatter.XData(:);
                        y_data = scatter.YData(:);
                        filtered_data{end+1} = table(x_data, y_data, 'VariableNames', {'X', 'Y'});
                    end
                end
            end
            
            % Close the figure after extracting data
            close(fig);
            
            % Save filtered data to a CSV file
            csv_filename = strrep(saved_filename, '.fig', '.csv');
            csv_filepath = fullfile('data', csv_filename);
            
            % Combine and save if data exists
            if ~isempty(filtered_data)
                final_table = vertcat(filtered_data{:});
                writetable(final_table, csv_filepath);
                fprintf('Filtered data saved to: %s\n', csv_filepath);
            else
                fprintf('No matching data found in the figure.\n');
            end
        end
    end
end