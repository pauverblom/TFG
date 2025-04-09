%% Define parameters
max_angle_values = [0, pi/2, pi];
corona_multiplicity_values = [1, 10, 50];
alpha_values = [0, 0.1, 1];
figures_folder = 'figures';

%% Loop through each max_angle to create combined figures
for m_idx = 1:length(max_angle_values)
    current_max_angle = max_angle_values(m_idx);
    
    % Create combined figure with tiled layout
    combined_fig = figure;
    tcl = tiledlayout(3, 3, 'TileSpacing', 'Compact');
    
    % Iterate over grid positions
    for row = 1:length(corona_multiplicity_values)
        for col = 1:length(alpha_values)
            current_rc = corona_multiplicity_values(row);
            current_alpha = alpha_values(col);
            
            % Load saved figure
            saved_filename = sprintf('corona_%d_angle_%.2f_alpha_%.1f.fig', ...
                current_rc, current_max_angle, current_alpha);
            saved_filepath = fullfile(figures_folder, saved_filename);
            
            if ~exist(saved_filepath, 'file')
                warning('File %s not found. Skipping.', saved_filepath);
                continue;
            end
            
            % Open figure and copy its content
            h_saved = openfig(saved_filepath, 'invisible');
            
            % --- KEY FIX: Directly copy children of the figure ---
            % Find all axes (excluding legends/colorbars)
            ax_src = findobj(h_saved, 'Type', 'axes', '-not', {'Tag', 'Colorbar', '-or', 'Tag', 'Legend'});
            
            if isempty(ax_src)
                warning('No valid axes in %s. Skipping.', saved_filename);
                close(h_saved);
                continue;
            end
            
            % Create target tile and copy content
            nexttile(tcl);  % Focus on current tile
            ax_dst = gca;   % Get handle to destination axes
            
            % Copy all children (lines, surfaces, etc.) from source to destination
            copyobj(ax_src.Children, ax_dst);
            
            % Copy axis properties (limits, labels, etc.)
            axis(ax_dst, 'tight');
            ax_dst.XLim = ax_src.XLim;
            ax_dst.YLim = ax_src.YLim;
            ax_dst.ZLim = ax_src.ZLim;
            ax_dst.XLabel.String = ax_src.XLabel.String;
            ax_dst.YLabel.String = ax_src.YLabel.String;
            ax_dst.Title.String = sprintf('$R_c = %d \\, R_s, \\quad \\theta_{\\max} = %.2f, \\quad \\alpha = %.1f$', current_rc, current_max_angle, current_alpha);
            ax_dst.Title.Interpreter = 'latex';
            
            % --- NEW: Preserve loglog scaling ---
            ax_dst.XScale = ax_src.XScale; % Copy X-axis scale (log/linear)
            ax_dst.YScale = ax_src.YScale; % Copy Y-axis scale (log/linear)
            
            close(h_saved); % Cleanup
        end
    end
end