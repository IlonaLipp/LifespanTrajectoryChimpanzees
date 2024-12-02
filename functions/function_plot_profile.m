function function_plot_profile(curr_profile, curr_profile_error, depth_sampled, contrast, node)
% FUNCTION_PLOT_PROFILE Creates a cortical depth profile plot with error bars
%
% Inputs:
%   curr_profile - Vector of profile values across cortical depth
%   curr_profile_error - Vector of error values for the profile
%   depth_sampled - Vector of sampled depth points
%   contrast - String indicating the contrast type (e.g. 'R1','MTsat','R2s')
%   node - Integer node/region number to display in title
%
% Creates a plot showing profile values from pial surface to white matter,
% with shaded error bars and consistent styling

    hold on;
    % Plot shaded error bars (error scaled by 0.5)
    se = shadedErrorBar(1:length(curr_profile), curr_profile, ...
                       0.5*curr_profile_error, 'b', 'transparent');
    
    % Add line plot of profile values
    plot(curr_profile, 'b-');
    hold off;
    
    % Configure x-axis to show pial and white matter boundaries
    set(gca, 'XTick', [1,length(depth_sampled)]);
    set(gca, 'XTickLabel', {'pial','WM'}, 'FontSize', 15, 'FontWeight', 'Bold');
    xtickangle(45);
    
    % Add labels and title
    xlabel('cortical depth', 'FontSize', 15, 'FontWeight', 'Bold');
    ylabel(contrast);
    title(['node ', sprintf('%.3d', node)]);
    
    % Set white background
    set(gcf, 'color', 'white');
end