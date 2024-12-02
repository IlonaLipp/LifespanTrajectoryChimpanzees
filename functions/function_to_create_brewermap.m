function function_to_create_brewermap(surface_file, vol_to_plot, vol_name, title_name, figure_file_to_save_prefix, display_range)
% FUNCTION_TO_CREATE_BREWERMAP Creates brain surface plots with Brewer colormap
%
% Inputs:
%   surface_file - Path to surface file in GIFTI format
%   vol_to_plot - Volume data to plot on surface
%   vol_name - Name of volume data for colorbar label
%   title_name - Title for the plot
%   figure_file_to_save_prefix - Prefix for output figure filenames
%   display_range - 2-element vector specifying colormap range [min max]
%
% Creates surface plots with lateral views of both hemispheres using Brewer
% Spectral colormap, with consistent lighting and styling

    % Convert volume data to single precision for plotting
    gg.cdata = single(vol_to_plot);            

    % Load surface mesh from GIFTI file
    gg_newmat = gifti(surface_file); 

    % Create lateral views (left and right hemispheres)
    rotation_angles = [90 -90];
    for ra = 1:length(rotation_angles)
        rotation_angle = rotation_angles(ra);

        % Create surface plot with title
        figure, plot(gg_newmat,gg);
        title(title_name);
        view([rotation_angle 0]);

        % Set lighting and material properties
        hfig = light; 
        set(hfig, 'position', [1 1 0.2]); 
        lighting gouraud; 
        material dull

        % Set colormap range and style
        caxis([display_range(1), display_range(end)]);
        colormap(brewermap(256, '*Spectral'));
        
        % Add and format colorbar
        hfig = colorbar();
        set(hfig,'XTick',display_range,'FontSize',10);
        ylabel(hfig, vol_name,'FontSize',11);
        
        % Set white background
        set(gcf,'color','w');

        % Save figure with error handling
        figure_file_to_save = [figure_file_to_save_prefix,'_',num2str(ra),'.jpg'];
        try
            saveas(gcf,figure_file_to_save);
        catch
            disp('Figure could not be saved');
        end
    end