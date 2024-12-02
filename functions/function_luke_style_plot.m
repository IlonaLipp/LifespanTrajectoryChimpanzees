function function_luke_style_plot(templ_file, surface_file, valuelabels, valuetype, figure_file_to_save_prefix, projection_file_to_save, display_range)
% FUNCTION_LUKE_STYLE_PLOT Creates surface plots of brain data with consistent styling
%
% Inputs:
%   templ_file - Template MGH file path to get matrix dimensions
%   surface_file - Surface mesh file (GIFTI format) for plotting
%   valuelabels - Data values to plot on surface
%   valuetype - Label for colorbar (e.g. 'R1 (s^{-1})')
%   figure_file_to_save_prefix - Prefix for output figure filenames
%   projection_file_to_save - Path to save projected data as MGH file
%   display_range - 2-element vector specifying colormap range [min max]
%
% Creates two figures showing lateral views (left and right) of the surface data
% with consistent lighting, colormaps and styling

    % Create both left and right lateral views
    rotation_angles = [90 -90];
    for ra = 1:length(rotation_angles)
        rotation_angle = rotation_angles(ra);

        % Load template to get matrix dimensions
        [vol, M, mr_parms, volsz] = load_mgh(templ_file);
        
        % Save data as MGH file
        save_mgh(valuelabels, projection_file_to_save, M, mr_parms);
        cdata = single(valuelabels);

        % Create surface plot
        gg_newmat = gifti(surface_file);
        gg.cdata = cdata';
        figure, plot(gg_newmat,gg);
        
        % Set view angle for lateral perspective
        view([rotation_angle 0])
        
        % Add lighting and material properties
        hfig = light; 
        set(hfig, 'position', [1 1 0.2]); 
        lighting gouraud; 
        material dull
        
        % Set colormap properties
        caxis([display_range(1), display_range(end)]);
        colormap(brewermap(256, '*Spectral'));
        
        % Add and format colorbar
        hfig = colorbar();
        set(hfig,'XTick',display_range,'FontSize',10);
        ylabel(hfig, valuetype,'FontSize',11);
        
        % Set figure properties and save
        set(gcf,'color','w');
        figure_file_to_save = [figure_file_to_save_prefix,'_',num2str(ra),'.jpg'];
        saveas(gcf,figure_file_to_save);
    end
end
