function function_make_roiwise_surface_plot(variable_vec, variable_name, plot_perc, projection_file_to_save, figure_file_to_save_prefix)
% FUNCTION_MAKE_ROIWISE_SURFACE_PLOT Creates surface plots of ROI-wise brain data
%
% Inputs:
%   variable_vec - Vector of values to plot for each ROI (76 values, 38 per hemisphere)
%   variable_name - Name of variable for colorbar label
%   plot_perc - 2-element vector specifying percentile range for colormap [min_perc max_perc]
%   projection_file_to_save - Path to save projected data as MGH file
%   figure_file_to_save_prefix - Prefix for output figure filenames
%
% Creates surface plots showing ROI values mapped onto brain hemispheres using BB38 atlas,
% with consistent styling and separate views for left and right hemispheres

    % Set paths to required FreeSurfer and atlas files
    freesurfer_folder = 'Cortical_analysis/group_freesurfer/fsaverage';
    atlas_folder = 'LS_FS_atlas/BB38chimp_forSharing_Ilona';
    
    % Calculate display range based on percentiles of non-zero values
    volidx = 1:length(variable_vec);
    display_range = [round(prctile(variable_vec(volidx),plot_perc(1)),2),round(prctile(variable_vec(volidx),plot_perc(2)), 2)];
    
    % Process each hemisphere separately
    hemispheres = {'lh','rh'};
    for hem = 1:2
        
        hemisphere = hemispheres{hem}
        % Split values between hemispheres (1-38 for left, 39-76 for right)
        if hem == 1
            values_to_assign = variable_vec(1:38);
        else
            values_to_assign = variable_vec(39:76);
        end

        % Load atlas annotation and surface files
        [vertices, labeling, colortable] = read_annotation([freesurfer_folder '/label/',hemisphere,'.BB38chimp.annot']);
        surface_file = ([freesurfer_folder '/surf/',hemisphere,'.white']);
        colors = colortable.table(1:end, 5);
        
        % Convert atlas labels to 0-based indices
        clear labeling_corr
        for v = 1:length(labeling)
                labeling_corr(v,1) = find(colors == labeling(v))-1; 
        end

        % Map ROI values to vertices based on labels
        clear val_vec
        for label = 1:length(labeling_corr)
            if labeling_corr(label) > 0
                val_vec(label) = values_to_assign(labeling_corr(label));
            else
                val_vec(label) = 0;
            end 
        end

        % Load template and save projected data
        templ_file = [freesurfer_folder '/surf/',hemisphere,'.white.avg.area.mgh']
        [vol, M, mr_parms, volsz] = load_mgh(templ_file);
        try
            save_mgh(val_vec, projection_file_to_save, M, mr_parms);
        end

        % Prepare surface data for plotting
        gg.cdata = single(val_vec)';
        gg_newmat = gifti(surface_file); 
         
        % Create lateral views (left and right)
        rotation_angles = [90 -90];
        for ra = 1:length(rotation_angles)
            rotation_angle = rotation_angles(ra);

            % Create surface plot with consistent styling
            colormap(brewermap(256, '*Spectral'));
            figure, plot(gg_newmat,gg); 
            
            % Set colormap range, falling back to data range if display_range fails
            try
                caxis([display_range(1), display_range(end)]);
            catch
                display_range = [min(val_vec), max(val_vec)];
                caxis([display_range(1), display_range(end)]);
            end
            
            % Set view angle and lighting
            view([rotation_angle 0])
            hfig = light; 
            set(hfig, 'position', [1 1 0.2]); 
            lighting gouraud; 
            material dull
            
            % Add and format colorbar
            colormap(brewermap(256, '*Spectral'));
            hfig = colorbar();
            try
                set(hfig,'XTick',display_range,'FontSize',10);
            end
            ylabel(hfig,variable_name,'FontSize',11,'interpreter','none');
            set(gcf,'color','w');
            
            % Save figure
            figure_file_to_save = [figure_file_to_save_prefix,'_',num2str(ra),'_',hemisphere,'.jpg'];
            try
                saveas(gcf,figure_file_to_save);
            end
        end
    end
end
