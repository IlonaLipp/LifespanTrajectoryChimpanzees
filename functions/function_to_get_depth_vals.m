function function_to_get_depth_vals(freesurfer_folder, contrasts, depth_sampled, atlas, nodenr, outdir, brain_undersc_id)
% FUNCTION_TO_GET_DEPTH_VALS Extracts and processes cortical depth profiles from FreeSurfer data
%
% Inputs:
%   freesurfer_folder - Path to FreeSurfer subject folder containing depth profiles
%   contrasts - Cell array of contrast names to process
%   depth_sampled - Vector of cortical depth values to sample
%   atlas - Name of atlas parcellation
%   nodenr - Number of nodes/regions in atlas
%   outdir - Output directory for saving profiles and plots
%   brain_undersc_id - Subject ID with underscores
%
% For each contrast, hemisphere, and atlas node:
% 1. Extracts depth profiles from MGH files
% 2. Calculates median values and IQR across vertices
% 3. Saves raw profiles and plots
% 4. Computes global profiles by averaging across nodes
% 5. Calculates partial profiles by regressing out global effects

    hemispheres = {'lh','rh'};
    
    % Process each contrast
    for con = 1:length(contrasts)
        contrast = contrasts{con};

        % Process each node
        for node = 1:nodenr
            
            % Process each hemisphere
            for hem = 1:2
                hemisphere = hemispheres{hem};

                % Define output filenames
                filename_profile = [outdir,'/',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv'];
                filename_profile_plot = [outdir,'/',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.png'];

                clear curr_profile curr_profile_error
                
                % Only process if plot doesn't exist
                if exist(filename_profile_plot) == 0
                
                    % Extract values at each depth
                    for d = 1:length(depth_sampled)
                        depth = depth_sampled(d);                 

                        % Try different filename formats
                        projectionfile = [freesurfer_folder,'/Profiles/equi_',hemisphere,'_',atlas,'-',sprintf('%.3d',node),'_',contrast,'_0p3_',sprintf('%.02f',depth),'.mgh'];
                        if exist(projectionfile) == 0
                              projectionfile = [freesurfer_folder,'/Profiles/equi_',hemisphere,'_',atlas,'-',sprintf('%.3d',node),'_',contrast,'_0p3_',sprintf('%.01f',depth),'.mgh'];   
                        end 
                        
                        % Load and process depth values if file exists
                        if exist(projectionfile) == 2
                            vals = load_mgh(projectionfile);
                            curr_profile(d) = nanmedian(vals(vals~=0));
                            curr_profile_error(d) = iqr(vals(vals~=0));
                        else % If no contrast file available
                            curr_profile(d) = NaN;
                            curr_profile_error(d) = NaN;
                        end
                        
                        % Exclude R1 for specific subjects
                        if strcmp(contrast(1:2),'R1')
                            if strcmp(brain_undersc_id(1:3), '004') || strcmp(brain_undersc_id(1:3), '032') || strcmp(brain_undersc_id(1:3), '033')
                                curr_profile(d) = NaN;
                                curr_profile_error(d) = NaN;
                            end
                        end
                    end

                    % Save profile data and create plot
                    csvwrite(filename_profile, [depth_sampled', curr_profile', curr_profile_error']);
                    f = figure(1);
                    set(f, 'visible','off');
                    function_plot_profile(curr_profile, curr_profile_error, depth_sampled, contrast, node);
                    saveas(gcf,filename_profile_plot);
                    close(f);
                end
            end
         end
    end
     
    % Calculate global profiles across nodes
    for con = 1:length(contrasts)
        contrast = contrasts{con};
        filename_globprofile_plot = [outdir,'/',brain_undersc_id,'_',atlas,'_global_',contrast,'.png'];
        filename_globprofile = [outdir,'/',brain_undersc_id,'_',atlas,'_global_',contrast,'.csv'];
        
        % Collect z-scored profiles from all nodes
        all_profiles.(contrast) = [];
        for node = 1:nodenr
            for hem = 1:2
                hemisphere = hemispheres{hem};
                filename_profile = [outdir,'/',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv'];
                profdata = csvread(filename_profile);
                curr_profile = profdata(:,2);
                curr_profile_z = (curr_profile - mean(curr_profile)) ./ std(curr_profile);
                all_profiles.(contrast) = [all_profiles.(contrast), curr_profile_z];
            end
        end
        
        % Calculate median and IQR of global profile
        global_profile.(contrast) = nanmedian(all_profiles.(contrast), 2);
        global_profile_error.(contrast) = iqr(all_profiles.(contrast)')';

        % Save global profile plot if it doesn't exist
        if exist(filename_globprofile_plot) == 0
            f = figure(1);
            set(f, 'visible','off');
            function_plot_profile(global_profile.(contrast), global_profile_error.(contrast), depth_sampled, contrast, 99);
            saveas(gcf,filename_globprofile_plot);
            close(f);
        end
        csvwrite(filename_globprofile, [depth_sampled', global_profile.(contrast), global_profile_error.(contrast)]);
    end
     
    % Calculate partial profiles by regressing out global effects
    for con = 1:length(contrasts)
        contrast = contrasts{con};
        all_profiles.(contrast) = [];
        for node = 1:nodenr
            for hem = 1:2
                hemisphere = hemispheres{hem};
                filename_profile = [outdir,'/',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv'];
                filename_partprofile = [outdir,'/',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'_partial.csv'];
                filename_partprofile_plot = [outdir,'/',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'_partial.png'];
                
                % Calculate partial profile if files don't exist
                if exist(filename_partprofile) == 0 || exist(filename_partprofile_plot) == 0
                    profdata = csvread(filename_profile);
                    curr_profile = profdata(:,2); 
                    curr_profile_error = profdata(:,3);
                    curr_profile_z = (curr_profile - mean(curr_profile)) ./ std(curr_profile);
                    
                    % Regress out global profile
                    [b, bint, r] = regress(curr_profile_z,[global_profile.(contrast),ones(length(curr_profile),1)]);
                    partial_profile = r;
                    
                    % Save partial profile data and plot
                    csvwrite(filename_partprofile, [depth_sampled', partial_profile, curr_profile_error]);
                    f = figure(1);
                    set(f, 'visible','off');
                    function_plot_profile(partial_profile, 0*curr_profile_error, depth_sampled, contrast, node);
                    saveas(gcf,filename_partprofile_plot);
                    close(f);
                end
            end
        end
    end
end