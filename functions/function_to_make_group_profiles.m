function [globalprofile] = function_to_make_group_profiles(freesurfer_folder, contrasts, depth_sampled, atlas, nodenr, indir, outdir)
% FUNCTION_TO_MAKE_GROUP_PROFILES Creates group-averaged cortical depth profiles
%
% Inputs:
%   freesurfer_folder - Path to FreeSurfer subject folder
%   contrasts - Cell array of contrast names to process
%   depth_sampled - Vector of cortical depth values
%   atlas - Name of atlas parcellation
%   nodenr - Number of nodes/regions in atlas
%   indir - Input directory containing individual profile files
%   outdir - Output directory for saving group profiles
%
% Outputs:
%   globalprofile - Structure containing global profile data
%
% For each contrast, node and hemisphere:
% 1. Loads individual subject profiles
% 2. Calculates z-scored and raw group profiles
% 3. Computes profile statistics (median, IQR, skewness)
% 4. Creates diagnostic plots and saves results

   hemispheres = {'lh','rh'};
   types = {'','_partial'}; % Raw and partial profiles
   
   % Process raw and partial profiles
   for t = 1:2
     type = types{t};
     
     % Process each contrast
     for con = 1:length(contrasts)
            contrast = contrasts{con};
            profcount = 0;
            
            % Process each node
            for node = 1:nodenr
                
                % Process each hemisphere
                for hem = 1:2
                    hemisphere = hemispheres{hem};
                    
                    % Define output filenames
                    filename_group_profile = [outdir,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,type,'.csv'];
                    filename_group_profile_raw = [outdir,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,type,'_notnormalised.csv'];
                    filename_group_profile_plot = [outdir,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,type,'.png'];
                    filename_group_profile_plot_raw = [outdir,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,type,'_notnormalised.png'];
                    filename_group_profile_troubleshooting_plot = [outdir,'/Group_profiles_troubleshooting_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,type,'.png'];
                    filename_group_medians = [outdir,'/Group_medians_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,type,'.csv'];

                    % Only process if files don't exist
                    if exist(filename_group_profile_plot) == 0 || exist(filename_group_medians) == 0
                        % Find and load all profile files for current node/hemisphere
                        filename_profiles = dir([indir,'/*_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,type,'.csv']);
                        
                        % Process each profile file
                        for f = 1:length(filename_profiles)
                           data = csvread([filename_profiles(f).folder,'/', filename_profiles(f).name]);
                           % Extract chimp ID from filename
                           chimp_number(f) = str2double(filename_profiles(f).name(1:3));
                           
                           % Extract values at each depth
                           for d = 1:length(depth_sampled)
                              idx = find(abs(data(:,1) - depth_sampled(d)) < .0001); % Find matching depth
                              all_profiles(d,f) = data(idx,2);
                           end
                           
                           % Z-score each profile
                           all_profiles_z(:,f) = (all_profiles(:,f) - mean(all_profiles(:,f)) ) ./ std(all_profiles(:,f)); 
                           if std(all_profiles(:,f)) == 0 % Handle zero variance case
                               all_profiles_z(:,f) = zeros(size(all_profiles_z(:,f)));
                           end
                        end

                        % Calculate group z-scored profile
                        curr_profile = nanmedian(all_profiles_z');
                        curr_profile_error = iqr(all_profiles_z');
                        csvwrite(filename_group_profile, [depth_sampled', curr_profile', curr_profile_error']);
                        
                        % Calculate raw group profile
                        curr_profile_nn = nanmedian(all_profiles');
                        curr_profile_nn_error = iqr(all_profiles');
                        csvwrite(filename_group_profile_raw, [depth_sampled', curr_profile_nn', curr_profile_nn_error']);
                        
                        % Calculate profile statistics for each subject
                        regional_median = median(all_profiles);
                        regional_iqr = iqr(all_profiles);
                        
                        % Calculate different skewness measures
                        for b = 1:size(all_profiles,2)

                            % Calculation as discussed with casey
                            regional_skewness(1,b) = function_profile_skewness(all_profiles(:,b));
                            % Comparison parameter from paquola paper
                            regional_skewness_old_paquola(1,b) = skewness(all_profiles(:,b));
                            regional_skewness_my_old(1,b) = old_implementation_function_profile_skewness(depth_sampled, all_profiles(:,b));
                            
                            if isnan(regional_skewness(1,b))
                                all_profiles(:,b)
                            end
                            
                            % Create legend entries with skewness values
                            legend_entry{b} = [num2str(chimp_number(b)), ': myold:', num2str(round(regional_skewness_my_old(1,b),2)) , ...
                                             '; CPold: ', num2str(round(regional_skewness_old_paquola(1,b),3)), ...
                                             '; CPnew: ', num2str(round(regional_skewness(1,b),3))];
                        end
                        
                        % Create diagnostic plots
                        colidx = find(~isnan(mean(all_profiles_z)));
                        reduced_matrix = all_profiles(:,colidx);
                        
                        % Plot individual profiles and summary statistics
                        f = figure(300);
                            set(f, 'position',[200 200 1500 500]);
                            
                        subplot(1,2,1)
                            cmap = colormap(jet(length(colidx)));
                            plot(reduced_matrix)
                            ylabel(contrast)
                            legend(legend_entry(colidx),'Location','NorthEastOutside')
                            
                        subplot(1,2,2)
                            plot(nanmedian(reduced_matrix'),'b')
                            hold on
                            plot(nanmean(reduced_matrix'),'r')
                            ylabel(contrast)
                            legend({'median','mean'})
                            
                        saveas(f,filename_group_profile_troubleshooting_plot);
                        close(f)

                        % Create and save profile plots
                        f = figure(1);
                            set(f, 'visible','off');
                            function_plot_profile(curr_profile, curr_profile_error, depth_sampled, contrast, node);
                            saveas(gcf,filename_group_profile_plot);
                            close(f)
                            
                        f = figure(1);
                            set(f, 'visible','off');
                            function_plot_profile(curr_profile_nn, curr_profile_nn_error, depth_sampled, contrast, node);
                            saveas(gcf,filename_group_profile_plot_raw);
                            close(f)    

                        % Save profile statistics
                        csvwrite(filename_group_medians, [chimp_number', regional_median', regional_iqr', ...
                                                        regional_skewness', regional_skewness_old_paquola']);
                    end
                end
            end
     end
   end
end
