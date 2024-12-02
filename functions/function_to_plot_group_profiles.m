function function_to_plot_group_profiles(freesurfer_folder, contrasts, depth_sampled, atlas, nodenr, outdir_group)
% FUNCTION_TO_PLOT_GROUP_PROFILES Creates plots comparing cortical depth profiles
%
% Inputs:
%   freesurfer_folder - Path to FreeSurfer subject folder
%   contrasts - Cell array of contrast names to plot
%   depth_sampled - Vector of cortical depth values
%   atlas - Name of atlas parcellation
%   nodenr - Number of nodes/regions in atlas
%   outdir_group - Output directory for saving group profile plots
%
% For each node:
% 1. Creates 2x2 subplot figure comparing left/right hemispheres
% 2. Shows both raw and partial profiles (global effects regressed out)
% 3. Plots profiles for all contrasts overlaid with different colors
% 4. Saves combined plot as PNG file

   hemispheres = {'lh','rh'};
   plotcols = {'r','b','g'}; % Colors for different contrasts
   
   % Process each node (including node 0 for whole brain)
    for node = 0:nodenr
        % Define output filename
        filename_profile_plot = [outdir_group,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_ALL.png'];   
        
        % Only create plot if it doesn't exist
        if exist(filename_profile_plot) == 0
            subpl_count = 1;
            
            % Create 2x2 subplot grid
            for hem = 1:2 % Loop through hemispheres
                for plottype = 1:2 % Loop through raw vs partial profiles
                    % Set up figure
                    f = figure(node+1);
                    set(f,'units','normalized','Position',[0 0 .9 .7]);
                    subplot(2,2,subpl_count); 
                    subpl_count = subpl_count + 1;
                    clear all_profiles;
                    
                    % Plot each contrast
                    for con = 1:3
                        contrast = contrasts{con};
                        hemisphere = hemispheres{hem};
                        
                        % Get appropriate profile data file
                        if plottype == 1
                            filename_group_profile = [outdir_group,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv'];
                        else
                            filename_group_profile = [outdir_group,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'_partial.csv'];
                        end
                        
                        % Load profile data
                        profdata = csvread(filename_group_profile);
                        profs_for_later(:,hem,plottype) = profdata(:,2);
                        curr_profile = profdata(:,2); 
                        curr_profile_error = profdata(:,2);
                        depth_sampled = profdata(:,1);
                        
                        % Store profiles for skewness calculation
                        if node > 0
                            all_profiles(:,con,node) = curr_profile;
                        end
                        
                        % Plot profile line
                        hold on;
                        plot(curr_profile,depth_sampled(end:-1:1),plotcols{con},'LineWidth',2);
                    end
                    hold off
                    
                    % Format axes and labels
                    set(gca,'YTick',[0 max(depth_sampled)]);
                    set(gca,'YTickLabel',{'WM','pial'},'FontSize',15,'FontWeight','Bold');
                    xtickangle(45);
                    ylabel('cortical depth','FontSize',15,'FontWeight','Bold');
                    xlabel('z');
                    set(gcf,'color','white');
                    
                    % Add title with skewness for non-zero nodes
                    if node > 0
                        group_skewness = function_profile_skewness(depth_sampled', mean(all_profiles(:,:,node),2));
                        title(['node ',sprintf('%.3d',node),' sk: ', sprintf('%.04f', group_skewness)]);
                    end
                    
                    % Add legend and set axis limits
                    legend(contrasts,'Location','southwestoutside')
                    if plottype == 1
                        xlim([-2 2]);
                    else
                        xlim([-2 2]); 
                    end
                end
            end
            
            % Save plot
            saveas(gcf,filename_profile_plot);
            close(f)
        end
    end
end
