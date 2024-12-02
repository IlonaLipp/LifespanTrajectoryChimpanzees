function function_to_plot_group_profiles_average_LR(freesurfer_folder, contrasts, depth_sampled, atlas, nodenr, outdir_group)
% FUNCTION_TO_PLOT_GROUP_PROFILES_AVERAGE_LR Creates plots comparing left/right averaged cortical depth profiles
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
% 1. Loads raw and partial profiles for left and right hemispheres
% 2. Averages profiles across hemispheres
% 3. Creates 1x2 subplot figure showing raw and partial profiles
% 4. Plots profiles for contrasts 2-3 overlaid with different colors
% 5. Saves combined plot as PNG file

   hemispheres = {'lh','rh'};
   
   % Define colors for different contrasts
   plotcols = {'r','b','g'}; % Basic colors
   plotcols = {[123 15 15]./255, [20 15 123]./255, [19 92 45]./255}; % RGB colors
   
   % Process each node (including node 0 for whole brain)
    for node = 0:nodenr
        filename_profile_plot = [outdir_group,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_ALL.png'];   
        
        % Only create plot if it doesn't exist
        if exist(filename_profile_plot) == 0
            % Load profile data for all hemispheres, plot types and contrasts
            clear profs_for_later profs_error_for_later profs_depths_sampled_for_later
            for hem = 1:2
                for plottype = 1:2 % 1=raw profiles, 2=partial profiles
                    for con = 2:3 % Only plot contrasts 2-3
                        contrast = contrasts{con};
                        hemisphere = hemispheres{hem};
                        
                        % Get appropriate profile data file
                        if plottype == 1
                            filename_group_profile = [outdir_group,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv'];
                        else
                            filename_group_profile = [outdir_group,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'_partial.csv'];
                        end
                        
                        % Load and store profile data
                        profdata = csvread(filename_group_profile);
                        curr_profile = profdata(:,2); 
                        curr_profile_error = profdata(:,3); 
                        depth_sampled = profdata(:,1);
                        
                        profs_for_later(:,con,hem,plottype) = curr_profile;
                        profs_error_for_later(:,con,hem,plottype) = curr_profile_error;
                        profs_depths_sampled_for_later(:,con,hem,plottype) = depth_sampled;
                    end
                end
            end
            
            % Create figure
            f = figure(node+1);
            set(f,'Position',[100 100 1200 300]);
            
            % Create 1x2 subplot grid
            subpl_count = 1;
            for plottype = 1:2 % Loop through raw vs partial profiles
                subplot(1,2,subpl_count); 
                subpl_count = subpl_count + 1;
                
                % Plot each contrast
                for con = 1:3
                    % Average profiles across hemispheres
                    curr_profile = mean(squeeze(profs_for_later(:,con,:,plottype)),2);
                    curr_profile_error = mean(squeeze(profs_error_for_later(:,con,:,plottype)),2);
                    depth_sampled = profs_depths_sampled_for_later(:,con,1,plottype);
                    
                    % Plot profile line
                    hold on;
                    plot(curr_profile,depth_sampled(end:-1:1),'Color',plotcols{con},'LineWidth',2);
                end
                hold off
                
                % Format axes and labels
                set(gca,'YTick',[0.05 0.25 0.5 0.75 max(depth_sampled)]);
                set(gca,'YTickLabel',{'WM','25%','50%','75%','pial'},'FontSize',16);
                xtickangle(45);
                ylabel('cortical depth','FontSize',16);
                xlabel('z');
                set(gcf,'color','white');
                legend(contrasts(2:3),'Location','southwestoutside')
                xlim([-2 2]);           
            end
            
            % Save plot
            saveas(gcf,filename_profile_plot);
            close(f)
        end
    end
end
