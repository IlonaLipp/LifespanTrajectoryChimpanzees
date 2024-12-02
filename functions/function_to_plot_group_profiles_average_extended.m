function function_to_plot_group_profiles_average_extended(freesurfer_folder, contrasts, depth_sampled, atlas, nodenr, outdir_group)
% FUNCTION_TO_PLOT_GROUP_PROFILES_AVERAGE_EXTENDED Creates plots comparing group-averaged cortical depth profiles
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
% 3. Plots profiles for contrasts 2-3 overlaid with different colors
% 4. Saves combined plot as PNG file

   hemispheres = {'lh','rh'};
   
   % Define colors for different contrasts
   plotcols = {'r','b','g'}; % Basic colors
   plotcols = {[123 15 15]./255, [20 15 123]./255, [19 92 45]./255}; % RGB colors
   
   % Process each node (including node 0 for whole brain)
    for node = 0:nodenr
        % Define output filename
        filename_profile_plot = [outdir_group,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_ALL_2.png'];   
        
        % Only create plot if it doesn't exist
        if exist(filename_profile_plot) == 0
            % Set up figure
            f = figure(node+1);
            set(f,'Position',[100 100 1000 800]);
            
            subpl_count = 0;
            % Create 2x2 subplot grid
            for hem = 1:2 % Loop through hemispheres
                for plottype = 1:2 % Loop through raw vs partial profiles
                    subpl_count = subpl_count + 1;
                    
                    % Plot contrasts 2-3 overlaid
                    for con = 2:3
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
                        curr_profile = profdata(:,2); 
                        curr_profile_error = profdata(:,3); 
                        depth_sampled = profdata(:,1);

                        % Create subplot
                        subplot(2,2,subpl_count); 
                        hold on;
                        plot(curr_profile,depth_sampled(end:-1:1),'Color',plotcols{con},'LineWidth',2);
                        hold off
                        
                        % Format axes and labels
                        set(gca,'YTick',[0 max(depth_sampled)]);
                        set(gca,'YTickLabel',{'WM','pial'},'FontSize',12);
                        xtickangle(45);
                        ylabel('cortical depth','FontSize',12);
                        xlabel('z');
                        set(gcf,'color','white');
                        legend(contrasts(2:con),'Location','southwestoutside')
                        xlim([-2 2]); 
                        title(hemisphere)
                    end
                end
            end
            
            % Save and close figure
            saveas(gcf,filename_profile_plot);
            close(f)
        end
    end      
end
