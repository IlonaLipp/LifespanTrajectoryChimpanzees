function function_to_make_age_grouped_profiles(freesurfer_folder, contrasts, depth_sampled, atlas, nodenr, indir, outdir)
% FUNCTION_TO_MAKE_AGE_GROUPED_PROFILES Creates cortical depth profiles grouped by age
%
% Inputs:
%   freesurfer_folder - Path to FreeSurfer subject folder
%   contrasts - Cell array of contrast names to process
%   depth_sampled - Vector of cortical depth values to sample
%   atlas - Name of atlas parcellation
%   nodenr - Number of nodes/regions in atlas
%   indir - Input directory containing profile data files
%   outdir - Output directory for saving grouped profile plots
%
% For each contrast, node and hemisphere:
% 1. Loads individual subject profiles from CSV files
% 2. Groups subjects into age bins (<6, 6-17, >17 years)
% 3. Calculates median profile and IQR for each age group
% 4. Creates plot showing profiles for all age groups

    % Load chimpanzee age data
    load_chimp_ages
    
    hemispheres = {'lh','rh'};
    
    % Process each contrast
    for con = 1:length(contrasts)
        contrast = contrasts{con};
        profcount = 0;
        
        % Process each node
        for node = 1:nodenr
            % Process each hemisphere
            for hem = 1:2
                hemisphere = hemispheres{hem};
                filename_group_profile_plot = [outdir,'/Group_age_grouped_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.png'];

                % Create plot if it doesn't exist
                if exist(filename_group_profile_plot) == 0 
                    % Find all profile files for current node/hemisphere
                    filename_profiles = dir([indir,'/*_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv']);
                    
                    % Load profile data from each file
                    for f = 1:length(filename_profiles)
                        data = csvread([filename_profiles(f).folder,'/', filename_profiles(f).name]);
                        % Extract chimp ID from filename
                        chimp_number(f) = str2double(filename_profiles(f).name(1:3));
                        
                        % Extract values at each depth
                        for d = 1:length(depth_sampled)
                            idx = find(abs(data(:,1) - depth_sampled(d)) < .0001); % Find matching depth
                            all_profiles(d,f) = data(idx,2);
                        end
                    end
                    
                    % Get ages for all chimps
                    ages_to_group = ages_from_database(chimp_number); 
                    
                    % Create plot
                    f = figure(1);
                    set(f, 'visible','off');
                    hold on;
                    
                    % Plot profiles for each age group
                    colstouse = {'r','g','b'}; % Colors for age groups
                    for group = 1:3
                        % Define age group boundaries
                        if group == 1
                            gridx = find(ages_to_group < 6);         % Young (<6 years)
                        elseif group == 2
                            gridx = find(ages_to_group >= 6 & ages_to_group < 17); % Adolescent (6-17 years)
                        elseif group == 3
                            gridx = find(ages_to_group >= 17);       % Adult (>17 years)
                        end
                        
                        % Calculate median profile and IQR for age group
                        curr_profile = nanmedian(all_profiles(:,gridx)');
                        curr_profile_error = iqr(all_profiles(:,gridx)');
                        
                        % Plot profile line
                        plot(depth_sampled',curr_profile',[colstouse{group},'-'],'LineWidth',3);
                    end
                    hold off;
                    
                    % Format plot
                    set(gca,'XTick',[0.05 0.25 0.5 0.75 max(depth_sampled)]);
                    set(gca,'XTickLabel',{'WM','25%','50%','75%','pial'},'FontSize',16);
                    xtickangle(45);
                    xlabel('cortical depth','FontSize',15,'FontWeight','Bold');
                    ylabel(contrast);
                    title(['node ',sprintf('%.3d',node),' ',hemisphere]);
                    legend({'<6','6-17','17+'},'Location','NorthEastOutside');
                    set(gcf,'color','white');
                    
                    % Save plot
                    saveas(gcf,filename_group_profile_plot);
                    close(f)    
                end
            end
        end
    end
end
