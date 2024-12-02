function global_profile = function_make_global_profile(contrasts, atlas, nodenr, outdir_group)
% FUNCTION_MAKE_GLOBAL_PROFILE Calculates global cortical profiles by averaging across nodes and hemispheres
%
% Inputs:
%   contrasts - Cell array of contrast names to process (e.g. 'R1', 'MTsat', 'R2s')
%   atlas - Name of brain atlas being used
%   nodenr - Number of nodes/regions in the atlas
%   outdir_group - Output directory containing group profile CSV files
%
% Output:
%   global_profile - Struct containing averaged profiles for each contrast
%
% The function reads cortical profiles from CSV files for each node and hemisphere,
% then averages them to create a global profile for each contrast type.

    % Define hemispheres and profile types
    hemispheres = {'lh','rh'};
    types = {'','_partial'};
    
    % Process only non-partial profiles for now (t=1)
    for t = 1
        type = types{t};
        
        % Iterate through contrasts
        for con = 1:length(contrasts)
            contrast = contrasts{con};
            count = 1;
            
            % Load and combine profiles from all nodes and hemispheres
            for node = 1:nodenr
                for hem = 1:2
                    hemisphere = hemispheres{hem};
                    
                    % Construct filename and load profile data
                    filename_group_profile = [outdir_group,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,type,'.csv'];
                    datamat = importdata(filename_group_profile);
                    
                    % Store profile data (2nd column contains profile values)
                    all_profiles(:,count) = datamat(:,2); 
                    count = count + 1;
                end
            end
            
            % Calculate mean profile across all nodes and hemispheres
            global_profile.(contrast) = mean(all_profiles,2);
        end
    end
end