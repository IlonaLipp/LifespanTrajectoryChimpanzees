function [all_medians all_medians_skewness all_table all_table_iqr all_table_skewness all_table_skewness_paquola] = function_make_average_projection(contrasts, atlas, nodenr, outdir)
% FUNCTION_MAKE_AVERAGE_PROJECTION Processes group median data from CSV files into tables
%
% Inputs:
%   contrasts - Cell array of contrast names to process
%   atlas - Name of brain atlas being used
%   nodenr - Number of nodes/regions in the atlas
%   outdir - Output directory containing the CSV files
%
% Outputs:
%   all_medians - Struct containing median values for each contrast
%   all_medians_skewness - Struct containing median skewness values
%   all_table - Struct containing full data matrices
%   all_table_iqr - Struct containing IQR values
%   all_table_skewness - Struct containing skewness values
%   all_table_skewness_paquola - Struct containing Paquola skewness values
%
% The function reads CSV files containing group median data for different brain regions
% and contrasts, processes them into organized data structures for further analysis.

    hemispheres = {'lh','rh'};
    
    % Process each contrast
    for con = 1:length(contrasts)
            contrast = contrasts{con};
            count = 1;
            
            % Iterate through hemispheres and nodes
            for hem = 1:2
                for node = 1:nodenr 
                    hemisphere = hemispheres{hem};
                    
                    % Construct filename and read CSV data
                    filename_group_medians = [outdir,'/Group_medians_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv'];
                    vals = csvread(filename_group_medians);
                    
                    % Extract median and skewness values
                    med_val(count) = median(vals(:,2));
                    mean_val(count) = mean(vals(:,2)); 
                    med_skewness(count) = median(vals(:,4));
                    med_skewness_paquola(count) = median(vals(:,5));
                    
                    % Store full matrices, handling potential size mismatches
                    try
                        full_matrix(:, count) = vals(:,2);
                    catch
                        % Silent error handling - consider adding warning or logging
                    end
                    
                    % Store additional metrics
                    full_matrix_iqr(:, count) = vals(:,3);
                    full_matrix_skewness(:, count) = vals(:,4);
                    full_matrix_skewness_paquola(:, count) = vals(:,5);
                    count = count + 1;
                end
            end
            
            % Store processed data in output structures
            all_medians.(contrast) = med_val;
            all_means.(contrast) = mean_val;
            all_table.(contrast) = full_matrix;
            all_table_iqr.(contrast) = full_matrix_iqr;
            all_table_skewness.(contrast) = full_matrix_skewness;
            all_table_skewness_paquola.(contrast) = full_matrix_skewness_paquola;
            all_medians_skewness.(contrast) = med_skewness;
            all_medians_skewness_paquola.(contrast) = med_skewness_paquola;
    end
end