function [all_medians all_medians_skewness all_table all_table_iqr all_table_skewness all_table_skewness_paquola] = function_make_average_projection(contrasts, atlas, nodenr, outdir)
    %%% loads informatino from Group_medians file and puts into tables
    
    hemispheres = {'lh','rh'};
     for con = 1:length(contrasts)
            contrast = contrasts{con};
            count = 1;
            for hem = 1:2
                for node = 1:nodenr 
                    hemisphere = hemispheres{hem};
                    filename_group_medians = [outdir,'/Group_medians_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv']
                    vals = csvread(filename_group_medians);
                    med_val(count) = median(vals(:,2));
                    mean_val(count) = median(vals(:,2));
                    med_skewness(count) = median(vals(:,4));
                    med_skewness_paquola(count) =  median(vals(:,5));
                    try
                    full_matrix(:, count) = vals(:,2);
                    catch
                       1+ 1 
                    end
                    
                    full_matrix_iqr(:, count) = vals(:,3);
                    full_matrix_skewness(:, count) = vals(:,4);
                    full_matrix_skewness_paquola(:, count) = vals(:,5);
                    count = count + 1;
                end
            end
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