function [value_labels_lh, value_labels_rh] = function_swap_labels(labeling_lh_corr, labeling_rh_corr, vector)
    %%% lh
    for l = 1:length(labeling_lh_corr)
        idx = labeling_lh_corr(l);
        if idx > 0
            value_labels_lh(l) = vector(idx);
        else
            value_labels_lh(l) = 0;
        end
    end
    %%% rh
    for r = 1:length(labeling_rh_corr)
       idx = labeling_rh_corr(r);
        if idx > 0
            value_labels_rh(r) = vector(idx);
        else
            value_labels_rh(r) = 0;
        end
    end
end