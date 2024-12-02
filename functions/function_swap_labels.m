function [value_labels_lh, value_labels_rh] = function_swap_labels(labeling_lh_corr, labeling_rh_corr, vector)
% FUNCTION_SWAP_LABELS Maps values from a vector to vertex labels for left and right hemispheres
%
% Inputs:
%   labeling_lh_corr - Vector of label indices for left hemisphere vertices
%   labeling_rh_corr - Vector of label indices for right hemisphere vertices 
%   vector - Vector of values to map to the vertices
%
% Outputs:
%   value_labels_lh - Vector of mapped values for left hemisphere vertices
%   value_labels_rh - Vector of mapped values for right hemisphere vertices
%
% For each hemisphere, assigns values from the input vector to vertices based on their
% label index. Vertices with label index 0 or less are assigned value 0.

    % Process left hemisphere vertices
    for l = 1:length(labeling_lh_corr)
        idx = labeling_lh_corr(l);
        % Map value from vector if valid label index, otherwise assign 0
        if idx > 0
            value_labels_lh(l) = vector(idx);
        else
            value_labels_lh(l) = 0;
        end
    end
    
    % Process right hemisphere vertices 
    for r = 1:length(labeling_rh_corr)
        idx = labeling_rh_corr(r);
        % Map value from vector if valid label index, otherwise assign 0
        if idx > 0
            value_labels_rh(r) = vector(idx);
        else
            value_labels_rh(r) = 0;
        end
    end
end