function [myskewness raw_data] = function_profile_skewness(curr_profile)
% FUNCTION_PROFILE_SKEWNESS Calculates skewness of a cortical depth profile
%
% Inputs:
%   curr_profile - Vector of profile values across cortical depth
%
% Outputs:
%   myskewness - Skewness value of the profile distribution
%   raw_data - Raw data points used for skewness calculation
%
% The function converts profile values to positive integers and creates a
% distribution by repeating depth positions weighted by the profile values.
% The skewness of this distribution represents the profile's asymmetry.

    % Use position indices rather than actual depth values
    depth_sampled = 1:length(curr_profile);
    
    % Rescale profile to positive integers (scaling doesn't affect skewness)
    % Add offset of 1 and multiply by 1000 to ensure sufficient resolution
    curr_profile_rescaled = round(1000*(min(curr_profile)*-1 + 1 + curr_profile));
    
    try
        % Create distribution by repeating depth positions weighted by profile values
        % This makes the mean correspond to the depth center of gravity
        raw_data = repelem(depth_sampled, curr_profile_rescaled);
        
        % Calculate skewness of the resulting distribution
        myskewness = skewness(raw_data);
    catch
        % Return NaN if calculation fails (typically due to NaN values)
        myskewness = NaN;
    end
end
