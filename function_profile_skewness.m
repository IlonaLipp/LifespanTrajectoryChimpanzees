
function [myskewness raw_data] = function_profile_skewness(curr_profile)
     depth_sampled = 1:length(curr_profile); %%% do use position in profile rather than actual depth
    % raw_data = repelem(depth_sampled, curr_profile);
     %%% profile needs to be privided in positive integers. different
     %%% scalings do not affect the skewness calculation, so how it is done
     %%% does not really matter
     curr_profile_rescaled = round(1000*(min(curr_profile)*-1 + 1 + curr_profile)); %%% make sure all values are positive
     try
        raw_data = repelem(depth_sampled, curr_profile_rescaled); %%% repeat elements of the profile position weighted by mpms. 
        %%% then mean of distribution will be corresponding to depth centre of gravity 
        myskewness = skewness(raw_data);
     catch %%% most likely happens if values are NaN
        myskewness = NaN;
     end

end
