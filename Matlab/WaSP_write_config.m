function CONFIG = WaSP_write_config( global_settings,views,config_file )
% Writes a config file for WaSP jpeg-pleno VM software
%
% 'global_settings' are four integers for which the first one is the most
% important since it defines the number of views to tbe encoded. The rest
% can be zero, although if at low rates you use the fixed weights you might
% want to use the STD_SEARCH and put a one (1) there.
%
% 'views' is a struct with encoding parameters defined for all views
% 
% Contact pekka.astola@tut.fi if you have any questions.

CONFIG = [];

CONFIG(1) = global_settings.number_of_views; 
CONFIG(2) = global_settings.YUV_TRANSFORM; % whether to use an experimental YUV transformation
% safe to use zero hero
CONFIG(3) = global_settings.YUV_RATIO_SEARCH; % for the experimental YUV transformation
CONFIG(4) = global_settings.STD_SEARCH; % whether or not to use a search for optimal 
% standard deviation parameter in fixed weight view merging, requires that
% views(ii).std is >0

for ii = 1:size(views,1)
    
    CONFIG(end+1) = views(ii).r; % row position of the view, for example for view 003_007 this would be 7
    CONFIG(end+1) = views(ii).c; % and this would be 3
    
    CONFIG(end+1) = views(ii).xx*100000; % camera displacement values for the current view,
    % xx horizontal, yy vertical
    % multiplication with a large number
    % to have enough precision in integer range
    % since currently we write all the config parameters as int32
    CONFIG(end+1) = views(ii).yy*100000;
    
    CONFIG(end+1) = views(ii).rate_color*100000;
    CONFIG(end+1) = views(ii).rate_depth*100000;
   
    CONFIG(end+1) = views(ii).Ms; % sparse filter order
    CONFIG(end+1) = views(ii).NNt; % sparse filter neighborhood -NNt:NNt x -NNt:NNt
    % for example, NNt = 3 means -3:3 x -3:3 neighborhood with 49 elements
    % in total
    
    CONFIG(end+1) = views(ii).std*100000; % put for example 10 for views(ii).std 
    % if at very low rate and there is not enough space to store the LS
    % merging coefficients
    CONFIG(end+1) = views(ii).minimum_depth; % for lytro dataset this should be the largest magnitude of the 
    % negative range. has to be a positive number since this number is
    % subsequently subtracted from the depth data. for HDCA dataset
    % this is always zero.
   
    CONFIG(end+1) = views(ii).number_of_color_references; % for view merging
    
    if views(ii).number_of_color_references>0
       
        % color references are referred using a "temporal indexing"
        % meaning, that the first view to be encoded has the index 0, the
        % second has the index 1, etc. when encoding the first view you
        % should always have views(1).number_of_color_references = 0 since
        % we have nothing to predict from.
        CONFIG(end+1:end+views(ii).number_of_color_references) = views(ii).color_references;

    end
    
    CONFIG(end+1) = views(ii).number_of_depth_references;
    
    if views(ii).number_of_depth_references>0
       
        % same logic here as with the color. we use different reference
        % views for depth prediction, since in many cases quality improves
        % if we DON'T predict from intermediate predicted depth views
        CONFIG(end+1:end+views(ii).number_of_depth_references) = views(ii).depth_references;
        
    end
    
    CONFIG(end+1) = views(ii).has_segmentation; %use zero here,
    % placeholder for future experiments with segmentation
    
end

fid = fopen( config_file, 'wb' );
fwrite(fid,CONFIG,'int32');
fclose(fid);


