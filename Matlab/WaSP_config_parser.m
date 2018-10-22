function [global_settings,views] = WaSP_config_parser( config_file )

fid = fopen( config_file, 'rb' );
CONFIG = fread(fid,'int32');
fclose(fid);

global_settings.number_of_views = CONFIG(1);
global_settings.YUV_TRANSFORM = CONFIG(2);
global_settings.YUV_RATIO_SEARCH = CONFIG(3);
global_settings.STD_SEARCH = CONFIG(4);

views = [];

jj = 5;

for ii = 1:global_settings.number_of_views
    
    view = [];
    view.r = CONFIG(jj); jj = jj+1;
    view.c = CONFIG(jj); jj = jj+1;
    
    view.xx = double(CONFIG(jj))/100000; jj = jj+1;
    view.yy = double(CONFIG(jj))/100000; jj = jj+1;
    
    view.rate_color = double(CONFIG(jj))/100000; jj = jj+1;
    view.rate_depth = double(CONFIG(jj))/100000; jj = jj+1;
   
    view.Ms = CONFIG(jj); jj = jj+1;
    view.NNt= CONFIG(jj); jj = jj+1;
    
    view.std = double(CONFIG(jj))/100000; jj = jj+1;
    view.minimum_depth = CONFIG(jj); jj = jj+1;
   
    view.number_of_color_references = CONFIG(jj); jj = jj+1;
    
    view.color_references = [];
    
    if view.number_of_color_references>0
       
        view.color_references = CONFIG(jj:jj+view.number_of_color_references-1);
        
        jj = jj+view.number_of_color_references;
        
    end
    
    view.number_of_depth_references = CONFIG(jj); jj = jj+1;
    
    view.depth_references = [];
    
    if view.number_of_depth_references>0
       
        view.depth_references = CONFIG(jj:jj+view.number_of_depth_references-1);
        
        jj = jj+view.number_of_depth_references;
        
    end
    
    view.has_segmentation = CONFIG(jj); jj = jj+1;
    
    views = [views; view];
    
end