function startup()
    try
        p = mfilename('fullpath');
        file_name = mfilename;
        current_path = p(1:end-1-length(file_name));
        addpath([current_path '/src']);
        addpath([current_path '/examples']);
        % addpath([current_path '/src/Advection1D']);
        % addpath([current_path '/src/Transport1D']);
        % addpath([current_path '/src/Bulk']);
        % addpath([current_path '/src/Tools']);
        % addpath([current_path '/database']);
        % addpath([current_path '/src/FVTool']);
        disp('Subsurface energy storage started successfully');
    catch 
        error('Something went wrong with the the subsurface energy storage start up. Please download the package again, extract it, and run the startup.m file.'); 
    end