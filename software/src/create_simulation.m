function [reservoir, model] = create_simulation(C)
    % create_simulation(C)
    % see also: read_input
    % creates a subsurface storage simulation problem based on an input
    % file read by read_input function
    % C is a cell array

    % Read the input files
    dim = sscanf(C{contains(C, 'DIMENSION')}, 'DIMENSION %d');
    L = sscanf(C{contains(C, 'Length')}, 'Length %f');
    H = sscanf(C{contains(C, 'THICKNESS')}, 'THICKNESS %f');
    W = sscanf(C{contains(C, 'WIDTH')}, 'WIDTH %f');
    depth = sscanf(C{contains(C, 'DEPTH')}, 'DEPTH %f');
    R = sscanf(C{contains(C, 'RADIUS')}, 'RADIUS %f');
    Nx = sscanf(C{contains(C, 'NX')}, 'NX %d');
    Ny = sscanf(C{contains(C, 'NY')}, 'NY %d');
    Nz = sscanf(C{contains(C, 'NZ')}, 'NZ %d');
    k = sscanf(C{contains(C, 'PERM')}, 'PERM %f');
    phi = sscanf(C{contains(C, 'POROS')}, 'POROS %f');
    k_stim = sscanf(C{contains(C, 'STIM_PERM')}, 'STIM_PERM %f');
    R_stim = sscanf(C{contains(C, 'STIM_RAD')}, 'STIM_RAD %f');
    Nx_stim = sscanf(C{contains(C, 'NX_STIM')}, 'NX_STIM %d');
    S_oil = sscanf(C{contains(C, 'OIL_SAT')}, 'OIL_SAT %f');
    S_water = sscanf(C{contains(C, 'WATER_SAT')}, 'WATER_SAT %f');
    S_gas = sscanf(C{contains(C, 'GAS_SAT')}, 'GAS_SAT %f');
    p0 = sscanf(C{contains(C, 'P_INIT')}, 'P_INIT %f');

    % sscanf(C{contains(C, '')}, ' %s'); % the sample code

end