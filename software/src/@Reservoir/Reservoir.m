% Creates a reservoir object that contains all the information of the 
% reservoir to be modeled.
classdef Reservoir
    properties
        name = 'Danish North Sea'
        domain_shape
        domain_length
        domain_width
        domain_thickness
        reservoir_depth
        permeability
        porosity
        oil_saturation
        water_saturation
        gas_saturation
        dykstra_parsons_coef
    end

    methods
        function phrm = PhreeqcRM(n_cells, n_threads)
            %{
            phrm = PhreeqcRM(n_cells, n_threads) 
            n_cells: number of reaction cells
            n_threads: number of CPU cores for OpenMP
            %}
            if ~libisloaded('libphreeqcrm')
                loadlibrary('libphreeqcrm','RM_interface_C.h');
            end
			phrm.ncells = n_cells;
			phrm.nthreads = n_threads;
        end
    end
end