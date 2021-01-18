% Creates a reservoir object that contains all the information of the 
% reservoir to be modeled.
classdef Reservoir
    properties
        name
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
        function res = Reservoir(name, ...
            domain_shape, ...
            domain_length, ...
            domain_width, ...
            domain_thickness, ...
            reservoir_depth, ...
            permeability, ...
            porosity, ...
            oil_saturation, ...
            water_saturation, ...
            gas_saturation, ...
            dykstra_parsons_coef)

            res.name = name;
            res.domain_shape = domain_shape; 
            res.domain_length = domain_length;
            res.domain_width = domain_width;
            res.domain_thickness = domain_thickness;
            res.reservoir_depth = reservoir_depth;
            res.permeability = permeability;
            res.porosity = porosity;
            res.oil_saturation = oil_saturation;
            res.water_saturation = water_saturation;
            res.gas_saturation = gas_saturation;
            res.dykstra_parsons_coef = dykstra_parsons_coef;
        end

    end
end