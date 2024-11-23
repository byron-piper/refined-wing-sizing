classdef Engine
    properties
        panels = {
            % Front
            [0 0 0 0; 0 0.5 -0.25 0; 0 0.5 -0.75 0; 0 0 -1 0],
            [0 0 0 0; 0 -0.5 -0.25 0; 0 -0.5 -0.75 0; 0 0 -1 0],
            % Back
            [1 0 0 0; 1 0.5 -0.25 0; 1 0.5 -0.75 0; 1 0 -1 0],
            [1 0 0 0; 1 -0.5 -0.25 0; 1 -0.5 -0.75 0; 1 0 -1 0],
            % Left side
            [0 0 0 1; 0 0.5 -0.25 1; 1 0.5 -0.25 1; 1 0 0 1],
            [0 0.5 -0.25 2; 0 0.5 -0.75 2; 1 0.5 -0.75 2; 1 0.5 -0.25 2],
            [0 0.5 -0.75 1; 0 0 -1 2; 1 0 -1 2; 1 0.5 -0.75 2],
            % Right side
            [0 0 -1 1; 0 -0.5 -0.75 1; 1 -0.5 -0.75 1; 1 0 -1 1],
            [0 -0.5 -0.75 2; 0 -0.5 -0.25 2; 1 -0.5 -0.25 2; 1 -0.5 -0.75 2],
            [0 -0.5 -0.25 1; 0 0 0 1; 1 0 0 1; 1 -0.5 -0.25 1]
        };
        mass
        engine_diameter
        engine_length
        cg = [0.5, 0.5, 0.5];
    end
    methods
        function obj = Engine(engine_diameter, engine_length, engine_mass)
            for i=1:length(obj.panels)
                obj.panels{i}(:, 1) = obj.panels{i}(:, 1) * engine_length;
                obj.panels{i}(:, 2) = obj.panels{i}(:, 2) * engine_diameter;
                obj.panels{i}(:, 3) = obj.panels{i}(:, 3) * engine_diameter;
            end

            obj.cg = [obj.cg(1)*engine_length, obj.cg(2)*engine_diameter, obj.cg(3)*engine_diameter];

            obj.engine_diameter = engine_diameter;
            obj.engine_length = engine_length;

            obj.mass = engine_mass;
        end
    end
end