classdef LandingGear
    properties
        panels = {
            % ~ ~ ~ 0 = front facing
            % ~ ~ ~ 1 = top facing
            % ~ ~ ~ 2 = side facing
            % Main column
            [0 0 0 1; 0 0.1 0 1; 0.1 0.1 0 1; 0.1 0 0 1],
            [0 0 0 0; 0 0 -1 0; 0 0.1 -1 0; 0 0.1 0 0],
            [0 0.1 0 2; 0 0.1 -1 2; 0.1 0.1 -1 2; 0.1 0.1 0 2],
            [0.1 0.1 0 0; 0.1 0.1 -1 0; 0.1 0 -1 0; 0.1 0 0 0],
            [0.1 0 0 2; 0.1 0 -1 2; 0 0 -1 3; 0 0 0 2],
            [0 0 -1 1; 0 0.1 -1 1; 0.1 0.1 -1 1; 0.1 0 -1 1],
            % Struct
            [0.05 0*cosd(135)-0.5*sind(135) -0.5*cosd(135)-0*sind(135)-0.353 0; 
             0.05 0.05*cosd(135)-0.5*sind(135) -0.5*cosd(135)-0.05*sind(135)-0.353 0; 
             0.05 0.05*cosd(135)+0*sind(135) 0*cosd(135)-0.05*sind(135)-0.353 0; 
             0.05 0 -0.353 0],
             % Trolley
             [-0.2 -0.1 -1 1; -0.2 0.2 -1 1; 0.3 0.2 -1 1; 0.3 -0.1 -1 1],
             % Wheels
             [-0.2-0.15*sqrt(2) 0.2 -1 2; 
              -0.2 0.2 -1-0.15*sqrt(2) 2;
              -0.2+0.15*sqrt(2) 0.2 -1 2;
              -0.2 0.2 -1+0.15*sqrt(2) 2],
             [0.3-0.15*sqrt(2) 0.2 -1 2; 
              0.3 0.2 -1-0.15*sqrt(2) 2;
              0.3+0.15*sqrt(2) 0.2 -1 2;
              0.3 0.2 -1+0.15*sqrt(2) 2],
             [-0.2-0.15*sqrt(2) -0.1 -1 2; 
              -0.2 -0.1 -1-0.15*sqrt(2) 2;
              -0.2+0.15*sqrt(2) -0.1 -1 2;
              -0.2 -0.1 -1+0.15*sqrt(2) 2],
             [0.3-0.15*sqrt(2) -0.1 -1 2; 
              0.3 -0.1 -1-0.15*sqrt(2) 2;
              0.3+0.15*sqrt(2) -0.1 -1 2;
              0.3 -0.1 -1+0.15*sqrt(2) 2],
        };
    end
    methods
        function obj = LandingGear(lg_height)
            for i=1:length(obj.panels)
                obj.panels{i}(:, 1) = obj.panels{i}(:, 1) * lg_height;
                obj.panels{i}(:, 2) = obj.panels{i}(:, 2) * lg_height;
                obj.panels{i}(:, 3) = obj.panels{i}(:, 3) * lg_height;
            end
        end
    end
end