classdef Wing < handle
    properties
        % Geometric
        panels % 2D matrix of panel mesh
        panel_colours % Array of hex colours for each panel
        rtp % Root thickness percentage
        ttp % Tip thickness percentage
        dihedral % Dihedral angle of wing, deg
        chord % Root chord length, mm
        tip_chord % Tip chord length, mm
        wing % Collection of wing panels and extrusions
        panel_thicknesses
        panel_angles

        % Engine
        engine_panels
        engine_data

        % Landing gear
        lg_panels

        % Structural metrics
        unit_loading
        shear
        moment
        deflection_angle
        deflection

        % Sizing
        taper_ratio = 1
        planform_area = 0
        aspect_ratio = 0
        GMC % Geometric mean chord
        MAC % Mean aerodynamic chord
        MAC_position % Mean aerodynamic chord position
    end
    methods
        function obj = Wing(panels, panel_colours, rtp, ttp, dihedral)
            obj.panels = panels;
            obj.panel_colours = panel_colours;
            obj.rtp = rtp;
            obj.ttp = ttp;
            obj.dihedral = dihedral;

            % Calculate root chord of the wing
            obj.chord = obj.panels{1}(4, 1) - obj.panels{1}(1, 1);
            obj.tip_chord = obj.panels{end}(3, 1) - obj.panels{end}(2, 1);

            % Add the panels into the wing array
            obj.wing = {panels};

            % Perform the advanced Airbus wing sizing procedure to
            % calculate final wing values
            obj.size_wing();
        end

        function obj = size_wing(obj)
            diBB = 0;
            DMX = 0;
            for i=1:length(obj.panels)
                panel = obj.panels{i};

                root_panel_chord = panel(4, 1) - panel(1, 1);
                tip_panel_chord = panel(3, 1) - panel(2, 1);

                panel_span = panel(2, 2) - panel(1, 2);
                panel_le_sweep = atand((panel(2, 1) - panel(1, 1)) / panel_span);
                panel_te_sweep = atand((panel(3, 1) - panel(4, 1)) / panel_span);

                delta_x_qrt_chord = (panel(2, 1) + 0.25*tip_panel_chord) - (panel(1, 1) + 0.25*root_panel_chord);
                panel_qrt_chord_sweep = delta_x_qrt_chord / panel_span;

                panel_taper_ratio = tip_panel_chord / root_panel_chord;
                panel_area = panel_span * ((tip_panel_chord + root_panel_chord) / 2);
                panel_aspect_ratio = panel_span^2 / panel_area;
                
                obj.taper_ratio = obj.taper_ratio * panel_taper_ratio;
                obj.planform_area = obj.planform_area + 2*panel_area;

                diBB = diBB + (root_panel_chord^2 * panel_span/(3*(panel_taper_ratio^2 + panel_taper_ratio + 1)));

                xBarLoc = delta_x_qrt_chord / (3*(1 + 2*panel_taper_ratio)) / (1 + panel_taper_ratio);
                DMX = DMX + panel_area * (xBarLoc + 0.25*root_panel_chord + panel(1, 1));
            end
            obj.aspect_ratio = obj.wingspan^2 / obj.planform_area;
            obj.GMC = obj.planform_area / obj.wingspan;
            obj.MAC = 2 * diBB / obj.planform_area;
            obj.MAC_position = 2 * DMX / obj.planform_area;
        end

        function obj = extrude_wing(obj)
            % Extrude the wing by modelling the top and bottom skin, as well
            % as adding extruded elements (rib-like elements) between both 
            % skins

            % Create copy of center panel and apply root and tip thickness
            % offsets
            top_skin = deal(obj.panels);
            top_skin = apply_panel_offsets(obj, top_skin, false);

            bottom_skin = deal(obj.panels);
            bottom_skin = apply_panel_offsets(obj, bottom_skin, true);

            % Calculate the thickness (z-height) of each panel at the
            % middle span of the panel
            obj.calc_panel_thicknesses();
            obj.calc_panel_angles();

            % Add extruded panels (rib-like elements) to the panels
            % definition
            %
            % Each panel now becomes
            % [x1 y1 z1 i1 j1 k1 u1 v1 w1]
            % [x2 y2 z2 i2 j2 k2 u2 v2 w2]
            % [x3 y3 z3 i3 j3 k3 u3 v3 w3]
            % [x4 y4 z4 i4 j4 k4 u4 v4 w4]
            %
            % Where xyz represents the wing panel, ijk represents the top
            % extrusion and uvw represents the bottom extrusion
            for i=1:length(obj.panels)
                % Get the middle span of the panel
                y_mid = (obj.panels{i}(2, 2) + obj.panels{i}(1, 2)) / 2;


                % Get the x-position at both the leading and trailing edge
                % at the panel's mid span
                [x_le, x_te] = get_x_on_boundaries(obj, obj.panels{i}, y_mid);

                % [i1 j1 k1]
                % [i2 j2 k2]
                % [i3 j3 k3]
                % [i4 j4 k4]
                top_extrusion = [x_le, y_mid, 0;
                                 x_le, y_mid, 0.5*obj.panel_thicknesses(i);
                                 x_te, y_mid, 0.5*obj.panel_thicknesses(i);
                                 x_te, y_mid, 0];
                    
                % [u1 v1 w1]
                % [u2 v2 w2]
                % [u3 v3 w3]
                % [u4 v4 w4]
                bottom_extrusion = [x_le, y_mid, 0;
                                    x_le, y_mid, -0.5*obj.panel_thicknesses(i);
                                    x_te, y_mid, -0.5*obj.panel_thicknesses(i);
                                    x_te, y_mid, 0];
                
                % Add the new panels to the existing wing panel
                obj.panels{i} = cat(2, obj.panels{i}, top_extrusion);
                obj.panels{i} = cat(2, obj.panels{i}, bottom_extrusion);
            end

            % Assign all panels to 'wing' property
            obj.wing = {obj.panels, top_skin, bottom_skin};
        end

        function obj = subdivide_wing(obj, n, uniform)
            % Subdivides the panels into 'n' number of sub-panels
            %
            % y-positions are calculated by doing a simple linear
            % interpolation between the root and tip y-displacement.
            %
            % x-positions are determined by getting the positions on both
            % the leading and trailing edge for the given y-positions
            %
            % A new 2D array is then created containing these new panel
            % definitions
            subdivided_panels = cell(length(obj.panels)*length(n), 1);
            subdivided_panel_colours = strings(1, length(obj.panels)*length(n));

            % If the user has given less than 3 parameters, set 'uniform'
            % to true. This makes this a default value and not required to
            % be passed
            if nargin < 3
                uniform = true;
            end

            % Calculate the number of subdivisions per panel if not
            % uniformly divided
            n_per_panels = n * ones(1, length(obj.panels));

            if uniform
                % Calculate number of subdivisions per panel based on their
                % span constribution to the wing's total semi-span
                n_per_panels = zeros(1, length(obj.panels));
                for i=1:length(obj.panels)
                    panel_span = obj.panels{i}(2, 2) - obj.panels{i}(1, 2);
                    span_contribution = panel_span / obj.semi_span;
                    n_per_panels(i) = round(span_contribution * n);
                end
            end

            count = 1;
            for i=1:length(obj.panels)
                % Get panel boundary information
                panel = obj.panels{i};
                panel_root_y = panel(1, 2);
                panel_tip_y = panel(3, 2);
                panel_span = panel_tip_y - panel_root_y;

                for j=1:n_per_panels(i)
                    % Get inboard and outboard y-positions for the panel
                    y2 = panel_root_y + j * (panel_span / (n_per_panels(i)));
                    y1 = panel_root_y + (j-1) * (panel_span / (n_per_panels(i)));

                    % Get the x-positions of each point based on their
                    % y-positions
                    [x1, x4] = get_x_on_boundaries(obj, panel, y1);
                    [x2, x3] = get_x_on_boundaries(obj, panel, y2);

                    z_i = 0;

                    subdivided_panels{count} = [x1, y1, z_i;
                                                x2, y2, z_i;
                                                x3, y2, z_i;
                                                x4, y1, z_i];
                    subdivided_panel_colours(count) = obj.panel_colours(i);
                    count = count + 1;
                end
            end

            % Update 'panels' and 'wing' with new panels
            obj.panels = subdivided_panels;
            obj.wing{1} = subdivided_panels;
            obj.panel_colours = subdivided_panel_colours;
        end

        function wingspan = wingspan(obj)
            % Calculates the entire wingspan of the wing
            wingspan = 2*obj.wing{1}{end}(2, 2);
        end

        function semi_span = semi_span(obj)
            % Calculates the semi span of the wing
            semi_span = obj.wing{1}{end}(2, 2);
        end
    
        function volume = volume(obj)
            % Calculates and returns the volume of the wing
            volume = 0;
            for i=2:length(obj.panels)
                root_thickness = obj.panel_thicknesses(i-1);
                tip_thickness = obj.panel_thicknesses(i);
                panel_span = obj.panels{i}(2, 2) - obj.panels{i}(1, 2);
                panel_root_chord = obj.panels{i}(4, 1) - obj.panels{i}(1, 1);
                volume = volume + 0.5*(root_thickness + tip_thickness) * panel_span* panel_root_chord;
            end
        end

        function [unit_loading, shear, moment, deflection_angle, deflection] = calc_structural_metrics(obj, W, E)
            y1 = cellfun(@(x) x(1, 2), obj.panels);
            y2 = cellfun(@(x) x(2, 2), obj.panels);
            
            q1 = ((4 * W) / (pi * obj.wingspan)) * sqrt(1 - ((2 * y1) / obj.wingspan).^2);
            q2 = ((4 * W) / (pi * obj.wingspan)) * sqrt(1 - ((2 * y1) / obj.wingspan).^2);

            unit_loading = (q2 + q1) / 2;
            shear = zeros(1, length(obj.panels));
            moment = zeros(1, length(obj.panels));
            deflection_angle = zeros(1, length(obj.panels));
            deflection = zeros(1, length(obj.panels));

            for i=length(shear)-1:-1:1
                engine_load = 0;
                if ~isempty(obj.engine_data) & i == obj.engine_data(1)
                    engine_load = obj.engine_data(end);
                end
                shear(i) = -engine_load + shear(i + 1) - unit_loading(i)*(y2(i) - y1(i));
            end

            for i=length(moment)-1:-1:1
                moment(i) = moment(i + 1) - ((shear(i+1) + shear(i))/2)*(y2(i) - y1(i));
            end

            for i=length(moment)-1:-1:1
                moment(i) = moment(i + 1) - ((shear(i+1) + shear(i))/2)*(y2(i) - y1(i));
            end

            for i=2:length(deflection_angle)
                deflection_angle(i) = deflection_angle(i-1) + 0.5*((moment(i) / (E * obj.second_moment(obj.panels{i}))) + (moment(i) / (E * obj.second_moment(obj.panels{i-1}))))*(y2(i) - y1(i));
            end

            for i=2:length(deflection)
                deflection(i) = deflection(i-1) + ((deflection_angle(i) + deflection_angle(i-1))/2)*(y2(i) - y1(i));
            end

            obj.unit_loading = unit_loading;
            obj.shear = shear;
            obj.moment = moment;
            obj.deflection_angle = deflection_angle;
            obj.deflection;

            obj.calc_panel_angles();
        end

        function I = second_moment(~, panel)
            if size(panel, 2) < 9
                error("Wing must be extruded before calculating second moment of area for given panel")
            end
            I = ((panel(2, 6) - panel(2, 9))^3 * (panel(4, 1) - panel(1, 1))) / 12;
        end


        function obj = add_engine(obj, engine, y_engine)
            % Add engine weight at specified x and y position

            obj.engine_panels = engine.panels;

            y_inboard = cellfun(@(x) x(1, 2), obj.panels);
            y_outboard = cellfun(@(x) x(2, 2), obj.panels);

            panel_idx = find(y_inboard <= y_engine & y_outboard >= y_engine);

            panel_idx = panel_idx(1);

            delta_x = 0;
            delta_y = y_engine;
            delta_z = obj.panels{panel_idx}(2, 9);

            %for i=1:length(obj.engine_panels)
            %    obj.engine_panels{i}(:, 2) = obj.engine_panels{i}(:, 2) + y_panel(panel_idx);
            %    obj.engine_panels{i}(:, 3) = obj.engine_panels{i}(:, 3) + panel_z_delta;
            %end

            obj.engine_data = [panel_idx, delta_x, delta_y, delta_z, engine.mass*9.80665];
        end
    
        function obj = add_landing_gear(obj, lg)
            obj.lg_panels = lg.panels;
        end
    end
    methods (Access = private)
        function panels = apply_panel_offsets(obj, panels, invert)
            % Loop through each panel and apply z-offsets based on a
            % calculated yk value
            for i=1:numel(panels)
                panels{i}(1, 3) = calc_yk(obj, panels{i}(1, 2), invert);
                panels{i}(2, 3) = calc_yk(obj, panels{i}(2, 2), invert);
                panels{i}(3, 3) = calc_yk(obj, panels{i}(3, 2), invert);
                panels{i}(4, 3) = calc_yk(obj, panels{i}(4, 2), invert);
            end
        end

        function zk = calc_yk(obj, y, invert)   
            % Calculates z-offset as linear interpolation between wing
            % root and tip thicknesses
            t0 = obj.rtp * obj.chord;
            t1 = obj.ttp * obj.tip_chord;
            zk = 0.5 * (t0 - ((2*y)/obj.wingspan)*(t0 - t1));

            if invert
                zk = -zk;
            end
        end

        function [x1, x2] = get_x_on_boundaries(~, panel, y)
            % This function calculates the x-value along the leading and
            % trailing edge of a panel for a given y-value
            le_m = ((panel(2, 1) - panel(1, 1)) / (panel(2, 2) - panel(1, 2)));
            le_c = panel(2, 1) - le_m * panel(2, 2);

            te_m = ((panel(3, 1) - panel(4, 1)) / (panel(3, 2) - panel(4, 2)));
            te_c = panel(3, 1) - te_m * panel(3, 2);

            x1 = le_m*y + le_c;
            x2 = te_m*y + te_c;
        end

        function obj = calc_panel_thicknesses(obj)
            % Calculates the "thickness" of each panel's max span, based on
            % linear interpolation of wing root and tip thicknesses
            obj.panel_thicknesses = zeros(1, length(obj.panels));

            tr = obj.rtp * obj.chord;
            tt = obj.ttp * obj.tip_chord;

            for i=1:length(obj.panels)
                y_mid = (obj.panels{i}(2, 2) + obj.panels{i}(1, 2)) / 2;
                obj.panel_thicknesses(1, i) = tr - ((2 * y_mid) / obj.wingspan)*(tr - tt);
            end
        end

        function obj = calc_panel_angles(obj)
            obj.panel_angles = obj.dihedral * ones(1, length(obj.panels));

            if ~isempty(obj.deflection_angle)
                obj.panel_angles = obj.panel_angles + obj.deflection_angle;
            end
        end

        function q_max = q_max(obj, Q)
            % Calculates the maximum elliptic load per unit span given the
            % entire aircraft's weight
            q_max = (4 * Q) / (pi * obj.wingspan);
        end
    
        function q = q_y(obj, Q, y)
            q = obj.q_max(Q)*sqrt(1 - ((2*y) / obj.wingspan)^2);
        end
    end
end