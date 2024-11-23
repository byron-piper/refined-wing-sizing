function [wing_surf] = draw_wing(aircraft)

% Loop through all panels and subpanels that make up wing and draw a
% trisurface
wing = aircraft.wing;
for a=1:numel(wing)
    for b=1:numel(wing{a})
        panel = wing{a}{b};

        panel_angle = 0;
        if ~isempty(aircraft.panel_angles)
            panel_angle = -aircraft.panel_angles(b);
        end

        x = panel(:, 1);
        y = panel(:, 2) * cosd(panel_angle) + panel(:, 3) * sind(panel_angle);
        z = panel(:, 3) * cosd(panel_angle) - panel(:, 2) * sind(panel_angle);
        
        tri_xy = delaunay(x, y);
        wing_surf = trisurf(tri_xy, x, y, z, 'FaceColor', aircraft.panel_colours(b), 'EdgeColor', '#EEEEEE');
        hold on

        if size(panel, 2) > 3

            i = panel(:, 4);
            j = panel(:, 5) * cosd(panel_angle) + panel(:, 6) * sind(panel_angle);
            k = panel(:, 6) * cosd(panel_angle) - panel(:, 5) * sind(panel_angle);

            u = panel(:, 7);
            v = panel(:, 8) * cosd(panel_angle) + panel(:, 9) * sind(panel_angle);
            w = panel(:, 9) * cosd(panel_angle) - panel(:, 8) * sind(panel_angle);

            tri_ik = delaunay(i, k);
            tri_uw = delaunay(u, w);

            trisurf(tri_ik, i, j, k, 'FaceColor', '#AAAAAA', 'EdgeColor', '#888888');
            trisurf(tri_uw, u, v, w, 'FaceColor', '#AAAAAA', 'EdgeColor', '#888888');
        end
    end
end

if ~isempty(aircraft.engine_panels)
    engine_panels = aircraft.engine_panels;

    attachment_panel = aircraft.panels{aircraft.engine_data(1)};
    delta_x = aircraft.engine_data(2);
    delta_y = aircraft.engine_data(3);
    delta_z = aircraft.engine_data(4);

    attachment_panel_angle = aircraft.panel_angles(aircraft.engine_data(1));
    engine_angle = attachment_panel_angle - aircraft.dihedral;

    for i=1:length(engine_panels)
        x = engine_panels{i}(:, 1);
        y = engine_panels{i}(:, 2) * cosd(-engine_angle) + engine_panels{i}(:, 3) * sind(-engine_angle) + ...
            (delta_y*cosd(-engine_angle) + delta_z*sind(-engine_angle));
        z = engine_panels{i}(:, 3) * cosd(-engine_angle) - engine_panels{i}(:, 2) * sind(-engine_angle) + ...
            (delta_z*cosd(-engine_angle) - delta_y*sind(-engine_angle) - delta_z);
        
        if engine_panels{i}(1, 4) == 0
            tri = delaunay(y, z);
        elseif engine_panels{i}(1, 4) == 1
            tri = delaunay(x, y);
        else 
            tri = delaunay(x, z);
        end
        
        trisurf(tri, x, y, z, 'FaceColor', '#0066FF', 'EdgeColor', '#0033AA');
    end
end

if ~isempty(aircraft.lg_panels)
    lg_panels = aircraft.lg_panels;

    for i=1:length(lg_panels)
        x = lg_panels{i}(:, 1);
        y = lg_panels{i}(:, 2);
        z = lg_panels{i}(:, 3);
        
        if lg_panels{i}(1, 4) == 0
            tri = delaunay(y, z);
        elseif lg_panels{i}(1, 4) == 1
            tri = delaunay(x, y);
        else 
            tri = delaunay(x, z);
        end
        
        trisurf(tri, x, y, z, 'FaceColor', '#00FF66', 'EdgeColor', '#00AA33');
    end
end

set(gca, 'YDir', 'reverse')
axis equal padded
xlabel("X (mm)")
ylabel("Y (mm)")
zlabel("Z (mm)")
hold off
end