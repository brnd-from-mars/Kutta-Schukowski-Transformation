clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry and flow definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 0.5;
center = [-0.07; 0.07];
%center = [0; 0];
radius = norm(center - [a;0]);

u = 2.5;
AoA = 20.0;

u_inf = u * cosd(AoA);
v_inf = u * sind(AoA);
p_atm = 30.0;
rho = 1.0;
Gamma = -6.0;
%Gamma = -4 * pi * a * v_inf;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate cylinder shape outline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = linspace(0, 2*pi, 30);
circle_x = center(1) + radius * cos(phi);
circle_y = center(2) + radius * sin(phi);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate world mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axis_x = linspace(-2, 2, 500);
axis_y = linspace(-2, 2, 500);
[mesh_x, mesh_y] = meshgrid(axis_x, axis_y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate flow around cylinder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[circle_Phi, circle_u, circle_v, circle_mask] = flow_cylinder(mesh_x, mesh_y, center, radius, u_inf, v_inf, Gamma);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Kutta-Schukowski-Transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[airfoil_x, airfoil_y] = kutta_schukowski_transformation(a, circle_x, circle_y); % get airfoil shape outline
airfoil_mask = double(~inpolygon(mesh_x, mesh_y, airfoil_x, airfoil_y)); % get mask for airfoil shape

[kst_mesh_x, kst_mesh_y] = kutta_schukowski_transformation(a, mesh_x, mesh_y); % transform coordinate mesh to kst
kst_mesh_x = kst_mesh_x .* circle_mask; % disgard coordinates previously inside cylinder
kst_mesh_y = kst_mesh_y .* circle_mask;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate flow around airfoil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kst_airfoil_Phi = circle_Phi; % potential stays constant after transformation
airfoil_Phi = griddata(kst_mesh_x, kst_mesh_y, kst_airfoil_Phi, mesh_x, mesh_y) .* airfoil_mask; % interpolate irregular kst mesh to grid mesh

[airfoil_u, airfoil_v] = gradient(airfoil_Phi, axis_x, axis_y); % velocity equals gradient in potential field 
airfoil_u = airfoil_u .* airfoil_mask; % appy mask for airfoil shape
airfoil_v = airfoil_v .* airfoil_mask;
airfoil_u = airfoil_u .* double(abs(airfoil_u) < 6); % disgard artifacts around airfoil outline
airfoil_v = airfoil_v .* double(abs(airfoil_v) < 6);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove circulation discontinuity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discontinuity due to potential jump in circulation flow
airfoil_x_min = min(airfoil_x); % calculate discontinuity length
disc_len = round(length(axis_x) * (airfoil_x_min - axis_x(1)) / (axis_x(end) - axis_x(1)) - 0.5);

delta_u = airfoil_u(2:end,1:disc_len) - airfoil_u(1:end-1,1:disc_len); % get discontinuity y position
delta_v = airfoil_v(2:end,1:disc_len) - airfoil_v(1:end-1,1:disc_len);
disc = (delta_u > 0.6 * u_inf) | (delta_v > 0.6 * v_inf);

for i = 1:disc_len
    bottom = find(disc(:,i),1,'first'); % get discontinuity start and end
    top = find(disc(:,i),1,'last');
    
    top = min(top, length(axis_y) / 2 + 10); % cap discontinuity position to avoid border effects
    bottom = max(bottom, length(axis_y) / 2 - 10);
    
    if (isempty(top) || isempty(bottom)) % continue if no discontinuity was found
        continue;
    end
    
    lin_u = linspace(airfoil_u(bottom-2,i),airfoil_u(top+2,i),top-bottom+5); % calculate linear gradient to fill discontinuity
    lin_v = linspace(airfoil_v(bottom-2,i),airfoil_v(top+2,i),top-bottom+5);
    
    airfoil_u(bottom-1:top+1,i) = lin_u(2:end-1)'; % overwrite discontinuity
    airfoil_v(bottom-1:top+1,i) = lin_v(2:end-1)';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate pressure field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
circle_p_dyn = (rho * (circle_u.^2 + circle_v.^2) / 2) .* circle_mask;
airfoil_p_dyn = (rho * (airfoil_u.^2 + airfoil_v.^2) / 2) .* airfoil_mask;

circle_p_stat = (p_atm - circle_p_dyn) .* circle_mask;
airfoil_p_stat = (p_atm - airfoil_p_dyn) .* airfoil_mask;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scan airfoil pressure distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
airfoil = [airfoil_x; airfoil_y]; % get segment geometry
segment_dir = [airfoil(:,2:end), airfoil(:,1)] - airfoil;
segment_length = vecnorm(segment_dir);
segment_pos = airfoil + 0.5 * segment_dir;
segment_normal = [segment_dir(2,:); -segment_dir(1,:)] ./ segment_length;

segment_pscan = segment_pos + 0.05 * segment_normal; % Scan pressure field
segment_pressure = interp2(mesh_x, mesh_y, airfoil_p_stat, segment_pscan(1,:), segment_pscan(2,:));
segment_force = - segment_pressure .* segment_length .* segment_normal;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate airfoil lift and CoL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
airfoil_lift = sum(segment_force, 2); % sum forces to get lift

l = length(segment_pos); % calculate moment of segments around origin
segment_moment = cross([segment_pos; zeros(1,l)], [segment_force; zeros(1,l)]);
segment_moment = segment_moment(3,:);

airfoil_moment = sum(segment_moment); % sum moments and calculate moment arm
h = airfoil_moment / norm(airfoil_lift);

airfoil_CoL = h * [airfoil_lift(2); -airfoil_lift(1)] / norm([airfoil_lift(2); -airfoil_lift(1)]); % apply moment arm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot flow diagrams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sgtitle('Demonstration der Kutta-Schukowski-Transformation');

subplot(2,2,1);
plot(circle_x, circle_y);
hold;
contour(mesh_x, mesh_y, circle_Phi);
colorbar;
streamline(mesh_x, mesh_y, circle_u, circle_v, -2*ones(1,30), linspace(-2,2,30));
legend('Zylinder', 'Potentiallinien', 'Stromlinien');
view(2);
axis equal;
shading interp;
xlim([-2,2]);
ylim([-2,2]);
hold;
title('Strom- und Potentiallinien bei Strömung um Zylinder');

subplot(2,2,2);
plot(airfoil_x, airfoil_y);
hold;
contour(mesh_x, mesh_y, airfoil_Phi);
colorbar;
streamline(mesh_x, mesh_y, airfoil_u, airfoil_v, -2*ones(1,30), linspace(-2,2,30));
%streamline(mesh_x, mesh_y, -airfoil_u, -airfoil_v, 2*ones(1,30), linspace(-2,2,30));
%quiver(segment_pos(1,:), segment_pos(2,:), segment_force(1,:), segment_force(2,:));
quiver(airfoil_CoL(1), airfoil_CoL(2), 0.1 * airfoil_lift(1), 0.1 * airfoil_lift(2));
legend('Airfoil', 'Potentiallinien', 'Stromlinien');
view(2);
axis equal;
shading interp;
xlim([-2,2]);
ylim([-2,2]);
hold;
title('Strom- und Potentiallinien bei Strömung um Airfoil');

subplot(2,2,3);
surf(mesh_x, mesh_y, circle_p_stat);
colorbar;
hold;
view(2);
axis equal;
shading interp;
xlim([-2,2]);
ylim([-2,2]);
hold;
title('Statischer Druck um Zylinder');

subplot(2,2,4);
surf(mesh_x, mesh_y, airfoil_p_stat);
colorbar;
hold;
view(2);
axis equal;
shading interp;
xlim([-2,2]);
ylim([-2,2]);
hold;
title('Statischer Druck um Airfoil');

