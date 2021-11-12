% flow_cylinder Flow around cylinder in 2-D
%   [Phi, u, v, mask] = flow_cylinder(mesh_x, mesh_y, center, R, u_inf)
%   returns the potential field, velocity field and a mask for the outisde
%   of the cylinder, given a mesh (generated with gridmesh), the dimensions
%   of the cylinder and the flow velocity at infinity.

function [Phi, u, v, mask] = flow_cylinder(mesh_x, mesh_y, center, R, u_inf, v_inf, Gamma)

    M_u = 2 * pi * u_inf * R^2; % dipol moment
    M_v = 2 * pi * v_inf * R^2;
    
    c_mesh_x = mesh_x - center(1); % transform mesh such that center is in the origin
    c_mesh_y = mesh_y - center(2);

    r_sq = c_mesh_x.^2 + c_mesh_y.^2; % distance to cylinder origin squared
    
    mask = double((sqrt(r_sq) - R) > 0); % mask for the outside of the cylinder

    Phi_u = (u_inf * c_mesh_x + M_u * c_mesh_x ./ r_sq / 2 / pi) .* mask; % potential field
    u_u = (u_inf * (1 + R^2 * (c_mesh_y.^2 - c_mesh_x.^2) ./ r_sq.^2)) .* mask; % velocity field
    v_u = (-2 * u_inf * R^2 * c_mesh_y .* c_mesh_x ./ r_sq.^2) .* mask;
    
    Phi_v = (v_inf * c_mesh_y + M_v * c_mesh_y ./ r_sq / 2 / pi) .* mask; % potential field
    u_v = (-2 * v_inf * R^2 * c_mesh_x .* c_mesh_y ./ r_sq.^2) .* mask; % velocity field
    v_v = (v_inf * (1 + R^2 * (c_mesh_x.^2 - c_mesh_y.^2) ./ r_sq.^2)) .* mask;
    
    
    phi = atan2(c_mesh_y, c_mesh_x);
    Phi_Gamma = (Gamma * phi / 2 / pi) .* mask;
    vphi_Gamma = (Gamma / 2 / pi ./ sqrt(r_sq)) .* mask;
    u_Gamma = (-sin(phi) .* vphi_Gamma) .* mask;
    v_Gamma = (cos(phi) .* vphi_Gamma) .* mask;
    
    
    Phi = Phi_u + Phi_v + Phi_Gamma;
    u = u_u + u_v + u_Gamma;
    v = v_u + v_v + v_Gamma;
end
