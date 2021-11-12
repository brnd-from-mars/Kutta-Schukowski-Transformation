% kutta_schukowski_transformation Transforms a set of x-y-coordinates
%   [transformed_x, transformed_y] = kutta_schukowski_transformation(a, x, y)
%   returns the x and y coodinates in world space after the
%   kutta-schukowski-transformation of x and y with the parameter a.

function [transformed_x, transformed_y] = kutta_schukowski_transformation(a, x, y)

    r_sq = x.^2 + y.^2; % distance to origin squared
    
    transformed_x = x .* (r_sq + a^2) ./ r_sq; % transform coordinates
    transformed_y = y .* (r_sq - a^2) ./ r_sq;
end
