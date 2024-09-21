function [X,Y,Z] = plotOrbit(kepEl, mu, deltaTh, stepTh)

% plotOrbit.m - Gives back an array of cartesian cordinates for the arc length deltaTh of the orbit described by
% kepEl.
%
% PROTOTYPE:
% plotOrbit(kepEl, mu, deltaTh, stepTh)
%
% DESCRIPTION:
% Plot the arc length of the orbit described by a set of orbital
% elements for a specific arc length.
%
% INPUT:
% kepEl [1x6] orbital elements [km,rad]
% mu [1x1] gravitational parameter [km^3/s^2]
% deltaTh [1x1] arc length [rad]
% stepTh [1x1] arc length step [rad]
%
% OUTPUT:
% X [1xn] X position [km]
% Y [1xn] Y position [km]
% Z [1xn] Z position [km]


if nargin == 2
    deltaTh = 2*pi;
    stepTh = deg2rad(1);
end

if nargin == 3
    stepTh = deg2rad(1);
end


theta_t = 0:stepTh:deltaTh;

rv = kep2car_mat(kepEl(1), kepEl(2), kepEl(3), kepEl(4), kepEl(5), theta_t,mu);

X = rv(1,:);
Y = rv(2,:);
Z = rv(3,:);



