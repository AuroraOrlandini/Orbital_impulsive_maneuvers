function [rv] = plotOrbit_tiniz(kepEl, mu, deltaTh, tiniz, stepTh)

% plotOrbit.m -  Gives back an array of cartesian cordinates for the arc length deltaTh of the orbit described by
% kepEl (from theta = tinz to theta = tiniz + deltaTh). 
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
% tiniz [1x1] initial anomaly [rad] 
%
% OUTPUT:
% rv [6 x n] position and velocity for every anomaly considered [km, km/s]


if nargin == 2
    deltaTh = 2*pi;
    stepTh = deg2rad(1);
    tiniz = 0;
end

if nargin == 3
    stepTh = deg2rad(1);
    tiniz = 0;
end

if nargin == 4
    stepTh = deg2rad(1);
end

if tiniz < 0
    tiniz = tiniz+2*pi;
end

if deltaTh < tiniz
    deltaTh = deltaTh+2*pi;
end

theta_t=tiniz:stepTh:deltaTh;

rv = kep2car_mat(kepEl(1), kepEl(2), kepEl(3), kepEl(4), kepEl(5), theta_t, mu);

