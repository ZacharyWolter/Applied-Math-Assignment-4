%Rate function describing Newton's law of gravitation
%INPUTS:
%t: the current time
%V: the vector of the position and velocity of the planet
% V = [x_p; y_p; dxdt_p; dydt_p]
%orbit_params: a struct describing the system parameters
% orbit_params.m_sun is the mass of the sun
% orbit_params.m_planet is the mass of the planet
% orbit_params.G is the gravitational constant
%OUTPUTS:
%dVdt: a column vector describing the time derivative of V:
% dVdt = [dxdt_p; dydt_p; d2xdt2_p; d2ydt2_p]
function dVdt = gravity_rate_func(t,V,orbit_params)
m_s = orbit_params.m_sun;
m_p = orbit_params.m_planet;
G = orbit_params.G;
x_p = V(1);
y_p = V(2);
r = [x_p, y_p];

F = -(m_p*m_s*G)/(norm(r))^3 * r;
d2rdt2 = F/m_p;

dVdt = [V(3); V(4); d2rdt2(1); d2rdt2(2)];

end