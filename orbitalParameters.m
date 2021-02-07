function [a,e,i,Omega,w,theta] = orbitalParameters(v,r,mu)

h_vec = cross(r,v); % angular momentum --- km^2/s

e_vec = (cross(v,h_vec)/mu) - r/norm(r); % eccentricity vector 
e = norm(e_vec); % eccentricity 

energyMech = (norm(v)^2)/2 - mu/norm(r); % mechanical energy --- kJ

a = -mu/(2*energyMech); % semi-major axis --- km

i = acosd(h_vec(3)/norm(h_vec)); % inclination angle --- rad

n_vec = cross([0,0,1],h_vec); % node vector 

Omega = acosd(n_vec(1) / norm(n_vec)); % right ascension of the ascending node
if n_vec(2) < 0 % if the y-component of the node vector is negative 
    Omega = 360 - Omega; 
end

w = acosd(dot(n_vec,e_vec)/(norm(n_vec)*norm(e_vec))); % argument of periapsis --- rad
if e_vec(3) < 0 % if z-component of eccentricity vector is negative
    w = 360 - w;
end
theta = acosd(dot(e_vec,r)/(norm(e_vec)*norm(r))); % true anomaly --- rad 
if dot(r,v) < 0 % if the position vector dotted with velocity vector is negative
    theta = 360 - theta;
end


end
