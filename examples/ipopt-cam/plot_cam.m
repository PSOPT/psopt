% Plots the optimal cam shape found by ipopt_cam, reading the radii from the
% file ipopt.m that the program writes on exit. The design arc is mirrored
% about the axis of symmetry and the fixed circular part is closed up.

ipopt;

r = x;

n = length(r);

d_theta = 2.0*pi/(5.0*(n+1));

theta = zeros(n,1);

for i=1:n

   theta(i) = (n-i+1)*d_theta;

end;

polarplot(theta,r,'b','LineWidth',2);

hold on;

theta = -theta;

polarplot(theta,r,'b','LineWidth',2);

hold on;

dtheta = abs(theta(2)-theta(1));

theta_circle = max(abs(theta)):dtheta:(2*pi-max(abs(theta)));

theta_circle = theta_circle(:);

r_circle = 1*ones(length(theta_circle),1);

polarplot(theta_circle,r_circle,'b','LineWidth',2);
