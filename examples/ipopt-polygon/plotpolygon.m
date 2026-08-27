% Draws the two optimal polygons side by side. It needs both solution files,
% so run the program twice first:
%
%     ./ipopt_polygon 16
%     ./ipopt_polygon 50
%
% which write ipopt16.m and ipopt50.m respectively.

clear;

ipopt16;

n = length(x)/2;

theta = x(n+1:end);
r = x(1:n);

x = [];

for i=1:n,
    x(i) = r(i)*cos(theta(i));
    y(i) = r(i)*sin(theta(i));
end;

x1 = x; y1 = y;

pgon1 = polyshape(x,y);

ipopt50;

n = length(x)/2;

theta = x(n+1:end);
r = x(1:n);

x = [];

for i=1:n,
    x(i) = r(i)*cos(theta(i));
    y(i) = r(i)*sin(theta(i));
end;

x2 = x; y2 = y;

pgon2 = polyshape(x,y);

figure;

subplot(1,2,1);
plot(pgon1);hold on;plot(x1,y1,'+');title('n=16');
subplot(1,2,2);
plot(pgon2);hold on;plot(x2,y2,'+');title('n=50');
