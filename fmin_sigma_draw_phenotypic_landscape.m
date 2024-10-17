%% Solve the variance and map the energy (multistable)
% clear
x1=[3.08400000000000,1.09900000000000,1.09940000000000,1.09530000000000,6.28470000000000,1.14220000000000,1.09940000000000,1.00720000000000,5.12540000000000,3.12370000000000,3.03360000000000,0.164300000000000,1.13710000000000,0];
x2=[3.15950000000000,1.09900000000000,1.09920000000000,1.09570000000000,5.49760000000000,1.14280000000000,1.09920000000000,1.44000000000000,5.12590000000000,3.12410000000000,3.03380000000000,0.164000000000000,1.13660000000000,0.432900000000000];
x3=[3.18550000000000,1.09900000000000,1.09850000000000,1.09580000000000,4.52740000000000,1.14440000000000,1.09850000000000,1.97320000000000,5.12620000000000,3.12460000000000,3.03410000000000,0.164400000000000,1.13630000000000,0.966100000000000];
e=0.06;
a=1.1;
b=2;
k=1;
n=3;
s=0.5;
d1 = [e,e,e,e,e,e,e,e,e,e/(k-b*s^n*(n*x1(10)^(n-1)/(s^n+x1(10)^n)^2)),e/(k-b*s^n*(n*x1(11)^(n-1)/(s^n+x1(11)^n)^2)),e, e, e/(k-a*(n*x1(14)^(n-1)*(s^n+x1(14)^n)-n*x1(14)^n*x1(14)^(n-1))/(s^n+x1(14)^n)^2)];
d2 = [e,e,e,e,e,e,e,e,e,e/(k-b*s^n*(n*x2(10)^(n-1)/(s^n+x2(10)^n)^2)),e/(k-b*s^n*(n*x2(11)^(n-1)/(s^n+x2(11)^n)^2)),e, e, e/(k-a*(n*x2(14)^(n-1)*(s^n+x2(14)^n)-n*x2(14)^n*x2(14)^(n-1))/(s^n+x2(14)^n)^2)];
d3 = [e,e,e,e,e,e,e,e,e,e/(k-b*s^n*(n*x3(10)^(n-1)/(s^n+x3(10)^n)^2)),e/(k-b*s^n*(n*x3(11)^(n-1)/(s^n+x3(11)^n)^2)),e, e, e/(k-a*(n*x3(14)^(n-1)*(s^n+x3(14)^n)-n*x3(14)^n*x3(14)^(n-1))/(s^n+x3(14)^n)^2)];

Z =  -log(1+0.6123*(1/sqrt(2*pi*d1(8))*exp(-(X-x1(8)).^2/(2*d1(8)))).*(1/sqrt(2*pi*d1(5))*exp(-(Y-x1(5)).^2/(2*d1(5)))) + 0.0829*(1/sqrt(2*pi*d2(8))*exp(-(X-x2(8)).^2/(2*d2(8)))).*(1/sqrt(2*pi*d2(5))*exp(-(Y-x2(5)).^2/(2*d2(5)))) + 0.3047*(1/sqrt(2*pi*d3(8))*exp(-(X-x3(8)).^2/(2*d3(8)))).*(1/sqrt(2*pi*d3(5))*exp(-(Y-x3(5)).^2/(2*d3(5)))));

%% plot phenotypic landscape
figure;
surfc(X, Y, Z);
hold on;
scatter3(x1(8), x1(5),0, 100, 'MarkerEdgeColor',[0 .5 .5], 'MarkerFaceColor','r');
hold on;
scatter3(x2(8), x2(5),0, 100, 'MarkerEdgeColor',[0 .6 .6], 'MarkerFaceColor','g');
hold on;
scatter3(x3(8), x3(5),0, 100, 'MarkerEdgeColor',[0 .7 .7], 'MarkerFaceColor','y');
xlabel('CHEK1');
ylabel('P53');

 %% Draw contour
figure;
contourf(X, Y, Z, 20);
hold on;
scatter(x1(8), x1(5), 100, 'MarkerEdgeColor',[0 .5 .5], 'MarkerFaceColor','r');
hold on;
scatter(x2(8), x2(5), 100, 'MarkerEdgeColor',[0 .6 .6], 'MarkerFaceColor','g');
hold on;
scatter(x3(8), x3(5), 100, 'MarkerEdgeColor',[0 .7 .7], 'MarkerFaceColor','y');
hold on;
syms x y 
p=-log(0.556*(1/sqrt(2*pi*d1(8))*exp(-(x-x1(8)).^2/(2*d1(8)))).*(1/sqrt(2*pi*d1(5))*exp(-(y-x1(5)).^2/(2*d1(5)))) + 0.097*(1/sqrt(2*pi*d2(8))*exp(-(x-x2(8)).^2/(2*d2(8)))).*(1/sqrt(2*pi*d2(5))*exp(-(y-x2(5)).^2/(2*d2(5)))) + 0.345*(1/sqrt(2*pi*d3(8))*exp(-(x-x3(8)).^2/(2*d3(8)))).*(1/sqrt(2*pi*d3(5))*exp(-(y-x3(5)).^2/(2*d3(5)))));
Dx=-e*diff(p,x);
Dy=-e*diff(p,y);
% Select draw area
xRange = linspace(0, 3, 20);
yRange = linspace(3, 8, 20);
[X, Y] = meshgrid(xRange, yRange);
% Calculate the probability flux at the grid points
DX = subs(Dx, {x, y}, {X, Y});
DY = subs(Dy, {x, y}, {X, Y});
quiver(X, Y, double(DX), double(DY), 'r');
% 
xlabel('CHEK1');
xlabel('CHEK1');
ylabel('P53');
colorbar;


