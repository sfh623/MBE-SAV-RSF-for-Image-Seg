
function test2 

clear all;clc

f = rand(5);

pa.dt = 0.1;
pa.delta = 1;

x = linspace(-1,1,64);
y = linspace(-1,1,64);
u = (x.^4)*(y.^6)';

f = 24*(x.^0)*(y.^6)' + 2 * 360*(x.^2)*(y.^4)' + 360*(x.^4)*(y.^2)';

f = u + pa.dt*pa.delta*f;
[ res ] = AinverseOpt1(f,pa);

norm(res-u,inf)


x = linspace(-pi,pi,64);
y = linspace(-pi,pi,64);
u = sin(x)'*sin(y);
f = 4*u;
f = u + pa.dt*pa.delta*f;
[res] = AinverseOpt1(f,pa);

err = res-u;
norm(err(:),inf)


