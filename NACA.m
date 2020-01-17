function coord = NACA(n3,n4,N)
%Half cosin spacing to make sure we have more points around the leading ege
x = (1-cos(linspace(0,pi,N)))/2;
%Maximum thickness
t_max = (10*n3 + n4)/100;

t = 5 * t_max * (0.2969*sqrt(x) - 0.1260.*x - 0.3516.*x.^2 + 0.2843.*x.^3 ...
    - 0.1015.*x.^4);
