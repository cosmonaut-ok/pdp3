function rho = rho_func(x)

rho = 0;
dr = 2/2000;

if (x < dr)
    rho = 1;
end