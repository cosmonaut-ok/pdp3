function dydx = mat4ode(x,y)

eps0 = 8.85e-12;

dydx = [   y(2)/x
         -x/eps0*rho_func(x) ];