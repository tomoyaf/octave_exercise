function rossler_equation_euler_and_rk4()
    a = 0.2
    b = 0.2
    c = 5.7
    
    x0 = 0.1
    y0 = 0.2
    z0 = 0.3
    
    dt = 0.05
    N = 10000
    
    [X1, Y1, Z1] = runge_kutta_method(a, b, c, x0, y0, z0, dt, N);
    [X2, Y2, Z2] = euler_method(a, b, c, x0, y0, z0, dt, N);
    
    T = 1 : N;
    
    figure;
    plot3(X1, Y1, Z1, X2, Y2, Z2);
    legend('RK4', 'Euler method', 'location', 'east');
end

function [X, Y, Z] = euler_method(a, b, c, x0, y0, z0, dt, N)
    X = 1 : N;
    Y = 1 : N;
    Z = 1 : N;
    
    X(1) = x0;
    Y(1) = y0;
    Z(1) = z0;
    
    for i = 1 : 1 : N - 1
        [X(i + 1), Y(i + 1), Z(i + 1)] = euler_method_1_step(a, b, c, X(i), Y(i), Z(i), dt);
    end
end

function [xa, ya, za] = euler_method_1_step(a, b, c, x, y, z, dt)
    xa = x + xt(a, b, c, x, y, z) * dt;
    ya = y + yt(a, b, c, x, y, z) * dt;
    za = z + zt(a, b, c, x, y, z) * dt;
end

function [X, Y, Z] = runge_kutta_method(a, b, c, x0, y0, z0, dt, N)
    X = 1 : N;
    Y = 1 : N;
    Z = 1 : N;
    
    X(1) = x0;
    Y(1) = y0;
    Z(1) = z0;
    
    for i = 1 : 1 : N - 1
        [X(i + 1), Y(i + 1), Z(i + 1)] = runge_kutta_method_1_step(a, b, c, X(i), Y(i), Z(i), dt);
    end
end

function [xr, yr, zr] = runge_kutta_method_1_step(a, b, c, x, y, z, h)
    ax = xt(a, b, c, x, y, z);
    ay = yt(a, b, c, x, y, z);
    az = zt(a, b, c, x, y, z);
    
    bx = xt(a, b, c, x + (h / 2) * ax, y + (h / 2) * ay, z + (h / 2) * az);
    by = yt(a, b, c, x + (h / 2) * ax, y + (h / 2) * ay, z + (h / 2) * az);
    bz = zt(a, b, c, x + (h / 2) * ax, y + (h / 2) * ay, z + (h / 2) * az);
    
    cx = xt(a, b, c, x + (h / 2) * bx, y + (h / 2) * by, z + (h / 2) * bz);
    cy = yt(a, b, c, x + (h / 2) * bx, y + (h / 2) * by, z + (h / 2) * bz);
    cz = zt(a, b, c, x + (h / 2) * bx, y + (h / 2) * by, z + (h / 2) * bz);
    
    dx = xt(a, b, c, x + h * cx, y + h * cy, z + h * cz);
    dy = yt(a, b, c, x + h * cx, y + h * cy, z + h * cz);
    dz = zt(a, b, c, x + h * cx, y + h * cy, z + h * cz);

    xr = x + (h / 6) * (ax + 2 * bx + 2 * cx + dx);
    yr = y + (h / 6) * (ay + 2 * by + 2 * cy + dy);
    zr = z + (h / 6) * (az + 2 * bz + 2 * cz + dz);
end

function res = xt(a, b, c, x, y, z)
    res = -y - z;
end

function res = yt(a, b, c, x, y, z)
    res = x + a * y;
end

function res = zt(a, b, c, x, y, z)
    res = b + z * (x - c);
end
