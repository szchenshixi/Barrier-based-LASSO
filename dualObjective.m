function G = dualObjective(v)
global g_y;   % nx1
G = -(1/4)*(v')*v - v'*g_y;
end