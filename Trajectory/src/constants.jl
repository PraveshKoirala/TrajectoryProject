FR  = 1   # flag radius
CR = 1    # collision radius
K = 5   # number of knot points x0, x1, x2, x3 back to x0 (not counted)

flags = [    0 4;
             -2 2;
             5 -3;
             3 1;
             8 2;
             2 2;]
flags[:, 1] .-= 3
start2 = [1, -2]
start1 = [-1, -2]