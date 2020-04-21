function H = my_hesse(x) % x is a column vector
x1 = [1,3]*x;
H = [2 - cos(x1), -3*cos(x1); -3*cos(x1), 2 - 9*cos(x1)];
end