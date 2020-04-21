y = [1;1];
x1 = -1.5:0.02:1;
x2 = -1:0.02:1;
[X, Y] = meshgrid(x1, x2);
Z = zeros(size(X));
for i = 1:size(x1,2)
    for j = 1:size(x2,2)
        Z(j,i) = y'*my_hesse([x1(i);x2(j)])*y;
    end
end
figure
mesh(X,Y,Z)