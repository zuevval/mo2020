x2 = -1.3:0.02:-0.7;
x1 = 0.2:0.02:0.6;
y = [0;1];
[X, Y] = meshgrid(x1, x2);
Z = zeros(size(X));
for i = 1:size(x1,2)
    for j = 1:size(x2,2)
        Z(j,i) = y'*my_hesse([x1(i);x2(j)])*y;
    end
end
figure
mesh(X,Y,Z)

min_z = Inf;
max_z = -Inf;
for fi = 0:0.05:2*pi
    y = [cos(fi);sin(fi)];
    for i = 1:size(x1,2)
        for j = 1:size(x2,2)
            z = y'*my_hesse([x1(i);x2(j)])*y;
            if z > max_z
                max_z = z
            end
            if z < min_z
                min_z = z
            end
        end
    end
end
