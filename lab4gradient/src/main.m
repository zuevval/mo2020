% ----  Running Python script to create `steps.txt`  ----
system('cd .. && python -m src.GradientDescent src/steps.txt'); % Windows
% status = system('cd .. ; python -m src.GradientDescent src/steps.txt'); % Linux

% ---- plotting function and points from `steps.txt` ----
f1 = @(X,Y)(X.^2 + Y.^2 - X + 2.* Y);
plotDescent(f1, 'steps.txt', -1:0.02:0.5, -1:0.02:-0.3)

% f2 = @(X, Y)(X.^2 + Y.^2 + cos(X+3.*Y) - X + 2.*Y); % our real function
