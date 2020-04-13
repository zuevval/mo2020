% ----  Running Python script to create `steps.txt`  ----
system('cd .. && python -m src.GradientDescent src/steps.txt'); % Windows
% status = system('cd .. ; python -m src.GradientDescent src/steps.txt'); % Linux

% ---- plotting function and points from `steps.txt` ----
f1 = @(X,Y)(X.^2 + Y.^2 - X + 2.* Y);
plotDescent(f1, 'steps.txt', -1:0.02:0.5, -1:0.02:-0.3)

%---- to visualize our real function, run `GradientDescent.py` ----
%---- and execute two lines above                              ----
% f2 = @(X, Y)(X.^2 + Y.^2 + cos(X+3.*Y) - X + 2.*Y); % our real function
% plotDescent(f2, 'real_function_steps.txt', -0.5:0.02:0.5, -1.3:0.02:0.2)
