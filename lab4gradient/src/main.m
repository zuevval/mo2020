% ---- Tune parameters here ---
x0 = 0.;
y0 = 0.;
alg = 'steepest'; % alg = 'steepest' or alg = 'newton'

% ----  Running Python script to create `steps.txt`  ----
cmd = 'cd .. && python -m src.grad_cli'; % Windows
% cmd = 'cd .. ; python -m src.grad_cli'; % Linux
full_cmd = [cmd sprintf(' --x0 %f ', x0) sprintf(' --y0 %f ', y0)...
    ' --out src/steps.txt --alg ' alg];
system(full_cmd); 

% ---- plotting function and points from `steps.txt` ----
f1 = @(X,Y)(X.^2 + Y.^2 - X + 2.* Y);
plotDescent(f1, 'steps.txt')

%---- to visualize our real function, run `GradientDescent.py` ----
%---- and execute two lines below                              ----
% f2 = @(X, Y)(X.^2 + Y.^2 + cos(X+3.*Y) - X + 2.*Y); % our real function
% plotDescent(f2, 'real_function_steps.txt', -0.5:0.02:0.5, -1.3:0.02:0.2)
