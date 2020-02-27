sv2 = [0 49/55 1/11 0]
Nk = [1, 2]
Nk_plus = Nk
AMN = [5 5 6 4; 4 5 -5 -4]
AMNk = AMN(:, 2:3)
BNkM = inv(AMNk)
cN = [1 2 3 4]'
cNk = [2 3]'
ykM = BNkM'*cNk
dkN = cN - (AMN)'*ykM

ukNk = BNkM*AMN(:,1)