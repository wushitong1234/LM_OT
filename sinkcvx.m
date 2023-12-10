function [vars,Err] = sinkcvx(P,F,Nub,Tx)
n = length(P); m = length(F);
% Generate random problem data.
% X1=P ==> Cp*X=P
Ap = eye(n);
Bp = ones(1,m);
Cp = kron(Ap,Bp);
% X^T1=F ==> Ct*X=F
At = ones(1,n);
Bt = eye(m);
Ct = kron(At,Bt);

% solve by cvx 
cvx_begin
    variables x(m*n) % x×ªÖÃºó×öreshape
    maximize sum(entr(x))
    % minimize sum(x.*log(x))
    subject to
        Cp*x == P;
        Ct*x == F;
        Nub*x >= Tx;
cvx_end
vars = x;
Err(1) = norm(Cp*x-P,1);
Err(2) = norm(Ct*x-F,1);
end
