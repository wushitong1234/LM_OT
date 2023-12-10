function [gamma,Err] = sinkhorn(P,F,Nu,Tx)
n = length(P); m = length(F);
% perform the sinkhorn 
niter_sinkhorn = 1000;
tol_newt = 1e-18;
tol_sink = 1e-16;

% function G and G'
Ga = @(ax,bx,lax) ax'*(Nu.*exp(lax*Nu))*bx-Tx;
Gad = @(ax,bx,lax) ax'*((Nu.*Nu).*exp(lax*Nu))*bx;
K = @(lax) exp(lax*Nu); 

% set up for sinkhorn
a = ones(n,1); b = ones(m,1); 
la = 1; Err = [];

% the main iteration
for i = 1:niter_sinkhorn
    b = F ./ ((K(la))'*a);
    Err(i,1) = norm(a.*(K(la)*b)-P, 1); 
    a = P ./ (K(la)*b);
    Err(i,2) = norm(b.*((K(la))'*a)-F, 1);
    Err(i,3) = norm(la*Ga(a,b,la)/(la+1));
    if Ga(a,b,0) > 0 % only for AWGN!!
        la = 0;
    else
        la = newton(a,b,la,tol_newt,Ga,Gad);
    end
%     if max(Err(i,:))<tol_sink
%         break;
%     end
end
gamma = diag(a)*K(la)*diag(b);
end

% Newton method 
function [lax] = newton(ax,bx,lax,tol,Ga,Gad)
for i = 1:20
    lax = lax - Ga(ax,bx,lax)/Gad(ax,bx,lax);
    if norm(Gad(ax,bx,lax)) < tol
        break;
    end
end
end
