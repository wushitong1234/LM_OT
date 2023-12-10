run('D:/project/queue_SDP/queue_solver/cvx/cvx_startup.m')
clc;clear
% set conition
n = 16; % number n^2 
m = 40; % number m^2
yint = -8; % m*ydelta=2*|yint|
ydelta= 2*(-yint)/(m-1); % interval length
y2delta = ydelta/2; % half interval
sigma = 1; % parameter in normal distribution
% set X and Y
X = zeros(n^2,2);
for i = 1:n
    for j = 1:n
        X((i-1)*n+j,1) = (n-1)-(i-1)*2;
        X((i-1)*n+j,2) = (n-1)-(j-1)*2;
    end
end
xdelta = sum(sum(X.*X))/(n*n);
X = X./sqrt(xdelta);
% set X and Y
Y = zeros(2,m^2);
for x = 1:m
    for j = 1:m
        Y(1,m*(x-1)+j) = yint+(x-1)*ydelta;
        Y(2,m*(x-1)+j) = yint+(j-1)*ydelta;
    end
end

Nu = X*Y; % 变量x转置后做reshape
%save('Nu_data.txt', 'Nu','-ascii');
Nub = reshape(Nu',[1,(m*n)^2]); % note: reshape 按照列排，所以要加转置 
%save('Nub_data.txt', 'Nub','-ascii');

% the marginal distribution
P = ones(n^2,1)/n^2;
%save('P_data.txt', 'P','-ascii');

% the product X and Y in each item 
et = 0.9; th = pi/18;
H1 = [1 0; 0 et]; H2 = [cos(th) -sin(th); sin(th) cos(th)]; H = H1*H2;
Nx = X*H*X'; Tix = diag(Nx);
Tx = sum(P.*Tix);
%save('Tx_data.txt', 'Tx','-ascii');

% the marginal distribution
F = zeros(m^2,1);
for x = 1:m^2
    for j = 1:n^2
        ss = normcdf(Y(1,x)+y2delta,X(j,1),sigma)-normcdf(Y(1,x)-y2delta,X(j,1),sigma);
        sy = normcdf(Y(2,x)+y2delta,X(j,2),sigma)-normcdf(Y(2,x)-y2delta,X(j,2),sigma);
        F(x) = F(x) + ss*sy*P(j); % 之前代码写错了F(i) = F(i) + ss*sy;
    end
end
F = F/sum(F);
%save('F_data.txt', 'F','-ascii');



% % run cvx and sinkhorn
% % the cvx time 
% for i = 1:10
%     ans_time_cvx = tic;
%     [vars,Err_cvx] = sinkcvx(P,F,Nub,Tx);
%     ans_time_cvx = toc(ans_time_cvx);
%     time_cvx(i) = ans_time_cvx;
%     ans_cvx(i) = -sum(entr(vars))-2*log(ydelta);
% end

% the sinkhorn time
for i = 1:100
    ans_time_sink = tic;
    [gamma,Err_sink] = sinkhorn(P,F,Nu,Tx);
    ans_time_sink = toc(ans_time_sink);
    time_sink(i) = ans_time_sink;
    ans_sink(i) = sum(sum(log(gamma).*gamma))-2*log(ydelta);
end

% save('new_LMs_4_10.mat','time_sink','ans_sink','time_cvx','ans_cvx')

save('sink_LMs_16_40.mat','time_sink','ans_sink')

% save('cvx_LMs_4_10.mat','time_cvx','ans_cvx')

