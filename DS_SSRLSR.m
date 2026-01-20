function [J,Z] =DS_SSRLSR( X , Par )
%X的横为N个样本数，纵为D维度
K=X'*X;
mu = 1e-4;
muMax=1e6;
iterMax=Par.maxIter;
[~,N] = size (X);
C = zeros (N, N);
J  = C;
Delta = zeros (N, N); % Y
Z=zeros(N,N);
%%
tol   = 1e-6;
iter    = 0;
%     err1(1) = inf; err2(1) = inf; err3(1) = inf; err4(1)=inf;
while  ( iter<iterMax )
%     Cpre = C;
%     Jpre = J;
%     Zpre = Z;
%    
 %% update C the coefficient matrix
    C=inv(2*Par.lambda_1*(Z'*Z)+mu*eye(N))*(Delta+mu*J+2*Par.lambda_1*Z'*Z);
     
    %% update J the data term matrix
    Q = (mu*C - Delta)/(Par.s*(2*Par.lambda_2+mu));
    J = SimplexProj(Q');
    J = Par.s*J';
    
    %%更新Z
    Ta=K+Par.lambda_1*eye(N);
    Tb=Par.lambda_1*(C*C'-C-C');
    Tc=K;
    Z=sylvester(Ta,Tb,Tc); 
    %% update Deltas the lagrange multiplier matrix
    Delta = Delta +mu * (J-C);
    mu=min(mu*Par.rho,muMax);
 
    %% computing errors
    iter = iter + 1;
%       if (  (err1(iter+1) <= tol && err1(iter+1)<=tol && err2(iter+1)<=tol && err3(iter+1)<=tol&& err4(iter+1)<=tol) || iter >= Par.maxIter )
%              terminate = true;
%             fprintf('  iter=%d\n',iter);
%       end
end

end
