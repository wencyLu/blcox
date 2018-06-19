function [f,Ax,L]=logl(X,b,cens)
%--------------------------------------------------------------------------
% logl.m: A two term Taylor expansion of l(b)
%--------------------------------------------------------------------------
% input: 
%      X -- the dataset X(all features)
%      b -- the regression coefficients 
%      cens --  y is survival time if cens = 1, and y is right censoring time if cens = 0
%   
% output:
%      f -- intermediate variables
%      Ax  -- intermediate variables
%      L  -- A two term Taylor expansion of l(b)
% author: Wenxin Lu
% date: 1 May 2018

[n,p] = size(X);
Xb = X*b;
r = exp(Xb);
risksum = cumsum(r, 'reverse');
% risksum = risksum(atrisk);
L = cens'*(Xb - log(risksum));
% Compute first derivative dL/db
for nn=1:n
    temp1=zeros(n,1);
    temp1(1:nn,1)=1;       %if ti<tnn ,temp1=1;
    temp2=r(nn)*temp1./risksum;
    A(nn,1)=temp2'*cens; 
    temp3=r(nn)^2*temp1./(risksum.^2);
    AA(nn,1)=temp3'*cens;
end
dl =cens-A;

% Compute second derivative d2L/db2
ddl=A-AA;
z=Xb+pinv(ddl)*dl;
AAA=diag(ddl);
% L1=(z-Xb)'*AAA*(z-Xb);
m=sqrt(AAA);
f=m*z;
Ax=m*X;
L2=(f-Ax*b)'*(f-Ax*b);
end
