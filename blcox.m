function [b,lambda] = blcox(X,y,cens,lambda)
%--------------------------------------------------------------------------
% blcox.m: BLasso-Cox model
%--------------------------------------------------------------------------
% input: 
%      X -- the dataset X(all features)
%      cens --  y is survival time if cens = 1, and y is right censoring time if cens = 0
%      lambda -- a parameter
% output:
%      b -- beta
% author: Wenxin Lu
% date: 1 May 2018

[~,p]=size(X);
[~,idx] = sort(y); 
X = X(idx,:);
b=zeros(p,1);
a=0.1;
bx=0.1;
for ii=1:15
%% updata b
[z,Ax,~]=logl(X,b,cens);
b=solve_byccd(z,Ax,lambda);

%% updata lambada
lambda=(length(find(b~=0))+a)/(norm(b,1)+bx);  %updata lambda

end
end
