function b=solve_byccd(z,X,lambda)
%--------------------------------------------------------------------------
% solve_byccd.m: the coordinate descent method to update beta
%--------------------------------------------------------------------------
% input: 
%      z -- intermediate variables
%      X -- the dataset X(all features)
%      lambda --  a parameter 
% output:
%      b -- new beta
% author: Wenxin Lu
% date: 1 May 2018
%%
[N,M]=size(X);
b=zeros(M,1);
%beta=zeros(M,1);
ddl=ones(N,1);

X_2=X.^2;
tempw=repmat(ddl,1,M);
norm_w=sum(tempw.*X_2)';

%%
DDL=repmat(ddl,1,M);
Z=repmat(z,1,M);
yuzhi=lambda;%*mean(ddl);

%%iterate until converge
for k=1:30
%     %calculate beta
%     for m=1:M
%         temp1=X*b-b(m)*X(:,m);
%         temp2=ddl.*X(:,m).*(z-temp1);
%         beta(m)=sum(temp2);
%     end

    B=repmat(b',N,1);
    temp1=repmat(X*b,1,M)-B.*X;
    temp2=DDL.*X.*(Z-temp1);
    beta=sum(temp2)';
    
    index1=find(abs(beta)<=yuzhi);
    index2=find(abs(beta)>yuzhi);
    b_new=zeros(M,1);
    b_new(index1,1)=0;
    b_new(index2,1)=sign(beta(index2)).*(abs(beta(index2))-yuzhi)./norm_w(index2);
    
    E=abs(b-b_new);
    [~,max_index]=max(E);
    b(max_index)=b_new(max_index);
    ftemp=ddl.*(z-X*b).^2;
    f_i=sum(ftemp)+lambda*norm(b,1);
    f(k)=f_i;
end
