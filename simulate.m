function [data,censored_percent,b]=simulate(n,q,bb,rho,k)
%--------------------------------------------------------------------------
% simulaiton.m: generate a simulation dataset
%--------------------------------------------------------------------------
% input: 
%      n -- the numble of the samples of dataset
%      q -- the numble of the features(genes) of dataset
%      bb --  coefficients of prognostic genes, for example:  b1= 0.8,
%             b2=-0.8, b3=0.8, b4=-0.8
%      rho -- the correlation coefficient, for example: rho= 0, 0.3, 0.6
%      k  --  a parameter used to adjust the censored_percent
%   
% output:
%      data -- the simulation dataset
%      censored_percent  -- the propotion of the samples with censoring
%                           time to all samples
%      b  -- coefficients of all features(genes)
% author: Wenxin Lu
% data: 1 May 2018

%% 
b=[bb zeros(1,q-size(bb,2))];
for i=1:q+1
    r(:,i)=normrnd(0,1,n,1);
end
for i=1:q
    X(:,i)=r(:,i+1)*sqrt(1-rho)+r(:,1)*sqrt(rho);
end
X=zscore(X);  %standardization of data

%% survival time T
a=0.1;
w=1;
u=unifrnd(0,100,n,1);
T=abs((1/a)*log(1-a*log(u)./(w*exp(X*b'))));


%% censoring time
h=unifrnd(0,100,n,1)*k;
t=zeros(n,1);
censor=zeros(n,1);
for i=1:n
    t(i)=min(h(i),T(i));
    if T(i)<=h(i)
       censor(i)=1;
    else
        censor(i)=0;
    end
end

%% T is survival time if censor = 1, and T is right censoring time if censor = 0
data=[X censor T];   
[n,p]=size(data);
censored_percent=size((find(data(:,q+1)==0)),1)/n;
