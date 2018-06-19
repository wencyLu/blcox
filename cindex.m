function c=cindex(X,y,Censored,b)
%--------------------------------------------------------------------------
% c_index.m: Calculate CI value
%--------------------------------------------------------------------------
% input: 
%      X -- the dataset X(all features)
%      y -- survival time
%      Censored --  y is survival time if cens = 1, and y is right censoring time if cens = 0
%      b -- beta
% output:
%      c -- the CI value
% author: Wenxin Lu
% date: 1 May 2018

[N,~]=size(X);
c1=0;
c2=0;
c_sum1 = 0;
c_sum2 = 0;
for i=1:N-1
    for j=i+1:N
        if (y(i)<y(j)&&Censored(i)==1)
            c_sum1=c_sum1+1;
        end
        if (y(j)<y(i)&&Censored(j)==1)
            c_sum2=c_sum2+1;
        end
        if (y(i)<y(j)&&(X(i,:)*b)>(X(j,:)*b)&&Censored(i)==1)
            c1=c1+1;   
        end
        if (y(j)<y(i)&&(X(j,:)*b)>(X(i,:)*b)&&Censored(j)==1)
            c2=c2+1;   
        end
    end
end
c=(c1+c2)/(c_sum1+c_sum2);
