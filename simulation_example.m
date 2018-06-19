clear;clc;
warning off
addpath('..\blcox\glmnet_matlab')

%%  generate a simulation dataset including 400 samples with 1000 genes
n=400; 
q=1000;
rho=0.3; % or rho=0 0.3 0.6
bb=[-0.5 0 0 0.8 0 0 -0.8 0 0 0.5];
k=0.4; % Adjust k, and make censored_percent approximately equal to 30%
[data,censored_percent,b]=simulate(n,q,bb,rho,k);

%%
t=10;  % repeat each experiment t times
c_sum_glmcox=0; %the sum of c_index of glmcox in 50 runs
c_sum_blcox=0; %the sum of c_index of blcox in 50 runs
c_n_glmcox=0; 
c_n_blcox=0;
t_sum_glmcox=0; %the sum of runtime of glmcox in 50 runs
t_sum_blcox=0;  %the sum of runtime of blcox in 50 runs
mse_sum_glmcox=0; %the sum of MSE of glmcox in 50 runs
mse_sum_blcox=0;  %the sum of MSE of blcox in 50 runs

f_select_glmcox=0; %the sum of features selected of glmcox in 50 runs
f_select_blcox=0; %the sum of features selected of blcox in 50 runs
real_f_glmcox=0;  %the sum of real features selected of glmcox in 50 runs
real_f_blcox=0;   %the sum of real features selected of blcox in 50 runs

%%% repeat each experiment t times
for jj=1:t
nn=300;   %% or n=150
train_index=randperm(400);
test_index=train_index(nn:nn+50);
train_index(nn+1:end)=[];

train_data = data(train_index,:);
X_train = fillna(train_data(:,1:q)); 
X_train = zscore(X_train); %%%%%归一化
Censored_train=train_data(:,q+1);
y_train=train_data(:,q+2);

test_data = data(test_index,:);
X_test = fillna(test_data(:,1:q));
X_test=zscore(X_test);
Censored_test=test_data(:,q+1);
y_test=test_data(:,q+2);

%% glmcox, which can be downloaded in http://web.stanford.edu/~hastie/glmnet_matlab/
t1=clock;
yy=[y_train Censored_train];
cvfit = cvglmnet(X_train, yy, 'cox');
b_glmcox=cvglmnetCoef(cvfit,'lambda_min');
c_glmcox=cindex(X_test,y_test,Censored_test,b_glmcox);
if c_glmcox>=0.5
   c_sum_glmcox=c_sum_glmcox+c_glmcox;
   c_n_glmcox=c_n_glmcox+1;
end

%% blcox
t2=clock;
lambda=0.1;
b_blcox=blcox(X_train,y_train,Censored_train,lambda);
c_blcox=cindex(X_test,y_test,Censored_test,b_blcox);
if c_blcox>=0.5
   c_sum_blcox=c_sum_blcox+c_blcox;
   c_n_blcox=c_n_blcox+1;
end
t3=clock;

t_sum_glmcox=t_sum_glmcox+etime(t2,t1);
t_sum_blcox=t_sum_blcox+etime(t3,t2);

mse_glmcox=sum((b_glmcox-b').^2)/q;
mse_blcox=sum((b_blcox-b').^2)/q;
mse_sum_glmcox=mse_sum_glmcox+mse_glmcox;
mse_sum_blcox=mse_sum_blcox+mse_blcox;

%% features selected &&  real features selected
for i=1:q
    if b(i)~=0
        b(i)=1;
    end
    if b_glmcox(i)~=0      
        b_glmcox(i)=1;
    end
    if b_blcox(i)~=0     
        b_blcox(i)=1;
    end
end

r_b=0;
fs_g=0;
rf_g=0;
fs_b=0;
rf_b=0;
for i=1:q
    if b(i)==1
        r_b=r_b+1;
     end
    if b_glmcox(i)==1
       fs_g=fs_g+1;
       if b_glmcox(i)==b(i)
          rf_g=rf_g+1;
       end

    end
    if b_blcox(i)==1
       fs_b=fs_b+1;
       if b_blcox(i)==b(i)
          rf_b=rf_b+1;
       end
    end
end
f_select_glmcox=f_select_glmcox+fs_g; 
f_select_blcox=f_select_blcox+fs_b; 
real_f_glmcox=real_f_glmcox+rf_g;  
real_f_blcox=real_f_blcox+rf_b;
end

%% rescult
c_mean_glmcox=c_sum_glmcox/c_n_glmcox;
c_mean_blcox=c_sum_blcox/c_n_blcox;

t_mean_glmcox=t_sum_glmcox/t;
t_mean_blcox=t_sum_blcox/t;

mse_mean_glmcox=mse_sum_glmcox/t;
mse_mean_blcox=mse_sum_blcox/t;

f_mean_g = f_select_glmcox/t;
f_mean_b = f_select_blcox/t;
r_mean_g = real_f_glmcox/t;  
r_mean_b = real_f_blcox/t; 
