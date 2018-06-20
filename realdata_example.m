clear;clc
%load ..\blcox\data\DLBCL-R2002-240-x-7399.mat
load ..\blcox\data\AML-Bullinger2004-116-x-6283.mat
%load ..\blcox\data\Veer2002-78-x-4751.mat
warning off
[N,D]=size(X);
train_size=floor(N*0.7);

%%
t=20;  % repeat each experiment t times
c_sum_glmcox=0; %the sum of c_index of glmcox in t runs
c_sum_blcox=0; %the sum of c_index of blcox in t runs
c_n_glmcox=0; 
c_n_blcox=0;
t_sum_glmcox=0; %the sum of runtime of glmcox in t runs
t_sum_blcox=0;  %the sum of runtime of blcox in t runs
f_select_glmcox=0; %the sum of features selected of glmcox in t runs
f_select_blcox=0; %the sum of features selected of blcox in t runs

for jj=1:t
index=randperm(N);
train_index=index(1:train_size);
test_index=index(train_size+1:N);

Censored=delta;
X=fillna(X);

X_train=X(train_index,:);
X_train=zscore(X_train);
Censored_train=Censored(train_index,:);
y_train=Y(train_index,:);

X_test=X(test_index,:);
X_test=zscore(X_test);
Censored_test=Censored(test_index,:);
y_test=Y(test_index,:);

y1=[y_train Censored_train];
for j=1:train_size
    if y1(j,1)==0
       y1(j,1)=0.01;
    end
end

%% glmcox, which can be downloaded in http://web.stanford.edu/~hastie/glmnet_matlab/
t1=clock;
cvfit = cvglmnet(X_train, y1, 'cox');
b_glmcox=cvglmnetCoef(cvfit,'lambda_min');
c_glmcox=cindex(X_test,y_test,Censored_test,b_glmcox);
if c_glmcox>=0.5
   c_sum_glmcox=c_sum_glmcox+c_glmcox;
   c_n_glmcox=c_n_glmcox+1;
end

%% blcox
t2=clock;
lambda=0.2;
b_blcox=blcox(X_train,y_train,Censored_train,lambda);
c_blcox = cindex(X_test,y_test,Censored_test,b_blcox);
if c_blcox>=0.5
   c_sum_blcox=c_sum_blcox+c_blcox;
   c_n_blcox=c_n_blcox+1;
end
t3=clock;
t_sum_glmcox=t_sum_glmcox+etime(t2,t1);
t_sum_blcox=t_sum_blcox+etime(t3,t2);

%% features selected &&  real features selected
fs_g=0;
fs_b=0;
for i=1:D
    if b_glmcox(i)~=0      
        fs_g=fs_g+1;
    end
    if b_blcox(i)~=0     
        fs_b=fs_b+1;
    end
end

f_select_glmcox=f_select_glmcox+fs_g; 
f_select_blcox=f_select_blcox+fs_b; 
end
%% rescult
c_mean_glmcox=c_sum_glmcox/c_n_glmcox;
c_mean_blcox=c_sum_blcox/c_n_blcox;

t_mean_glmcox=t_sum_glmcox/t;
t_mean_blcox=t_sum_blcox/t;

f_mean_g = f_select_glmcox/c_n_glmcox;
f_mean_b = f_select_blcox/c_n_blcox;
