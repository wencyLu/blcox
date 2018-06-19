function X_n=fillna(X)
% Fill missing values with mean
[N,D]=size(X);
MeanX=[];
for i=1:D
   tempx=[];
    for j=1:N
        if isnan(X(j,i))
            a=0;
        else
            tempx=[tempx;X(j,i)];
        end
    end
    MeanX=[MeanX mean(tempx)];
end
for i=1:D
    for j=1:N
        if isnan(X(j,i))
            X(j,i)=MeanX(i);
        end
    end
end
X_n = X;
