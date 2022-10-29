data=load('ALLF.csv');
Target=load('at.csv');
X = data;
y = Target;

%% 80 :20
%456
rand_num = randperm(2222);
X_train = X(rand_num(1:1777),:);
y_train = y(rand_num(1:1777),:);

X_test = X(rand_num(1777:end),:);
y_test = y(rand_num(1777:end),:);

%%

k =1;
%w = [0.3; 0.3; 0.2; 0.2];
%'euclidean'
%,'Standardize',1
% chebychev     good
% mahalanobis good accuracy
% correlation 79%
% spearman 72%
%clf=fitcknn(X_train,y_train,'NumNeighbors',1,'distance','euclidean');
%%
%% CV partition
c = cvpartition(y_train,'k',5);
%% feature selection
opts=statset('display','iter');
fun=@(train_data,train_labels,test_data,test_labels)...
     sum(predict(fitcknn(train_data,train_labels),test_data) ~= test_labels)
[fs,history]=sequentialfs(fun, X_train,y_train,'cv',c,'options',opts,'nfeatures',30)
%% Best heperparameters
X_train_w_best_features = X_train(:,fs);
%%
clf= fitcknn(X_train_w_best_features,y_train,'IncludeTies',false,'Distance','correlation',...
    'NumNeighbors',k);
%clf = knnsearch(X_train,y_train,'k',10,'IncludeTies',true,'Distance','cityblock');
rng(1); % For reproducibility
Cl = crossval(clf);
classError = kfoldLoss(Cl)
x_test_w_best_features = X_test(:,fs);
A=predict(clf,x_test_w_best_features)
%Imagery
cmpfe
TS=TS(:,fs);
B=predict(clf,TS);
%% accuracy
accuracy= sum(A==y_test)/length(y_test)*100
acccuracyi=sum(B==ti)/length(ti)*100
T3=Data3(:,fs);
%T3=T3(:,fs);
w=ones(1111,1);
B=predict(clf,T3);
acccuracy_Ri=sum(B==w)/length(w)*100
TSw=Data4(:,fs);
%TSw=TSw(:,fs);
B2=predict(clf,TSw);
ww=-1*ones(1111,1);
acccuracy_Fi=sum(B2==ww)/length(ww)*100





