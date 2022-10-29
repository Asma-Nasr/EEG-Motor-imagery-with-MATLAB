%%
data=load('ALLF.csv');
Target=load('at.csv');
X = data;
y = Target;

%% 80 :20
rand_num = randperm(2222);
X_train = X(rand_num(1:1777),:);
y_train = y(rand_num(1:1777),:);

X_test = X(rand_num(1777:end),:);
y_test = y(rand_num(1777:end),:);
%% CV partition
c= cvpartition(y_train,'k',5);
%% feature selection
opts=statset('display','iter');
fun=@(train_data,train_labels,test_data,test_labels)...
     sum(predict(fitcsvm(train_data,train_labels,'KernelFunction','linear'),test_data) ~= test_labels)
[fs,history]=sequentialfs(fun, X_train,y_train,'cv',c,'options',opts,'nfeatures',25)
%% Best heperparameters
X_train_w_best_features = X_train(:,fs);
Md1 = fitcsvm(X_train_w_best_features,y_train,'KernelFunction','linear',...
    'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus','ShowPlots',true));
%% test
x_test_w_best_features = X_test(:,fs);
accuracy= sum(predict(Md1,x_test_w_best_features)==y_test)/length(y_test)*100
cmpfe

%% accuracy
TS5=TS(:,fs);
B=predict(Md1,TS5);
acccuracyi=sum(B==ti)/length(ti)*100
T3=Data3(:,fs);
%T3=T3(:,fs);
w=ones(1111,1);
B1=predict(Md1,T3);
acccuracy_Ri=sum(B1==w)/length(w)*100
TSw=Data4(:,fs);
%TSw=TSw(:,fs);
B2=predict(Md1,TSw);
ww=-1*ones(1111,1);
acccuracy_Fi=sum(B2==ww)/length(ww)*100

%%

%Results = confusionmat(y_test,predict(Md1,x_test_w_best_features));
%save('Md1','Md1')

