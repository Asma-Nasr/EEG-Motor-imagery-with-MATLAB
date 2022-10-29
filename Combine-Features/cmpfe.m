Data1=load('MRF.csv')';
Data2=load('MLF.csv')';
Data3=load('IRF.csv')';
Data4=load('ILF.csv')';

data=[Data1;Data2];

csvwrite('AllF.csv',data)

%%
%a=[ones(456,1) ones(684,1);ones(684,1) -1*ones(684,1);-1*ones(684,1) ones(684,1);-1*ones(684,1) -1*ones(684,1)];
%t=[-10*ones(684,1);-5*ones(684,1);5*ones(684,1);10*ones(684,1)];
%t=[-2*ones(456,1);-1*ones(456,1);1*ones(456,1);2*ones(456,1)];
rr=[1*ones(1111,1);-1*ones(1111,1)];
csvwrite('at.csv',rr)
%%
%TEST imagery
%TS=[Data3;Data7;Data9;Data12];
TS=[Data3;Data4];
csvwrite('TI.csv',TS)
ti=[ones(1111,1);-1*ones(1111,1)];
%ti=[ones(228,1) ones(228,1);ones(228,1) -1*ones(228,1);-1*ones(228,1) ones(228,1);-1*ones(228,1) -1*ones(228,1)];

%csvwrite('TL2.csv',TS1)

