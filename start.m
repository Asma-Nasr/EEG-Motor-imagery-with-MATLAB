% import the signals
importsignal;
%
[r1,c1]=size(P4);
[r2,c2]=size(O2);
[r3,c3]=size(F8);
[r4,c4]=size(T6);
[r5,c5]=size(C4);
[r6,c6]=size(CZ);
[r7,c7]=size(FZ);
[r8,c8]=size(F4);
[r9,c9]=size(FP2);
[r10,c10]=size(F7);
[r11,c11]=size(FP1);
[r12,c12]=size(T3);
[r13,c13]=size(F3);
[r14,c14]=size(C3);
[r15,c15]=size(T4);
[r16,c16]=size(T5);
[r17,c17]=size(P3);
[r18,c18]=size(O1);
[r19,c19]=size(PZ);
t=[ r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11 r12 r13 r14 r15 r16 r17 r18 r19];
%t=[c1 c2 c3 c4 c5 c6 c7 c8];
rr=min(t);

FP1=FP1(500:rr);
P4=P4(500:rr);
O2=O2(500:rr);
T3=T3(500:rr);
T4=T4(500:rr);
C4=C4(500:rr);
CZ=CZ(500:rr);
FZ=FZ(500:rr);
F4=F4(500:rr);
FP2=FP2(500:rr);
T5=T5(500:rr);
T6=T6(500:rr);
F3=F3(500:rr);
C3=C3(500:rr);
F7=F7(500:rr);
F8=F8(500:rr);
P3=P3(500:rr);
O1=O1(500:rr);
PZ=PZ(500:rr);
% Filter the signal

[filFP1,filFP2,filF3,filF4,filC3,filC4,filP3,filP4,filO1,filO2,filF7,filF8,filT3,filT4,filT5,filT6,filFZ,filCZ,filPZ]=filterX(FP1,FP2,F3,F4,C3,C4,P3,P4,O1,O2,F7,F8,T3,T4,T5,T6,FZ,CZ,PZ);

%---- stpectrogarm.......... channels

%----KDE----
%FeatureExtraction_PSD
myMFCC
%FeatureExtraction_PSD
csvwrite('ILF.csv',FE_H)
%--------- set targets
%[LBK,row,col] = setTargets(FE_H,10);
%H=extALL_H;
%save('LBK','LBK');
%----------------------CLASSIFICATION---------------------------------
