function [filFP1,filFP2,filF7,filF3,filFZ,filF4,filF8,filT5,filC3,filCZ,filC4,filT6,filT3,filP3,filPZ,filP4,filT4,filO1,filO2]=filterX(FP1,FP2,F7,F3,FZ,F4,F8,T5,C3,CZ,C4,T6,T3,P3,PZ,P4,T4,O1,O2)
% function [c3,filc4]=filterX(c3,c4)
% TO ELEIMINATE THE FIRST 2 SECOONDS
%number of instants will be deleted
%------------------
%FP1 FP2 F3 F4 C3 C4 P3 P4 O1 O2 F7 F8 T3 T4 T5 T6 FZ CZ PZ
nFP1= normc(FP1) ;
nFP2= normc(FP2) ;

nF7= normc(F7) ;
nF3= normc(F3) ;
nFZ= normc(FZ) ;
nF4= normc(F4) ;
nF8= normc(F8) ;

nT5= normc(T5) ;
nC3= normc(C3) ;
nCZ= normc(CZ) ;
nC4= normc(C4) ;
nT6= normc(T6) ;

nT3= normc(T3) ;
nP3= normc(P3) ;
nPZ= normc(PZ) ;
nP4= normc(P4) ;
nT4= normc(T4) ;

nO1= normc(O1) ;
nO2= normc(O2) ;

%TO FILTER
[N, Wp] = ellipord([8 30]/250,[7 34]/250,.1,40);
%[N, Wp] = ellipord([8/250 35/250],[7/250 37/250],0.1,60);
% theta- beta, delta-gamma
[B,A] = ellip(N,.1,40,Wp);
filFP1=filter(B,A,nFP1);
filFP2=filter(B,A,nFP2);

filF7=filter(B,A,nF7);
filF3=filter(B,A,nF3);
filFZ=filter(B,A,nFZ);
filF4=filter(B,A,nF4);
filF8=filter(B,A,nF8);

filT5=filter(B,A,nT5);
filC3=filter(B,A,nC3);
filCZ=filter(B,A,nCZ);
filC4=filter(B,A,nC4);
filT6=filter(B,A,nT6);

filT3=filter(B,A,nT3);
filP3=filter(B,A,nP3);
filPZ=filter(B,A,nPZ);
filP4=filter(B,A,nP4);
filT4=filter(B,A,nT4);

filO1=filter(B,A,nO1);
filO2=filter(B,A,nO2);