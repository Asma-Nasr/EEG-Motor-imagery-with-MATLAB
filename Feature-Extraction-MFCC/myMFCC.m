%     melcepst( sig  ,fs ,W  ,nc,p,n,inc
nc=12;p= 16;n=128;inc=64;

c_FP1=melcepst(filFP1,250,'M',nc,p,n,inc);
c_FP2=melcepst(filFP2,250,'M',nc,p,n,inc);
c_F3=melcepst(filF3,250,'M',nc,p,n,inc);
c_F4=melcepst(filF4,250,'M',nc,p,n,inc);

c_C3=melcepst(filC3,250,'M',nc,p,n,inc);
c_C4=melcepst(filC4,250,'M',nc,p,n,inc);
c_P3=melcepst(filP3,250,'M',nc,p,n,inc);
c_P4=melcepst(filP4,250,'M',nc,p,n,inc);

c_O1=melcepst(filO1,250,'M',nc,p,n,inc);
c_O2=melcepst(filO2,250,'M',nc,p,n,inc);
c_F7=melcepst(filF7,250,'M',nc,p,n,inc);
c_F8=melcepst(filF8,250,'M',nc,p,n,inc);

c_T3=melcepst(filT3,250,'M',nc,p,n,inc);
c_T4=melcepst(filT4,250,'M',nc,p,n,inc);
c_T5=melcepst(filT5,250,'M',nc,p,n,inc);
c_T6=melcepst(filT6,250,'M',nc,p,n,inc);

c_FZ=melcepst(filFZ,250,'M',nc,p,n,inc);
c_CZ=melcepst(filCZ,250,'M',nc,p,n,inc);
c_PZ=melcepst(filPZ,250,'M',nc,p,n,inc);


FE_H=[c_FP1,c_FP2,c_F3,c_F4,c_C3,c_C4,c_P3,c_P4,c_O1,c_O2,c_F7,c_F8,c_T3,c_T4,c_T5,c_T6,c_FZ,c_CZ,c_PZ]';

