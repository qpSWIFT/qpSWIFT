clc;
clearvars;
close all;

H_QP=1;
f_QP=-0.231;
G_QP=1;
h_QP=0;
opts.VERBOSE = 0;
% opts.MAXITER=5;[Xpd,FVALpd] = qpSWIFT_pd_united(sparse(H_QP),[f_QP],sparse([(G_QP)]),h_QP,opts);
[X,FVAL] = qpSWIFT(sparse(H_QP),[f_QP],sparse([(G_QP)]),h_QP,opts);
[x,s,fval,iter] = qpSWIFT_matlab(H_QP,f_QP,G_QP,h_QP,[],[]);