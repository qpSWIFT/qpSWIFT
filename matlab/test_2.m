clc;
clear all;
close all;

load('casestudies.mat');

for i = 1:length(MM)

    P = MM(i).P;
    q = MM(i).q;
    A = MM(i).A;
    b = MM(i).b;
    G = MM(i).G;
    h = MM(i).h;
    
    options = optimoptions('quadprog','Display','off');
    [X,FVAL,EXITFLAG] = quadprog(P,q,G,h,A,b,[],[],[],options);
    
    
    opts.REALTOL = 1e-7;
    [sol1,basic_info1,adv_info1] = qpSWIFT(sparse(P),q,sparse(A),b,sparse(G),h,opts);
    
    [sol2,basic_info2,adv_info2] = qpSWIFT_pd_united(sparse(P),q,sparse(A),b,sparse(G),h,opts);
    
    [x,s,fval,jj,x_sol,s_sol,z_sol] = qpSWIFT_matlab(P,q,G,h,A,b);
        
        disp(['CaseStudy :',num2str(i)]);
    if EXITFLAG == 1
        disp('QUADPROG is succesful');
    else
        disp('QUADPROG is not succesful');
    end
    
    if basic_info1.Iterations ~= 100
        disp('qpSWIFT is successful');
    else
        disp('qpSWIFT is not successful');
    end
    
    if basic_info2.Iterations ~= 100
        disp('Modfified qpSWIFT is successful');
    else
        disp('Modfified qpSWIFT is not successful');
    end
    
    if jj ~= 100
        disp('Custom QP is successful');
        fprintf("\n\n");
    else
        disp('Custom QP is not successful');
        fprintf("\n\n");
    end
    
end
    
