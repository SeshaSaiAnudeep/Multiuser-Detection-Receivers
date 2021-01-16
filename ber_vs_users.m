clc;clear variables;
N = 8;P_delta = 3;P=1000;Nsim=1000;
for kk=2:N
for ii=1:Nsim
    S = (2*randi([0,1],N,kk)-1);
    S = S/sqrt((S(:,1)'*S(:,1)));
    R = S'*S;
    while rank(R)<kk
        S = (2*randi([0,1],N,kk)-1);
        S = S/sqrt((S(:,1)'*S(:,1)));%S=S/sqrt(N)
        R = S'*S;
    end
    B = 2*randi([0,1],kk,P)-1;
    SNRdB = 6;
    pkk_dB = (0:kk-1)*P_delta +SNRdB;
    pkk_lin = 10.^(pkk_dB/10);
    A = sqrt(pkk_lin);
    aa = A/sqrt(var(A));
    aa_sd = sqrt(var(A)) -0.1*sqrt(var(A));
    a = aa*aa_sd;
%     A_inv = 1./(pkk_lin);
    A_inv = 1./(a.^2);
    A_inv = diag(A_inv);
    A = diag(A);
    n = randn(N,P)+1j*randn(N,P);
    n = n/sqrt(2);
    y = S*A*B+n;
    z = S'*y;
    b_est_SUMF = sign(real(z));
    ber_SUMF(kk,ii) = (1/(kk*P))*sum(sum(b_est_SUMF~=B));
    b_est_Deco = sign(real(R\z));
    ber_Deco(kk,ii) = (1/(kk*P))*sum(sum(b_est_Deco~=B));
    b_est_MMSE = sign(real((R+A_inv)\z));
    ber_MMSE(kk,ii) = (1/(kk*P))*sum(sum(b_est_MMSE~=B));
    b_est_SIC(kk,:) = sign(real(z(kk,:)));
    for mm = 1:kk-1
        r = R(:,mm);
        b_est_SIC(mm,:)= sign(real(z(mm,:)')-B(mm+1:end,:)'*A([mm+1:end], [mm+1:end])*r(mm+1:end));
    end
    ber_SIC(kk,ii) = (1/(kk*P))*sum(sum(b_est_SIC~=B));
    R_hat = R-eye(kk,kk);
    b = b_est_SUMF;
    for i=1:25
        b = sign(real(z-R_hat*A*b));
    end
    ber_PIC(kk,ii) = (1/(kk*P))*sum(sum(b~=B));
      % DF-IC
    FF = chol(R);
    Ft = FF';
    z_hat = (FF')\z;
    b_est_DC(kk,:) = sign(real(z_hat(kk,:)));
    for mm = 1:kk-1
        r = Ft(:,mm);
        b_est_DC(mm,:)= sign(real(z_hat(mm,:)')-B(mm+1:end,:)'*A([mm+1:end], [mm+1:end])*r(mm+1:end));
    end
    ber_DC(kk,ii) = (1/(kk*P))*sum(sum(b_est_DC~=B));
    ber_theo(kk,ii) = (1/kk)*sum(qfunc(sqrt(2)*(diag(A))));
end
end
ber_SUMF = sum(ber_SUMF,2)/Nsim;
ber_Deco = sum(ber_Deco,2)/Nsim;
ber_MMSE = sum(ber_MMSE,2)/Nsim;
ber_SIC  = sum(ber_SIC,2)/Nsim;
ber_PIC  = sum(ber_PIC,2)/Nsim;
ber_DC   = sum(ber_DC,2)/Nsim;
ber_theo = sum(ber_theo,2)/Nsim;
semilogy((2:N),ber_SUMF(2:N),'b-o','MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor',[0.5 .5 1],'linewidth',1.5);hold on;
semilogy((2:N),ber_Deco(2:N),'r->','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor',[1 .5 0.5],'linewidth',1.5);hold on;
semilogy((2:N),ber_MMSE(2:N),'-<','color',[0 0.4 0],'MarkerSize',5,'MarkerEdgeColor',[0 0.5 0],'MarkerFaceColor',[0.6 1 0.6],'linewidth',1.5);
semilogy((2:N),ber_PIC(2:N),'-^','color',[0 0 0],'MarkerSize',5,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[68,69,69]*(1/255),'linewidth',1.5);
semilogy((2:N),ber_SIC(2:N),'m-v','MarkerSize',5,'MarkerEdgeColor',[1,8/255,127/255],'MarkerFaceColor',[255,113,181]*(1/255),'linewidth',1.5);
semilogy((2:N),ber_DC(2:N),'-d','color',[84 22 120]/255,'MarkerSize',5,'MarkerEdgeColor',[112,39,195]/255,'MarkerFaceColor',[150,0,205]*(1/255),'linewidth',1.5);
semilogy((2:N),ber_theo(2:N),'-p','color',[0 204 204]/255,'MarkerSize',5,'MarkerEdgeColor',[0,204,204]/255,'MarkerFaceColor',[0,255,255]*(1/255),'linewidth',1.5);grid on;
% ylim([1e-4,1e0])
legend('SUMF','Decorr','MMSE','PIC','SIC','DF-IC','Theory','location','northwest');
xlabel('Number of Users K');ylabel('Average bit error probability')
title('P_{\Delta}=3dB across consecutive users; SNR of user 1 is 6dB.')