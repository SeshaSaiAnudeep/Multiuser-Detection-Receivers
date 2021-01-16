clc;clear variables;close all;
N = 8;P_delta = 3;P=1000;Nsim=1000;K=3;
%delta is difference between two concecutive users
snrdb = 0:12;
for ksnr = 1:length(snrdb)
for ii=1:Nsim
    S = (2*randi([0,1],N,K)-1);
    S = S/sqrt((S(:,1)'*S(:,1)));
    R = S'*S;
    while rank(R)<K
        S = (2*randi([0,1],N,K)-1);
        S = S/sqrt((S(:,1)'*S(:,1)));%S=S/sqrt(N)
        R = S'*S;
    end
    B = 2*randi([0,1],K,P)-1;
    SNRdB = snrdb(ksnr);
    pK_dB = (0:K-1)*P_delta +SNRdB;
    pK_lin = 10.^(pK_dB/10);
    A = sqrt(pK_lin);
    A_inv = 1./(pK_lin);
    A_inv = diag(A_inv);
    A = diag(A);
    n = randn(N,P)+1j*randn(N,P);
    n = n/sqrt(2);
    y = S*A*B+n;
    z = S'*y;
    % SINGLE USER MATCHED FILTER
    b_est_SUMF = sign(real(z));
    ber_SUMF(ksnr,ii) = (1/(K*P))*sum(sum(b_est_SUMF~=B));
    % DECORRELATING RECEIVER
    b_est_Deco = sign(real(R\z));
    ber_Deco(ksnr,ii) = (1/(K*P))*sum(sum(b_est_Deco~=B));
    % MMSE RECEIVER
    b_est_MMSE = sign(real((R+A_inv)\z));
    ber_MMSE(ksnr,ii) = (1/(K*P))*sum(sum(b_est_MMSE~=B));
    % SUCCESSIVE INTERFERNCE CANCELLATION
    b_est_SIC(K,:) = sign(real(z(K,:)));
    for mm = 1:K-1
        r = R(:,mm);
        b_est_SIC(mm,:)= sign(real(z(mm,:)')-B(mm+1:end,:)'*A([mm+1:end], [mm+1:end])*r(mm+1:end));
    end
    ber_SIC(ksnr,ii) = (1/(K*P))*sum(sum(b_est_SIC~=B));
    % PARALLEL INTERFEENCE CANCELLATION
    R_hat = R-eye(K,K);
    b = b_est_SUMF;
    for i=1:10
        b = sign(real(z-R_hat*A*b));
    end
    ber_PIC(ksnr,ii) = (1/(K*P))*sum(sum(b~=B));
    % DF-IC
    FF = chol(R);
    Ft = FF';
    z_hat = (FF')\z;
    b_est_DC(K,:) = sign(real(z_hat(K,:)));
    for mm = 1:K-1
        r = Ft(:,mm);
        b_est_DC(mm,:)= sign(real(z_hat(mm,:)')-B(mm+1:end,:)'*A([mm+1:end], [mm+1:end])*r(mm+1:end));
    end
    ber_DC(ksnr,ii) = (1/(K*P))*sum(sum(b_est_DC~=B));
    ber_theo(ksnr,ii) = (1/K)*sum(qfunc(sqrt(2)*(diag(A))));
end
end
ber_SUMF = sum(ber_SUMF,2)/Nsim;
ber_Deco = sum(ber_Deco,2)/Nsim;
ber_MMSE = sum(ber_MMSE,2)/Nsim;
ber_SIC  = sum(ber_SIC,2)/Nsim;
ber_PIC = sum(ber_PIC,2)/Nsim;
ber_DC = sum(ber_DC,2)/Nsim;
ber_theo = sum(ber_theo,2)/Nsim;

semilogy(snrdb,ber_SUMF,'b-o','MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor',[0.5 .5 1],'linewidth',1.5);hold on;
semilogy(snrdb,ber_Deco,'r->','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor',[1 .5 0.5],'linewidth',1.5);hold on;
semilogy(snrdb,ber_MMSE,'-<','color',[0 0.4 0],'MarkerSize',5,'MarkerEdgeColor',[0 0.5 0],'MarkerFaceColor',[0.6 1 0.6],'linewidth',1.5);grid on;
semilogy(snrdb,ber_PIC,'-^','color',[0 0 0],'MarkerSize',5,'MarkerEdgeColor',[0.25 0.25 0.25],'MarkerFaceColor',[0.55 0.55 0.55],'linewidth',1.5);grid on;
semilogy(snrdb,ber_SIC,'m-v','MarkerSize',5,'MarkerEdgeColor',[1,8/255,127/255],'MarkerFaceColor',[255,113,181]*(1/255),'linewidth',1.5);grid on;
semilogy(snrdb,ber_DC,'-d','color',[84 22 120]/255,'MarkerSize',5,'MarkerEdgeColor',[112,39,195]/255,'MarkerFaceColor',[150,0,205]*(1/255),'linewidth',1.5);grid on;
semilogy(snrdb,ber_theo,'-p','color',[0 204 204]/255,'MarkerSize',5,'MarkerEdgeColor',[0,204,204]/255,'MarkerFaceColor',[0,255,255]*(1/255),'linewidth',1.5);grid on;
% ylim([1e-5,1e-1])
legend('SUMF','Decorr','MMSE','PIC','SIC','DF IC','Theory','location','southwest');
xlabel('SNR for 1^{st} user');ylabel('Average bit error probability')
title('P_{\Delta}=3dB across consecutive users; K=3 users.')