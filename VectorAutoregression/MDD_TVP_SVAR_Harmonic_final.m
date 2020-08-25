% Log-Marginal Data Likelihood for TVP-SVARs

% p(Y)= p(Y|Theta_star)*p(Theta_star)/p(Theta_star|Y)

% March 22th, 2012

% Fernando Pérez Forero
% fernandojose.perez@upf.edu

% This program is for post-estimation analysis.
% It follows Geweke (1999)'s harmonic mean estimator
% It is necessary to first load a workspace with Gibbs-Sampling Output

% This program is very time consuming, it takes 4 min. in a 3GHz 16GB RAM

% ====================================================================

clc

%% 1. Compute posterior mean and variance of hyper-parameters

Ndraws=nrep/nthin;

omega_param=0;

if TVP_Beta==1
    omega_param=omega_param+K*(K+1)/2;
    omega_param=omega_param+K;
end

if TVP_Alpha==1
    omega_param=omega_param+numa*(numa+1)/2;
    omega_param=omega_param+numa;
end

if TVP_Sigma==1
    omega_param=omega_param+M;
    omega_param=omega_param+M;
end

% Mean and Variance

Omega_post=zeros(Ndraws,omega_param);
tic
for k=1:Ndraws
    Bt_temp     = Bt_post{k};
    Bt_temp     = Bt_temp(:,1);
    At_temp     = At_post{k};
    At_temp     = At_temp(:,1);
    Sigt_temp   = Sigt_post{k};
    Sigt_temp     = Sigt_temp(:,1);
    Q_temp=reshape(Q_post(:,k),K,K);
    Q_temp=Q_temp.*tril(ones(K,K));
    Q_temp=reshape(Q_temp,K^2,1);
    indexQ=find(Q_temp~=0);
    Q_temp=Q_temp(indexQ);
    
    Sa_temp=reshape(Sa_post(:,k),numa,numa);
    Sa_temp=Sa_temp.*tril(ones(numa,numa));
    Sa_temp=reshape(Sa_temp,numa^2,1);
    indexSa=find(Sa_temp~=0);
    Sa_temp=Sa_temp(indexSa);     
    
    if omega_param>0    
        Omega_temp=[Bt_temp;At_temp;Sigt_temp;Q_temp;Sa_temp;W_post(:,k)];
        Omega_post(k,:)=Omega_temp';
    end
end   
toc

tic

Omega_hat=mean(Omega_post);
%V_hat=cov(Omega_post);
V_hat=zeros(omega_param,omega_param);
 for k=1:Ndraws
     V_hat=V_hat+(1/Ndraws)*(Omega_post(k,:)-Omega_hat)'*(Omega_post(k,:)-Omega_hat);
 end
toc

%%
% Additional settings
warning off all

%inv_Vhat=eye(omega_param)/V_hat;
inv_Vhat=pinv(V_hat);
ldet_Vhat=-logdet(inv_Vhat);

%% 2. Compute the mean harmonic estimator

mdd=0;
Ndraws_star=nrep/nthin;
h=20;
tic
for k=1:Ndraws_star
    
    R=round(unifrnd(1,Ndraws));
    
    Bt_temp     = Bt_post{R};
    At_temp     = At_post{R};
    Sigt_temp   = Sigt_post{R};
    Qdraw_temp  = reshape(Q_post(:,R),K,K);
    Sadraw_temp = reshape(Sa_post(:,R),numa,numa);        
    Wdraw_temp  = diag(W_post(:,R));
        
    % State space form

    % Y_t = X_t'*B_t  + inv(A_t)*Sigma_t*e_t (1)  e_t ~ N(0,I)
    % B_t = T_t*B_t_1 + u_t                  (2)  u_t ~ N(0,Q)

    % Initialize the Filter
    LogL_Yk = 0;
    Bt_t_1 = Bt_temp(:,1);
    Pt_t   = VB_OLS;
    T_t    = eye(K);
    Pt_t_1 = T_t*Pt_t*T_t'+Qdraw_temp;
    Z_B    = Z(1:M,:);        % From measurement equation
    Var_e  = (reshape(S_A*At_temp(:,1)+s_A,M,M))\diag(Sigt_temp(:,1));
    Var_e  = Var_e*Var_e';
    P_inv  = eye(size(Z_B,1))/(Z_B*Pt_t_1*Z_B'+Var_e);
    K_B    = Pt_t_1'*Z_B'*P_inv; % Kalman gain
    Pt_t   = Pt_t_1-K_B*(Z_B*Pt_t_1*Z_B'+Var_e)*K_B';
    Pred_e = y(:,1)-Z_B*Bt_t_1;
    Bt_t   = Bt_t_1 + K_B*Pred_e;
    LogL_Yk = LogL_Yk-0.5*M*log(2*pi)+0.5*logdet(P_inv)...
         -0.5*Pred_e'*P_inv*Pred_e;

    % Iterate

    for i=2:t
    % Prediction:        
      Pt_t_1 = T_t*Pt_t*T_t'+Qdraw_temp;
      Bt_t_1 = T_t*Bt_t;        
    % Updating:
      Var_e  = (reshape(S_A*At_temp(:,i)+s_A,M,M))\diag(Sigt_temp(:,i));
      Var_e  = Var_e*Var_e';
      Z_B    = Z((i-1)*M+1:i*M,:);
      P_inv  = eye(size(Z_B,1))/(Z_B*Pt_t_1*Z_B'+Var_e);
      K_B    = Pt_t_1'*Z_B'*P_inv; % Kalman gain
      Pt_t   = Pt_t_1-K_B*(Z_B*Pt_t_1*Z_B'+Var_e)*K_B';
      Pred_e = y(:,i)-Z_B*Bt_t_1;
      Bt_t   = Bt_t_1 + K_B*Pred_e;
      LogL_Yk = LogL_Yk-0.5*M*log(2*pi)+0.5*logdet(P_inv)...
         -0.5*Pred_e'*P_inv*Pred_e;         
    end
    
    lB_prpdf=-0.5*K*log(2*pi)...
             -0.5*logdet(B_0_prvar)...
             -0.5*(Bt_temp(:,1)-B_0_prmean)'*inv(B_0_prvar)*(Bt_temp(:,1)-B_0_prmean);
    
    lA_prpdf=-0.5*numa*log(2*pi)...
             -0.5*logdet(A_0_prvar)...
             -0.5*(At_temp(:,1)-A_0_prmean)'*inv(A_0_prvar)*(At_temp(:,1)-A_0_prmean);
         
    lSig_prpdf=-0.5*K*log(2*pi)...
             -0.5*logdet(B_0_prvar)...
             -0.5*(log(2*Sigt_temp(:,1))-sigma_prmean)'*inv(sigma_prvar)*(log(2*Sigt_temp(:,1))-sigma_prmean);         
         
         
    if TVP_Beta==1
        [Q_prpdf lQ_prpdf]=inv_wishpdf(Qdraw_temp,Q_prmean,Q_prvar);
    else
        Q_prpdf=0;
        lQ_prpdf=1;
    end
    
    if TVP_Alpha==1
        [Sa_prpdf lSa_prpdf]=inv_wishpdf(Sadraw_temp,Sa_prmean,Sa_prvar);
    else
        Sa_prpdf=0;
        lSa_prpdf=1;
    end
    
    W_prpdf=1;
    lW_prpdf=0;
    
    if TVP_Sigma==1
        for l=1:M
        [W_prpdft lW_prpdft]=inv_wishpdf(Wdraw_temp(l,l),W_prmean,W_prvar);
        W_prpdf=W_prpdf*W_prpdft;
        lW_prpdf=lW_prpdf+lW_prpdft;    
        end
    end
    prior_pdf=Q_prpdf*Sa_prpdf*W_prpdf;    
    lprior_pdf=lB_prpdf+lA_prpdf+lSig_prpdf+lQ_prpdf+lSa_prpdf+lW_prpdf;
    
    Q_temp=reshape(Q_post(:,R),K,K);
    Q_temp=Q_temp.*tril(ones(K,K));
    Q_temp=reshape(Q_temp,K^2,1);
    indexQ=find(Q_temp~=0);
    Q_temp=Q_temp(indexQ);
    
    Sa_temp=reshape(Sa_post(:,R),numa,numa);
    Sa_temp=Sa_temp.*tril(ones(numa,numa));
    Sa_temp=reshape(Sa_temp,numa^2,1);
    indexSa=find(Sa_temp~=0);
    Sa_temp=Sa_temp(indexSa); 
        
    Omega_temp=[Bt_temp(:,1);At_temp(:,1);Sigt_temp(:,1);Q_temp;Sa_temp;W_post(:,R)];
    tau=0.9;
    
    if omega_param>0
        exp_omega=(Omega_temp-Omega_hat')'*inv_Vhat*(Omega_temp-Omega_hat');
    else
        exp_omega=0;
    end
    
    if exp_omega<=chi2inv(tau,omega_param)
        lg_omega=-0.5*exp_omega-0.5*ldet_Vhat...
            -log(tau)-(omega_param/2)*log(2*pi);
    else
        lg_omega=1;
    end
    
    mdd=mdd+exp_approx(lg_omega,0,h)/(exp(LogL_Yk)*exp_approx(lprior_pdf,0,h));    
    
end   
toc
mdd=(mdd/Ndraws_star)^(-1);

disp('Log-MDD')
lmdd=log(mdd)
      