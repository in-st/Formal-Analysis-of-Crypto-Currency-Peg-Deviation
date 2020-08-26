% This code produces the SVAR estimates.
% This code is based on "Estimating overidentified, non-recursive, time varying coefficients structural VARs"... 
% ...by Fabio Canova and Fernando J. Pérez Forero (Quantitative Economics, Volume 6, Issue 2 (July 2015)).

% Modified to run on GNU Octave 5.2.0 by
% Huachen Li
% -------------------------------------------------------------------------
pkg load statistics
clear all
close all
clc
warning off all
seed=10101010101010101;
randn('state',seed);
rand('twister',seed);

%% 1. LOAD DATA
%====================================
load Data  %1:Peg Deviation (%); 2:dln(BTC); 3:dln(volume); 4:Fees Per Transaction (USD)

% Sample range: 4/14/2017 - 7/1/2020

ydata=data(:,1:end);
date = 2017+104/365:1/365:2020+182/365;;
yearlab = date';
variables={'x','p','v','f'};
labels=strvcat(variables);

tcodey=[1,1,1,1];

for i_y = 1:size(ydata,2)
    ydata(:,i_y) = transx(ydata(:,i_y),tcodey(i_y));
end

% Scaling data
t2 = size(ydata,1);
standardize=1;

switch(standardize)
    case 1 % Standardized
        ydata = (ydata- repmat(mean(ydata,1),t2,1))./repmat(std(ydata,1),t2,1);
    case 2 % Non-standardized
        ydata(:,[3,6])=0.01*ydata(:,[3,6]);
end

Y=ydata(:,[4 2 1 3]); %4:Fees Per Transaction (USD); 2:dln(BTC); 1:Peg Deviation (%);  3:dln(volume); 

% Number of observations and variables in Y
t=size(Y,1); % t is the time-series observations of Y
M=size(Y,2); % M is the dimensionality of Y

%% 2. Estimation Setup
% =====================

% 2.1. General options
% --------------------
TVP_Beta = 1;   % 0. Constant, 1. Time-varying
TVP_Alpha = 1;  % 0. Constant, 1. Time-varying
TVP_Sigma = 1;  % 0. Constant, 1. Time-varying
cc = 1.01;  % Beta volatility
tau = 90; % tau is the size of the training sample
p=2; % Lag order
ident_A    = 2*triu(ones(M,M),+1)+eye(M);
ident_name = {'Cholesky'};
constant=1; % 0. None. 1: Intercept. 2: Trend+intercept.
nrep=200;  % Number of replications
nburn=200;  % Number of burn-in-draws
it_print=100;  % Print in the screen every "it_print"-th iteration
nthin=1;  % Thinning factor
prior=1; % 1. ML-based prior. 0. None
df=5;    % Degrees of freedom for t-student proposal distribution
stat=2; % 1: Unrestricted B^T draws % 2: Truncated B^T for stationary draws
sign_restr=0; % 1. Activate Sign Restrictions for MP shocks
lr_restr=0; % 1. Activate Long-run Restrictions for MP shocks
draft=0; % 1. Start from a previous draft. 0. Start from the beginning

if TVP_Beta==1
    single_beta=2; % 1. Single Move. 2.Multi-Move.
    red_dim=0; % 1. Reduced-dimension TVP-SVAR
    Xi_rand=0; % 1. Hierarchical model. 0. Re-parametrized SUR model
else
    single_beta=0; % 1. Single Move. 2.Multi-Move.
    red_dim=0; % 1. Reduced-dimension TVP-SVAR
    Xi_rand=0; % 1. Hierarchical model. 0. Re-parametrized SUR model    
end

if TVP_Sigma==1
    diag_sigma=1; % 1. Diagonal Sigma, 2. Garch    
else
    diag_sigma=0; % 1. Diagonal Sigma
end

index_smallSA=find(ident_A==1|ident_A==-1);
s_A=zeros(size(ident_A,1),size(ident_A,2));
s_A(index_smallSA)=ident_A(index_smallSA);
s_A=reshape(s_A,size(ident_A,1)*size(ident_A,2),1);
index_bigSA=find(ident_A==2|ident_A==-2);
S_A=zeros(size(ident_A,1)*size(ident_A,2),size(index_bigSA,1));

for i=1:size(index_bigSA,1)
    S_A(index_bigSA(i,1),i)=ident_A(index_bigSA(i,1));
end

if isempty(index_bigSA)==1
else
    S_A=S_A/2;
end

numa =numel(index_bigSA);

index2a=cell(M,1);
index2a2=cell(M,1);
for i=1:M
    index2a{i}=find(ident_A(i,:)==2|ident_A(i,:)==-2);
end

% 2.2. Set the VAR matrices
% -------------------------

ylag = mlag2(Y,p); % Lags
% Form RHS matrix X_t = [1 y_t-1 y_t-2 ... y_t-k] for t=1:T
ylag = ylag(p+tau+1:t,:);

K = M*(constant) + p*(M^2); % K is the number of elements in the state vector

% Create Z_t matrix.
Z = zeros((t-tau-p)*M,K);

for i = 1:t-tau-p
    switch(constant)
        case 0
            ztemp = [];
        case 1
            ztemp = eye(M);
        case 2
            ztemp = [eye(M),(i+p)*eye(M)];        
    end
    
    for j = 1:p        
        xtemp = ylag(i,(j-1)*M+1:j*M);
        xtemp = kron(eye(M),xtemp);
        ztemp = [ztemp xtemp];
    end
    Z((i-1)*M+1:i*M,:) = ztemp;
end

% Set the final sample for estimation

y = Y(tau+p+1:t,:)';
yearlab = yearlab(tau+p+1:t);
t=size(y,2);   %Time series observations

if TVP_Beta==1
    if red_dim==1
        if Xi_rand==1 % Hierarchical Model
            % Set number of factors
            dim_Theta=15;
            Xi_draw=zeros(K,dim_Theta);        
            Xi_draw(1:dim_Theta,1:dim_Theta)=eye(dim_Theta);        
        else % Re-parametrized SUR model
            % Factors' setup
    
            % Intercepts
    
             Xi0=[eye(M);
              zeros(K-M,M)];
      
            % Variable effect
      
            Xi1=zeros(M,M);
    
            temp_0=eye(M);

            for i=1:M*p
                Xi1=[Xi1;temp_0];
            end

            % Lag effect
    
            Xi2=[zeros(M,p);
                 kron(eye(p),ones(M^2,1))];
     
            % Common component
     
            Xi3=[zeros(M,1);ones(K-M,1)];

            % Equation effect
    
            Xi4=zeros(M,M);     
     
            for i=1:p
                  Xi4=[Xi4;kron(eye(M),ones(M,1))];    
            end     
    
            Xi_draw=[Xi0,Xi1,Xi2,Xi3,Xi4];
            dim_Theta=size(Xi_draw,2);    
        end
    end
end

% 3.3. Select Priors
% ------------------

if prior==1
    identification=1; 
    % ML estimates with tau observations
     load([sprintf('%s_%d_%d_%d_%d_%d_%d_%d.mat','prior_ML_final',tau,...
        identification,constant,p,standardize,tcodey(1),red_dim)]);
     % load('prior');
    if TVP_Beta==1 && Xi_rand==1
        Theta_OLS   =   zeros(dim_Theta,1);
        VTheta_OLS  =   eye(dim_Theta);        
    end
else
    % Uninformative values
    
    if TVP_Beta==1 && red_dim==1
        B_OLS = zeros(K,1);
        VB_OLS = eye(K);
   
        Theta_OLS   =   zeros(dim_Theta,1);
        VTheta_OLS  =   eye(dim_Theta);
        if Xi_rand==1
            Xi_0_prmean = zeros(K,dim_Theta);
            k_Xi=sqrt(1e-4);
            Xi_0_prvar  = (k_Xi)^2*eye(dim_Theta);
        end
    else   
        B_OLS     = zeros(K,1);
        VB_OLS    = eye(K);
    end
    
    A_OLS     = zeros(numa,1);
    VA_OLS    = eye(numa);    
    sigma_OLS = zeros(M,1);    
end

% Set some hyperparameters here
if TVP_Beta==1 
    if red_dim==1
        if Xi_rand==1
            k_Q = sqrt(1e-4/2);
            k_R = sqrt(1e-2);
        else
            k_Q = sqrt(1e-4/2);
        end
    else
        if lr_restr==1
            k_Q = sqrt(1e-4/2);        
            %k_Q = sqrt(1e-5);        
            %k_Q = sqrt(1e-4);
        else
            k_Q = sqrt(1e-4/2);        
        end        
    end
else
   k_Q = 0; 
end

if TVP_Alpha==1 
    if lr_restr==1    
        k_S = sqrt(1e-3);
    else
        k_S = sqrt(1e-3);
    end
else
    k_S = 0;
end
if TVP_Sigma==1     
    if lr_restr==1
        k_W = sqrt(1e-3);
    else
        k_W = sqrt(1e-4);
    end
else
    k_W = 0;
end

% Prior hyperparameters
if red_dim==1
    if Xi_rand==1
    sizeQ = K; % Size of matrix Q
    sizeR = dim_Theta; % Size of matrix Q
    else
        sizeQ =dim_Theta; % Size of matrix Q
    end
else
    sizeQ = K; % Size of matrix Q
end
sizeW = M; % Size of matrix W

% Set prior means and variances (_prmean / _prvar)
% These are the Kalman filter initial conditions for the time-varying
% parameters B(t), A(t) and (log) SIGMA(t).

% VAR Coefficients
% B_0 ~ N(B_OLS, 4Var(B_OLS))
if red_dim==1
    
    if Xi_rand==1
       Xi_0_prmean = zeros(K,dim_Theta);
       k_Xi=sqrt(1e-4);
       Xi_0_prvar  = (k_Xi)^2*eye(dim_Theta);       
       Theta_0_prmean = Theta_OLS;
       Theta_0_prvar = VTheta_OLS;       

    else
       Theta_0_prmean = Theta_OLS;
       Theta_0_prvar = VTheta_OLS;          
    end    
else
    B_0_prmean = B_OLS;
    B_0_prvar = 4*VB_OLS;
end

Beta_counter=0; % Counter for stationary draws

% Contemporaneous coefficients matrix
% A_0 ~ N(A_OLS, I_n)
A_0_prmean = A_OLS;
A_0_prvar = eye(numa);

% Shock Volatilities
if diag_sigma<2 % Stochastic volatility
    % log(sigma_0) ~ N(log(sigma_OLS),I_n)
    sigma_prmean = log(sigma_OLS);
    if lr_restr==1
        sigma_prvar = 10*eye(M);
    else
        sigma_prvar = 10*eye(M);
    end
else % Garch Model
    % sigma_0 ~ N(sigma_OLS,I_n)
    sigma_prmean = sigma_OLS;
    sigma_prvar = 10*eye(M);
end

% Q is the covariance of B(t), S is the covariance of A(t) and W is the
% covariance of (log) SIGMA(t)

if TVP_Beta==1

    % Q ~ IW(k2_Q*Var(B_OLS),1+dim(B))
    if red_dim==1
        if Xi_rand==1
            Q_prmean = ((k_Q).^2)*diag(diag(VB_OLS));        
            R_prmean = ((k_R).^2)*VTheta_OLS;
    
            Q_prvar = 1 + 1;
            R_prvar = 1 + 1;    
        else
            Q_prmean = ((k_Q).^2)*VTheta_OLS;        
            Q_prvar = 1 + 1;
        end
    else
        Q_prmean = ((k_Q).^2)*VB_OLS;
        Q_prvar = 1 + sizeQ;
    end
end

if TVP_Alpha==1
    % S ~ IW(k2_S*Var(A_OLS),(1+dim(alpha)))
    Sa_prmean = ((k_S)^2)*VA_OLS;
    Sa_prvar = 1 + numa;
end

if TVP_Sigma==1
    % W ~ IW(k2_W*I_n,(1+dim(sigma)))

    switch(diag_sigma)
        case 1
            % Diagonal Var-cov matrix
            W_prmean = 1*((k_W)^2);
            W_prvar = 1 + 1;
        case 2
            % Unrestricted Var-cov matrix
            W_prmean = ((k_W)^2)*eye(M);
            W_prvar = 1 + sizeW;
    end
end

% 3.4. Additional settings
% -----------------------
%Omori et al (2007)
q_s = [0.00609; 0.04775; 0.13057; 0.20674; 0.22715; 0.18842; 0.12047; 0.05591; 0.01575; 0.00115];     % probabilities
m_s = [1.92677; 1.34744; 0.73504; 0.02266; -0.85173; -1.97278; -3.46788; -5.55246; -8.68384; -14.65];% means
u2_s = [0.11265; 0.17788; 0.26768; 0.40611; 0.62699; 0.98583; 1.57469; 2.54498; 4.16591; 7.33342];    % variances
%========= INITIALIZE MATRICES:
if draft==1 % Start from a previous work
    
    perc_start=0.5;
    irep=perc_start*(nrep+nburn);
    
    % Gibbs sampling settings (B^T block)
    % ----------------------------------------

    if TVP_Beta==1 && single_beta==1
        c_b=1;% This constant re-scales the variance of proposal density of B(t)
        Rsimul=25; % Koop and Potter (2011)
    else
        c_b=0;
        rho_b=1;
    end
    
    % Metropolis Hastings Settings (A^T block)
    % ----------------------------------------
    if TVP_Alpha==1
        c_a=1e-3; % This constant re-scales the variance of proposal density        
    else
        c_a=1e-3; % This constant re-scales the variance of proposal density        
    end

    if sign_restr==1
        horiz=3;
    else
        horiz=0;
    end    
    
    load([sprintf('%s_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%s_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d.mat',...
            'results_TVP_VAR_MH_revised_draft',M,p,irep,nrep+nburn,nburn,nthin,c_b,c_a,k_Q,k_S,k_W,...
            prior,stat,diag_sigma,char(ident_name),single_beta,TVP_Beta,TVP_Alpha,...
            TVP_Sigma,standardize,tcodey(1),sign_restr,horiz,red_dim,Xi_rand,lr_restr)]);
        
    it_start=irep;
else % Start from the beginning
    it_start=1;
    % Specify covariance matrices for measurement and state equations
    consQ = 0.0001;
    consSa = 0.0001;
    consH = 0.0001;
    consW = 0.01;
    Ht = kron(ones(t,1),consH*eye(M));   % Initialize Htdraw, a draw from the VAR covariance matrix
    Htsd = kron(ones(t,1),sqrt(consH)*eye(M)); % Cholesky of Htdraw defined above
    Htchol = kron(ones(t,1),sqrt(consH)*eye(M)); % Cholesky of Htdraw defined above
    Qdraw = consQ*eye(sizeQ);   % Initialize Qdraw, a draw from the covariance matrix Q
    Sadraw = consSa*eye(numa);  % Initialize Sdraw, a draw from the covariance matrix S
    Wdraw = consW*eye(M);    % Initialize Wdraw, a draw from the covariance matrix W
    if red_dim==1
        Thetatdraw = zeros(dim_Theta,t);     % Initialize Thetatdraw, a draw of the mean factor coefficients, Theta(t)
        Btdraw = Xi_draw*Thetatdraw; 
        if Xi_rand==1
            Rdraw = consQ*eye(sizeR);   % Initialize Rdraw, a draw from the covariance matrix R
        end
    else
        Btdraw = zeros(K,t);     % Initialize Btdraw, a draw of the mean VAR coefficients, B(t)
    end
    Atdraw = zeros(numa,t);  % Initialize Atdraw, a draw of the non 0 or 1 elements of A(t)
    if TVP_Alpha==0
       alpha_draw= zeros(numa,1);
       P_draw=eye(numa);       % Initial value of var(alpha) (unconditional = 1)    
    end
    
    % Recover matrix A(t)
    
    capAt = zeros(M*t,M);
    for i = 1:t
        capAt((i-1)*M+1:i*M,:) = reshape(S_A*Atdraw(:,i)+s_A,M,M);
    end    
    
    if diag_sigma<2
        Sigtdraw = zeros(M,t);   % Initialize Sigtdraw, a draw of the log-diagonal of SIGMA(t)
    else
        Sigtdraw = 0.5*ones(M,t);   % Initialize Sigtdraw, a draw of the diagonal of SIGMA(t)
    end
    sigt = kron(ones(t,1),0.01*eye(M));   % Matrix of the exponent of Sigtdraws (SIGMA(t))
    statedraw = round(unifrnd(1,M))*ones(t,M);       % initialize the draw of the indicator variable 
                               % (of 7-component mixture of Normals approximation)
    Zs = kron(ones(t,1),eye(M));
    prw = zeros(numel(q_s),1);
      
    bigj = zeros(M,M*p);
    bigj(1:M,1:M) = eye(M);    

    % Storage matrices for posteriors and stuff

    %----------------
    % Store results
    %----------------
    
    At_post=cell(1,nrep/nthin);
    Bt_post=cell(1,nrep/nthin);
    Sigt_post=cell(1,nrep/nthin);
    Ht_post=cell(1,nrep/nthin);
    if TVP_Alpha==1
        alpha_a_ekf=zeros(numa,t);
        alpha_a_ekf_sm=zeros(numa,t);
        Pt_t_1_a=zeros(numa,numa,t);
        P_ekf_a=zeros(numa,numa,t);
        P_ekf_sm_a=zeros(numa,numa,t);
        P_ekf_sm_star_a=zeros(numa,numa,t);
    end

    
    if TVP_Beta==1 && red_dim==1
       Thetat_post=cell(1,nrep/nthin);
       Xi_post=cell(1,nrep/nthin);
    end

    % Pre-allocating memory
    % ---------------------

    for i=1:nrep/nthin
        At_post{i}=Atdraw;
        Bt_post{i}=Btdraw;
        Sigt_post{i}=Sigtdraw;     
        if TVP_Beta==1 && red_dim==1
            Thetat_post{i}=Thetatdraw;
            if Xi_rand==1
                Xi_post{i}=Xi_draw;
                R_post=zeros(sizeR^2,nrep/nthin);
            end
        end
        Ht_post{i}=Ht;
    end

    if TVP_Beta==1    
        Q_post=zeros(sizeQ^2,nrep/nthin);
    end
    if TVP_Alpha==1
        Sa_post=zeros(numa^2,nrep/nthin);
    end
    if TVP_Sigma==1
        W_post=zeros(M,nrep/nthin);
    end

    % Gibbs sampling settings (B^T block)
    % ----------------------------------------
    if TVP_Beta==1 && single_beta==1
        c_b=1;% This constant re-scales the variance of proposal density of B(t)
        Rsimul=25; % Koop and Potter (2011)
        Beta_counter=zeros(1,t); % Counter for stationary draws
        Q_counter=0;
    else
        c_b=0;
        rho_b=1;
    end

    % Metropolis Hastings Settings (A^T block)
    % ----------------------------------------
       
    acc_counter_a=0;  % Counter for accepted draws
    out_counter=0;    % Counter for out-of-bounds draws
    if TVP_Alpha==1
        c_a=0.1*1e-3; % This constant re-scales the variance of proposal density        
    else
        c_a=1e-3; % This constant re-scales the variance of proposal density        
    end

    if sign_restr==1
        horiz=3;
        Bt_old=Btdraw;
        At_old=Atdraw;
        Sigt_old=Sigtdraw;
        sign_counter=0;
    else
        horiz=0;
    end

    if lr_restr==1
        Bt_old=Btdraw;
        At_old=Atdraw;
        Sigt_old=Sigtdraw;
        Ht_old=Ht;
        Q_old=Qdraw;
        Sa_old=Sadraw;
        W_old=Wdraw;        
        lr_counter=0;
    end    
    
    Bound=20;
    if TVP_Alpha==1
        LB =-Bound*ones(numa,t);
        UB = Bound*ones(numa,t);
    else
        LB =-Bound*ones(numa,1);
        UB = Bound*ones(numa,1);        
    end

    if diag_sigma==2    
        if TVP_Sigma==1
            LB_s = zeros(M,t);
            UB_s = Bound*ones(M,t);
            acc_counter_s=zeros(M,1);  % Counter for accepted draws
            out_counter_s=zeros(M,1);  % Counter for accepted draws                        
            c_s=0.5*1e-3; % This constant re-scales the variance of proposal density               
            delta=0.5; % Geweke & Tanizaki (2001)            
        else
            LB_s =-Bound*ones(M,1);
            UB_s = Bound*ones(M,1);        
        end
    else
        delta=0;
    end
end

%% 4. Metropolis within Gibbs Sampling
%=======================================
tic;
disp('Number of iterations');

rep=0;

% Main loop
for irep = it_start:nrep + nburn
    
    % Print iterations
    if mod(irep,it_print) == 0
        disp(irep);toc;
    end
    % -----------------------------------------------------------------------------------------
    %   STEP I: Sample B from p(B|y,A,Sigma,V)
    % -----------------------------------------------------------------------------------------
    
    if TVP_Beta==1    
        if red_dim==1
            if Xi_rand==1
                draw_beta_conditional_hierarchical
            else
                draw_beta_conditional_reparametrized
            end
        else    
            switch(single_beta)
                case 1            
                    draw_beta_conditional_single
                case 2
                    if lr_restr==1
                       draw_beta_conditional_lr
                    else
                       draw_beta_conditional
                    end
            end     
        end
    else
        draw_beta_conditional_constant
    end
        
    %-------------------------------------------------------------------------------------------
    %   STEP II: Draw A(t) from p(At|y,B,Sigma,V)
    %-------------------------------------------------------------------------------------------
    
    if TVP_Alpha==1
        draw_alpha_a_MH_TVP_dec13
    else
        draw_alpha_a_MH_constant
    end

    %------------------------------------------------------------------------------------------
    %   STEP III: Draw diagonal VAR covariance matrix SIGMA(t)
    %------------------------------------------------------------------------------------------
    
    if TVP_Sigma==1
        if diag_sigma==2
            draw_sigma_TVP_corrected_GARCH2
        else
            draw_sigma_TVP_corrected
        end
    else
        draw_sigma_constant_corrected
    end
    
    %------------------------------------------------------------------------------------------
    %   STEP IV: Evaluate Sign Restrictions
    %------------------------------------------------------------------------------------------    
    
    if sign_restr==1
       draw_sign_restrictions
    end
           
    %------------------------------------------------------------------------------------------
    %   STEP V: Draw the Reduced-Form Covariance Matrix H(t)
    %------------------------------------------------------------------------------------------        

    Ht = zeros(M*t,M);
    Htsd = zeros(M*t,M);
    for i = 1:t                
        stem = sigt((i-1)*M+1:i*M,:);
        if red_dim==1 && Xi_rand==0
            H_temp = stem;
        else
            a = capAt((i-1)*M+1:i*M,:);
            H_temp = a\stem;
        end
        
        Hsd=H_temp;        
        Hdraw = Hsd*Hsd';
        Ht((i-1)*M+1:i*M,:) = Hdraw;  % H(t) : Var-cov of reduced-form VAR
        Htsd((i-1)*M+1:i*M,:) = Hsd;  
    end
    
    %----------------------------SAVE AFTER-BURN-IN DRAWS AND IMPULSE RESPONSES -----------------
    if irep > nburn && mod((irep-nburn),nthin)==0 
        rep=rep+1;    
        Bt_post{(irep-nburn)/nthin} = Btdraw; % regression coefficients B(t)
        At_post{(irep-nburn)/nthin} = Atdraw; % contemp- coeffs matrix A(t)
        if diag_sigma<2
            Sigt_post{(irep-nburn)/nthin} = exp(0.5*Sigtdraw); % diagonal std matrix SIGMA(t)
        else
            Sigt_post{(irep-nburn)/nthin} = Sigtdraw; % diagonal std matrix SIGMA(t)
        end
        if TVP_Beta==1            
            Q_post(:,(irep-nburn)/nthin)=reshape(Qdraw,sizeQ^2,1); % covariance matrix Q of B(t)        
        end
        if TVP_Alpha==1
            Sa_post(:,(irep-nburn)/nthin)=reshape(Sadraw,numa^2,1); % covariance matrix S of A(t)
        end        
        if TVP_Sigma==1        
            W_post(:,(irep-nburn)/nthin)=diag(Wdraw); % covariance matrix W of Sigma(t)
        end        
        if TVP_Beta==1 && red_dim==1
            Thetat_post{(irep-nburn)/nthin} = Thetatdraw; % factors Theta(t)
            if Xi_rand==1
                Xi_post{(irep-nburn)/nthin} = Xi_draw; % Loadings Xi
                R_post(:,(irep-nburn)/nthin)=reshape(Rdraw,sizeR^2,1); % covariance matrix R of Theta(t)        
            end
        end
        
        Ht_post{(irep-nburn)/nthin} = Htsd;
        
        %----------------
        % Get time-varying correlations and variances
        stemp6 = zeros(M,1);
        stemp5 = [];
        stemp7 = [];
        for i = 1:t
            stemp8 = corrvc(Ht((i-1)*M+1:i*M,:));
            stemp7a = [];
            ic = 1;
            for j = 1:M
                if j>1;
                    stemp7a = [stemp7a ; stemp8(j,1:ic)'];
                    ic = ic+1;
                end
                stemp6(j,1) = sqrt(Ht((i-1)*M+j,j));
            end
            stemp5 = [stemp5 ; stemp6'];
            stemp7 = [stemp7 ; stemp7a'];
        end
    end % END saving after burn-in results 
 
end %END main Gibbs loop (for irep = 1:nrep+nburn)
clc;
toc; % Stop timer and print total time
%=============================SAMPLING ENDS HERE==================================
%% 5. Print acceptance rates and save output

disp('Estimation of SVAR is done')

disp('   ')

if numa+M == M*(M+1)/2
    disp('The SVAR model is exactly identified')
elseif numa+M < M*(M+1)/2
    sprintf('The SVAR model is over-identified by %d restrictions',M*(M+1)/2-(numa+M));
end

disp ('  ')


sprintf('Average Acceptance rate in Metropolis Step sampling (A^T): %s percent',num2str(mean(100*acc_counter_a/(nrep + nburn))))

sprintf('Rate of out-of-bounds draws in Metropolis Step sampling (A^T): %s percent',num2str(mean(100*out_counter/(nrep + nburn))))

sprintf('Average Acceptance rate in Gibbs Sampling (stationary draws of B^T): %s percent',num2str(mean(100*Beta_counter/(nrep + nburn))))

if diag_sigma==2
    sprintf('Average Acceptance rate in Metropolis Step sampling (Sigma^T): %s percent',num2str(mean(100*acc_counter_s/(nrep + nburn))))
end

if single_beta==1
    sprintf('Acceptance rate in Gibbs Sampling Q: %s percent',num2str(mean(100*Q_counter/(nrep + nburn))))
end

if sign_restr==1
    sprintf('Acceptance rate in Sign Restrictions: %s percent',num2str(mean(100*sign_counter/(nrep + nburn))))
end

%%
save([sprintf('%s_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%s_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d.mat',...
    'results_TVP_VAR_MH_revised',M,p,nrep+nburn,nburn,nthin,c_b,c_a,k_Q,k_S,k_W,...
    prior,stat,diag_sigma,char(ident_name),single_beta,TVP_Beta,TVP_Alpha,...
    TVP_Sigma,standardize,tcodey(1),sign_restr,horiz,red_dim,Xi_rand,df,delta,rho_b,lr_restr)]);
