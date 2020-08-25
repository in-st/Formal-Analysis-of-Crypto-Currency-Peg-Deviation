% "Estimating overidentified, nonrecursive time varying coefficients SVARs"
% ------------------------------------------------------------------------------------
% This code implements the SVAR model as in Canova and Pérez Forero (2014)
% It contains general options for time variation and sampling in different parameter
% blocks

% This a modified version of the TVP-VAR code from Koop and Korobilis
% http://personal.strath.ac.uk/gary.koop/bayes_matlab_code_by_koop_and_korobilis.html
% ************************************************************************************
% This code has been modified in order to allow for non-recursive identification schemes.
% It uses Sims & Zha (2006)'s identification strategy

% It also incorporates the corrigendum of Del Negro and Primiceri (2013) for
% sampling stochastic volatility

% Modified by
% Fernando Pérez Forero
% fernando.perez@bcrp.gob.pe
% fernandojose.perez@upf.edu
% fernandoperezforero8@gmail.com

% This version: February, 2014.
% -------------------------------------------------------------------------

clear all
close all
clc
warning off all
seed=1;

randn('state',seed);
rand('twister',seed);

disp('Estimating SVAR model')

%% 1. LOAD DATA
%====================================
% US quarterly data

load Data_US_Q

ydata=Data_US_Q(:,1:end-1);
yearlab=Data_US_Q(:,end);

% -------------------
% Select data horizon
% -------------------
% Data availability: 1959.Q1-2005.Q4

% Initial date
% ------------
    yy=1959;   %year
    qq=1;      %quarter     
    ini=yy+round2(qq/4,0.0001);    
  
% Final date
% ------------
    yy=2005;   %year
    qq=4;      %quarter
    fin=yy+round2(qq/4,0.0001);
    
ydata=ydata(find(yearlab==ini):find(yearlab==fin),:);    
yearlab=yearlab(find(yearlab==ini):find(yearlab==fin),:);    

%% 2. Data Transformation 
%=======================================

variables={'GDP','P','U','R','M','Pcom'};
labels=strvcat(variables);

% Pcom : Commodity Price Index
% M    : M2
% R    : Federal Funds Rate (FFR)
% GDP  : Aggregate Gross Domestic Product index (Volume, base 2005=100)
% P    : GDP deflator
% U    : Unemployment Rate

% Transformation codes:
%  1 : None                 : (Y_t)
%  4 : Log Levels           : Log(Y_t) 
%  5 : Quarterly growth     : Log(Y_t)-Log(Y_t_1)
% 20 : Year-to-year growth  : Log(Y_t)-Log(Y_t_4)

%tcodey=[5,5,1,5,5,1];  % Quarterly growth
tcodey=[20,20,1,20,20,1]; % Year-to-year growth

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

% Adjust the number of observations after transforming data

switch(tcodey(1))
    case 1
        ydata=ydata(:,:);
        yearlab=yearlab(:,:);
    case 4
        ydata=ydata(:,:);
        yearlab=yearlab(:,:);
    case 5
        ydata=ydata(2:end,:);
        yearlab=yearlab(2:end,:);
    case 20
        ydata=ydata(5:end,:);
        yearlab=yearlab(5:end,:);
end

% Set variables in the following order:

% GDP  : Aggregate Gross Domestic Product index (Volume, base 2005=100)
% P    : GDP deflator
% U    : Unemployment Rate
% R    : Federal Funds Rate (FFR)
% M    : M2
% Pcom : Commodity Price Index

Y=ydata(:,[4 5 6 3 2 1]);

% Number of observations and variables in Y
t=size(Y,1); % t is the time-series observations of Y
M=size(Y,2); % M is the dimensionality of Y

%% 3. Estimation Setup
% =====================

% 3.1. General options
% --------------------

% Set which parameters are time-varying. The order is:
% 1. Beta: mean (auto)-regression coefficients
% 2. Alpha: Contemporaneous coefficients
% 3. Sigma: Log volatilities
TVP_Beta = 1;   % 0. Constant, 1. Time-varying
TVP_Alpha = 1;  % 0. Constant, 1. Time-varying
TVP_Sigma = 1;  % 0. Constant, 1. Time-varying

tau = 40; % tau is the size of the training sample
crit_lag=4; % 1. FPE. 2. AIC. 3. HQ. 4. SIC.
identification=2; % 1. Cholesky. 2. Sims & Zha (2006)
constant=1; % 0. None. 1: Intercept. 2: Trend+intercept.
test=1; % 1: Short chain. 2: Long chain
prior=1; % 1. ML-based prior. 0. None
df=5;    % Degrees of freedom for t-student proposal distribution
stat=2; % 1: Unrestricted B^T draws % 2: Truncated B^T for stationary draws
sign_restr=0; % 1. Activate Sign Restrictions for MP shocks
lr_restr=1; % 1. Activate Long-run Restrictions for MP shocks
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

disp('   ')

if TVP_Alpha==1
    disp('Block (A^T) is Time varying')
else    
    disp('Block (A^T) is constant')
end

if TVP_Sigma==1
    disp('Block (Sigma^T) is Time varying')
    if diag_sigma==2
        disp('Block (Sigma^T) is GARCH-type')
    end
else
    disp('Block (Sigma^T) is constant')
end

if TVP_Beta==1
    disp('Block (B^T) is Time varying')
    if single_beta==1
        disp('Block (B^T) is Single-Move')
    else
        disp('Block (B^T) is Multi-Move')
    end
    if red_dim==1
        if Xi_rand==1
            disp('Estimating Hierarchical model')
        else
            disp('Estimating Re-parametrized SUR model')        
        end
    end    
else    
    disp('Block (B^T) is constant')
end

disp('   ')

% 3.2. Lag order selection
% ------------------------
p_bar=12;

Sigma_u_tilde=zeros(M,M,p_bar);
IC=zeros(4,p_bar);

DETSIGVEC=zeros(1,p_bar);

for p=1:p_bar

    Z=[ones(1,t-p)];
    for i=1:p
        Z=[Z;Y(p+1-i:t-i,:)'];
    end
    Sigma_u_tilde(:,:,p)=(1/t)*Y(p+1:t,:)'*(eye(t-p)-Z'*(eye(size(Z,1))/(Z*Z'))*Z)*Y(p+1:t,:); 
    IC(1,p)=(((t+M*p+1)/(t-M*p-1))^M)*det(Sigma_u_tilde(:,:,p));                 % FPE
    IC(2,p)=log(det(Sigma_u_tilde(:,:,p)))+2*p*(M^2)/t;                          % AIC
    IC(3,p)=log(det(Sigma_u_tilde(:,:,p)))+2*log(log(t))*p*(M^2)/t;              % HQ
    IC(4,p)=log(det(Sigma_u_tilde(:,:,p)))+log(t)*p*(M^2)/t;                     % SIC
    clear Z
end
 
p_opt=zeros(1,4);
for i=1:4
    p=0;
    diff=IC(i,2)-IC(i,1);
    if diff<0
        while diff<0
             p=p+1;
             diff=IC(i,p+1)-IC(i,p);
        end
        p_opt(i)=p;
    else
        p_opt(i)=p;
    end
end

switch(crit_lag)
    case 1        
        disp('Optimal Lag length using FPE criterion')
        p=p_opt(crit_lag)
    case 2        
        disp('Optimal Lag length using AIC criterion')
        p=p_opt(crit_lag)
    case 3        
        disp('Optimal Lag length using HQ criterion')
        p=p_opt(crit_lag)
    case 4        
        disp('Optimal Lag length using SIC criterion')
        p=p_opt(crit_lag)
    case 5
        disp('Lag length')
        p=2
end

% 3.3. Identification
% -------------------

switch(identification)
    case 1        
        ident_A    = 2*triu(ones(M,M),+1)+eye(M);
        ident_name = {'Cholesky'};
        disp('Identification scheme: Cholesky')
        disp('   ')
    case 2
        load ident_SZ06_alt.prn
        ident_A    = ident_SZ06_alt;
        ident_name = {'SZ06'};
        disp('Identification scheme: Sims & Zha (2006)')
        disp('   ')
end

if sign_restr==1
   disp('Using Sign Restrictions for Monetary Policy shocks')
   disp('   ')
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

if numa+M == M*(M+1)/2
    disp('The SVAR model is exactly identified')
elseif numa+M < M*(M+1)/2
    sprintf('The SVAR model is over-identified by %d restrictions',M*(M+1)/2-(numa+M));
else
    error('The SVAR model is under-identified')
end

index2a=cell(M,1);
index2a2=cell(M,1);
for i=1:M
    index2a{i}=find(ident_A(i,:)==2|ident_A(i,:)==-2);
end

% 3.4. Set the VAR matrices
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
% 3.5. Set the number of draws for simulation (memory!)
% ------------------------------------------------------
switch(test)
    case 1
        nrep     =   2000;  % Number of replications
        nburn    =   1000;  % Number of burn-in-draws
        it_print =    100;  % Print in the screen every "it_print"-th iteration
        nthin    =      4;  % Thinning factor
    case 2
        nrep     =  50000;  % Number of replications
        nburn    = 100000;  % Number of burn-in-draws
        it_print =   1000;  % Print in the screen every "it_print"-th iteration
        nthin    =    100;  % Thinning factor          
end

% 3.6. Select Priors
% ------------------

if prior==1
    % ML estimates with tau observations
    load([sprintf('%s_%d_%d_%d_%d_%d_%d_%d.mat','prior_ML_final',tau,...
       identification,constant,p,standardize,tcodey(1),red_dim)]);
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

% 3.7. Additional settings
% -----------------------
% Parameters of the 7 component mixture approximation to a log(chi^2)
% density: from Kim et al (1998)
q_s = [0.00730; 0.10556; 0.00002; 0.04395; 0.34001; 0.24566; 0.25750];     % probabilities
m_s = [-10.12999; -3.97281; -8.56686; 2.77786; 0.61942; 1.79518; -1.08819];% means
u2_s = [5.79596; 2.61369; 5.17950; 0.16735; 0.64009; 0.34023; 1.26261];    % variances

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
        c_a=0.5*1e-3; % This constant re-scales the variance of proposal density        
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
    
    
    
    % Save a preliminar draft every 10%    
    if test>1
        if mod(irep,0.1*(nrep+nburn))==0
        save([sprintf('%s_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%s_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d.mat',...
            'results_TVP_VAR_MH_revised_draft',M,p,irep,nrep+nburn,nburn,nthin,c_b,c_a,k_Q,k_S,k_W,...
            prior,stat,diag_sigma,char(ident_name),single_beta,TVP_Beta,TVP_Alpha,...
            TVP_Sigma,standardize,tcodey(1),sign_restr,horiz,red_dim,Xi_rand,lr_restr)]);
        end
    end
end %END main Gibbs loop (for irep = 1:nrep+nburn)
clc;
toc; % Stop timer and print total time
%=============================SAMPLING ENDS HERE==================================
%% 5. Print acceptance rates and save output

disp('Estimation of SVAR is done')

disp('   ')

if TVP_Alpha==1
    disp('Block (A^T) is Time varying')
else    
    disp('Block (A^T) is constant')
end

if TVP_Sigma==1
    disp('Block (Sigma^T) is Time varying')
else
    disp('Block (Sigma^T) is constant')
end

if TVP_Beta==1
    disp('Block (B^T) is Time varying')
    if single_beta==1
        disp('Block (B^T) is Single-Move')
    else
        disp('Block (B^T) is Multi-Move')
    end
    if red_dim==1
        if Xi_rand==1
            disp('Hierarchical model')
        else
            disp('Re-parametrized SUR model')        
        end
    end    
else    
    disp('Block (B^T) is constant')
end

disp ('  ')
switch(identification)
    case 1        
        disp('Identification scheme: Cholesky')
    case 2
        disp('Identification scheme: Sims & Zha (2006)')
end

if numa+M == M*(M+1)/2
    disp('The SVAR model is exactly identified')
elseif numa+M < M*(M+1)/2
    sprintf('The SVAR model is over-identified by %d restrictions',M*(M+1)/2-(numa+M));
end

disp ('  ')

disp('Lag length:')
p

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
