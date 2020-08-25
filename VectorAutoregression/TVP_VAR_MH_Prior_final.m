% "Estimating overidentified, non-recursive, time varying coefficients structural VARs"

% Computes the ML estimates for the Training sample using different
% starting values.

% The code uses both Cholesky and Sims & Zha (2006)'s identification strategy

% Modified by
% Fernando Pérez Forero
% fernando.perez@bcrp.gob.pe
% fernandojose.perez@upf.edu
% fernandoperezforero8@gmail.com

% This version: February, 2014
% -------------------------------------------------------------------------

clear all;
clc;
seed=12; rand( 'state' , seed ); randn('state',seed ); % set the random seed

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
% P    : GDP de?ator
% U    : Unemployment Rate

% Transformation codes:
%  1 : None             : (Y_t)
%  4 : Log Levels       : Log(Y_t) 
%  5 : Quarterly growth : Log(Y_t)-Log(Y_t_1)
% 20 : Annual growth    : Log(Y_t)-Log(Y_t_4)

%tcodey=[5,5,1,5,5,1];
tcodey=[20,20,1,20,20,1];

for i_y = 1:size(ydata,2)
    ydata(:,i_y) = transx(ydata(:,i_y),tcodey(i_y));
end

% Demean and standardize data
t2 = size(ydata,1);
scales=std(ydata,1);
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

Y=ydata(:,[4 5 6 3 2 1]);

% Number of observations and variables in Y
t=size(Y,1); % t is the time-series observations of Y
M=size(Y,2); % M is the dimensionality of Y


%% 3. Estimation Setup
% =====================

% 3.1. General options
% --------------------

tau = 40; % tau is the size of the training sample
crit_lag=5; % 1. FPE. 2. AIC. 3. HQ. 4. SIC.
identification=2; % 1. Cholesky. 2. Sims & Zha (2006)
constant=1; % 0. None. 1: Intercept. 2: Trend+intercept.
red_dim=0; % Reduced-dimension as in Canova and Ciccarelli (2009)
Nbar=10; % Number of initial points

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
        p_opt(i)=p+1;
    else
        p_opt(i)=p+1;
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
    case 2
        load ident_SZ06_alt.prn
        ident_A    = ident_SZ06_alt;
        ident_name = {'SZ06'};
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

if red_dim==1
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
    
    Xi=[Xi0,Xi1,Xi2,Xi3,Xi4];

    dim_Theta=size(Xi,2);

end

%% 3. Maximum Likelihood Estimation
% =================================

Y_prior=Y(1:p+tau,:);

reduced_form = vare(Y_prior,p); % Using spatial-econometrics' routine
t_prior=reduced_form(1).nobs;

B_OLS1=[];
for i=1:M
    B_OLS1=[B_OLS1,reduced_form(i).beta];
end

B_OLS=[];
B_OLS=[B_OLS;B_OLS1(end,:)'];
B_OLStemp=B_OLS1(1:end-1,:);

if mod(p,2)==0
   B_OLStemp=B_OLStemp([(1:2:M*p-1),(2:2:M*p)],:);
else
   B_OLStemp=B_OLStemp([(1:2:M*p),(2:2:M*p-1)],:);
end

B_OLS=[B_OLS;reshape(B_OLStemp(1:fix(M*p/2),:),numel(B_OLStemp(1:fix(M*p/2),:)),1)];
B_OLS=[B_OLS;reshape(B_OLStemp(fix(M*p/2)+1:end,:),numel(B_OLStemp(fix(M*p/2)+1:end,:)),1)];

Zprior = [];
for i = p+1:tau+p
    ztemp = eye(M);
    for j = 1:p;
        xlag = Y_prior(i-j,1:M);
        xtemp = zeros(M,M*M);
        for jj = 1:M;
            xtemp(jj,(jj-1)*M+1:jj*M) = xlag;
        end
        ztemp = [ztemp   xtemp];
     end
    Zprior = [Zprior ; ztemp];
end

sse2 = zeros(M,M);
for i = 1:tau
    zhat1 = Zprior((i-1)*M+1:i*M,:);
    sse2 = sse2 + (Y_prior(i+p,:)' - zhat1*B_OLS)*(Y_prior(i+p,:)' - zhat1*B_OLS)';
end

VB_OLS = zeros(K,K);
for i = 1:tau
    zhat1 = Zprior((i-1)*M+1:i*M,:);
    VB_OLS = VB_OLS + zhat1'*inv(sse2)*zhat1;
end
VB_OLS = inv(VB_OLS);

if isstruct(reduced_form)
   for i = 1:M
       err(:,i)=reduced_form(i).resid;
   end
   omega1 = err'*err/rows(err); %Var-Cov matrix of VAR residuals
else
   error('YOU MUST INPUT A VAR RESULTS STRUCTURE')
end

if red_dim==1
    VTheta_OLS = zeros(dim_Theta,dim_Theta);
    xhy = zeros(dim_Theta,1);
    for i = 1:tau
        zhat1 = Zprior((i-1)*M+1:i*M,:)*Xi;
        VTheta_OLS = VTheta_OLS + zhat1'*zhat1;
        xhy = xhy + zhat1'*(Y_prior(p+i,:)');
    end

    VTheta_OLS = pinv(VTheta_OLS);
    Theta_OLS = VTheta_OLS*xhy;

    sse2 = zeros(M,M);
    for i = 1:size(Y_prior,1)-p
        zhat1 = Zprior((i-1)*M+1:i*M,:)*Xi;
        sse2 = sse2 + (Y_prior(i+p,:)' - zhat1*Theta_OLS)*(Y_prior(i+p,:)' - zhat1*Theta_OLS)';
    end

    VTheta_OLS = zeros(dim_Theta,dim_Theta);
    for i = 1:size(Y_prior,1)-p
        zhat1 = Zprior((i-1)*M+1:i*M,:)*Xi;
        VTheta_OLS = VTheta_OLS + zhat1'*inv(sse2)*zhat1;
    end
    VTheta_OLS = pinv(VTheta_OLS);        
    
    VTheta_OLS = diag(diag(VTheta_OLS));        

end

theta_0T=zeros(Nbar,numa+M);
theta_finT=zeros(Nbar,numa+M);
LoglT=zeros(Nbar,1);
tic
for k=1:Nbar

% Settings for csminwel usage

 theta0=mvnrnd(zeros(numa,1),eye(numa),1);
 temp=lognrnd(zeros(1,M),ones(1,M));
 theta0=[theta0,temp];

 tol_crit=1e-4;
 max_iter=1e100;

    [grad0 badg0]=numgrad(@loglike_svar_AB,theta0,M,numa,omega1,t_prior,s_A,S_A,0);
    
    H0=1*eye(numa+M);

    [LogL,theta_fin,grad_theta,h_grad,itct,fcount,retcodeh] = csminwel(@loglike_svar_AB,theta0,H0,[],tol_crit,max_iter,M,numa,omega1,t_prior,s_A,S_A,0);
    
    theta_0T(k,:)=theta0;
    theta_finT(k,:)=theta_fin;
    LoglT(k,1)=-LogL;
end
toc

[Logl_fin,ind_max]=max(real((LoglT)));

theta_fin=theta_finT(ind_max,:);
theta0=theta_0T(ind_max,:);

A_OLS=theta_fin(1:numa)';
VA_OLS=diag(abs(A_OLS));    
sigma_OLS=exp(theta_fin(numa+1:end))';

if red_dim==1
    keep tau theta_fin theta0 Logl_fin B_OLS VB_OLS Theta_OLS VTheta_OLS A_OLS VA_OLS sigma_OLS ident_name constant identification tcodey p standardize omega1 red_dim
else
    keep tau theta_fin theta0 Logl_fin B_OLS VB_OLS A_OLS VA_OLS sigma_OLS ident_name constant identification tcodey p standardize omega1 red_dim
end

save([sprintf('%s_%d_%d_%d_%d_%d_%d_%d.mat','prior_ML_final',tau,...
       identification,constant,p,standardize,tcodey(1),red_dim)]);
