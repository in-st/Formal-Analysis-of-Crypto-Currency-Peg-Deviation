% This code produces the training sample prior.
% This code is based on "Estimating overidentified, non-recursive, time varying coefficients structural VARs"... 
% ...by Fabio Canova and Fernando J. Pérez Forero (Quantitative Economics, Volume 6, Issue 2 (July 2015)).

% Modified to run on GNU Octave 5.2.0 by
% Huachen Li

% -------------------------------------------------------------------------
pkg load statistics
clear all;
clc;
seed=101010101010101;
rand('twister',seed); 
randn('state',seed); 

%% 1. Estimation Setup
%====================================

load Data  %1:Peg Deviation (%); 2:dln(BTC); 3:dln(USDT volume); 4:Fees Per Transaction (USD)

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

t2 = size(ydata,1);
scales=std(ydata,1);
standardize=1;

switch(standardize)
    case 1 % Standardized
        ydata = (ydata- repmat(mean(ydata,1),t2,1))./repmat(std(ydata,1),t2,1);
    case 2 % Non-standardized
        ydata(:,[3,6])=0.01*ydata(:,[3,6]);
end


Y=ydata(:,[4 2 1 3]); % Variable ordering

% Number of observations and variables in Y
t=size(Y,1); % t is the time-series observations of Y
M=size(Y,2); % M is the dimensionality of Y


tau = 90; % tau is the size of the training sample
identification=1; % Cholesky
ident_A    = [2*triu(ones(M,M),+1)+eye(M)]';
ident_name = {'Cholesky'};
p = 2; % Lag order
constant=1; % 0. None. 1: Intercept. 2: Trend+intercept.
red_dim=0; % Reduced-dimension as in Canova and Ciccarelli (2009)
Nbar=40; % Number of initial points

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

ylag = mlag2(Y,p); 
ylag = ylag(p+tau+1:t,:);
K = M*(constant) + p*(M^2); 
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

y = Y(tau+p+1:t,:)';
yearlab = yearlab(tau+p+1:t);
t=size(y,2);  

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

%% 2. Maximum Likelihood Estimation
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
LoglT1 = -LoglT;
A_OLS=theta_fin(1:numa)';
VA_OLS=diag(abs(A_OLS));    
sigma_OLS=exp(theta_fin(numa+1:end))';

if red_dim==1
    keep tau theta_fin theta0 Logl_fin B_OLS VB_OLS Theta_OLS VTheta_OLS A_OLS VA_OLS sigma_OLS ident_name constant identification tcodey p standardize omega1 red_dim
else
    keep LoglT1 tau theta_fin theta0 Logl_fin B_OLS VB_OLS A_OLS VA_OLS sigma_OLS ident_name constant identification tcodey p standardize omega1 red_dim
end

save([sprintf('%s_%d_%d_%d_%d_%d_%d_%d.mat','prior_ML_final',tau,...
       identification,constant,p,standardize,tcodey(1),red_dim)]);
