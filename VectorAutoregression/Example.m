% A general algorithm for estimating SVARs

% Example: Static SVAR

% Fernando Pérez Forero (fernandojose.perez@upf.edu)
% This version: May 5th, 2013

%% 0. Data simulation

clear all
clc
seed=0; rand( 'state' , seed ); randn('state',seed ); % set the random seed

% Identification

ident_A=[1 0 2;
         2 1 0;
         0 2 1];

alpha_true=[0.8 0.5 0.5]';
     
M=size(ident_A,1);    
     
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

% Initial value

t=500;

A_alpha=reshape(S_A*alpha_true+s_A,M,M);

y=zeros(M,t);

for i=1:t
    y(:,i)=mvnrnd(zeros(M,1),inv(A_alpha)*inv(A_alpha)',1)';
end

%% 1. Re-parametrization

ytildeA=zeros(M,t);
Z = zeros(t*M,numa);

for i = 1:t
    
    ytildeA(:,i)=reshape(s_A,M,M)*y(:,i);

    Stemp=-kron(y(:,i)',eye(M))*S_A;

    Z((i-1)*M+1:i*M,:) =  Stemp;
end

clear Stemp

ZZ=zeros(numa,numa);
Zy=zeros(numa,1);

for i = 1:t    
    ZZ=ZZ+Z((i-1)*M+1:i*M,:)'*Z((i-1)*M+1:i*M,:);
    Zy=Zy+Z((i-1)*M+1:i*M,:)'*ytildeA(:,i);    
end

alpha_OLS=ZZ\Zy;

sse_OLS=zeros(M,M);

for i=1:t
    err_alpha=ytildeA(:,i)-Z((i-1)*M+1:i*M,:)*alpha_OLS;    
    sse_OLS=sse_OLS+(err_alpha*err_alpha');
end

P_OLS=zeros(numa,numa);
for i=1:t
    P_OLS=P_OLS+Z((i-1)*M+1:i*M,:)'*inv(sse_OLS)*Z((i-1)*M+1:i*M,:);
end

P_OLS=inv(P_OLS);

clear ZZ Zy 

%% 2. Metropolis-Hastings 

test=2;

switch(test)
    case 1
        nrep     =  10000;  % Number of replications
        nburn    =   5000;  % Number of burn-in-draws
        it_print =   1000;  % Print in the screen every "it_print"-th iteration
        nthin    =     20;  % Thinning factor
    case 2
        nrep     =  50000;  % Number of replications
        nburn    = 100000;  % Number of burn-in-draws
        it_print =   1000;  % Print in the screen every "it_print"-th iteration
        nthin    =    100;  % Thinning factor          
    case 3
        nrep     =  50000;  % Number of replications
        nburn    = 450000;  % Number of burn-in-draws
        it_print =   1000;  % Print in the screen every "it_print"-th iteration
        nthin    =    100;  % Thinning factor                  
end

% Pre-allocating memory
% ---------------------
alpha_post=zeros(numa,nrep/nthin);

acc_counter=0;
Bound=20;
LB =-Bound*ones(numa,1);
UB = Bound*ones(numa,1);
out_counter=0;

prop=2; % 1. Normal, 2. T-Student

% Scale of the proposal density
switch(prop)
    case 1 % Normal
        df=0;
        %r=1;        
        %r=2*1e-2; % 8%
        %r=1e-2; % 17%
        r=0.5*1e-2; % 28%
        %r=1e-3; % 75%
        %r=1e-4; % 91%
    case 2 % T-student
        df=4;
        %df=5;
        r=0.5*1e-2; % 24%
end

alpha_draw=alpha_OLS;
P_draw=P_OLS;
sse_draw=sse_OLS;

tic;
disp('Number of iterations');

rep=0;
lpostdraw = -9e+200;

% Main loop
for irep = 1:nrep + nburn
    
    % Print iterations
    if mod(irep,it_print) == 0
        disp(irep);toc;
    end
    
    % Draw from proposal density 
    switch(prop)
        case 1 % alpha_{i+1} ~ N(alpha_{i},r*P(alpha_{i}))
            alpha_can = mvnrnd(alpha_draw,r*P_draw,1)';
        case 2 % alpha_{i+1} ~ T(alpha_{i},r*P(alpha_{i}))
            alpha_can = alpha_draw+sqrt(r)*mvtrnd(P_draw,df,1)';
    end

    if min(alpha_can> LB)==1
        if min(alpha_can < UB)==1
            
            lpostcan=0;
            
            for i=1:t
                alpha_err=ytildeA(:,i)-Z((i-1)*M+1:i*M,:)*alpha_can;
                lpostcan=lpostcan-0.5*alpha_err'*alpha_err;
            end
           
            % Evaluate the jacobian term 
            lpostcan=lpostcan+t*log(det(reshape(S_A*alpha_can+s_A,M,M)))...
                            -0.5*M*t*log(2*pi);            

            % Acceptance probability
            laccprob = lpostcan-lpostdraw;
            
        else
            laccprob=-9e+200;
            out_counter=out_counter+1;
        end
    else
        laccprob=-9e+200;
        out_counter=out_counter+1;
    end

    %Accept candidate draw with log prob = laccprob, else keep old draw
    if log(rand)<laccprob
        lpostdraw=lpostcan;
        alpha_draw=alpha_can;
        acc_counter=acc_counter+1;     
    end
    
    sse_draw=zeros(M,M);

    for i=1:t
        err_alpha=ytildeA(:,i)-Z((i-1)*M+1:i*M,:)*alpha_draw;
        sse_draw=sse_draw+(err_alpha*err_alpha');
    end
    
    P_draw=zeros(numa,numa);
    for i=1:t
        P_draw=P_draw+Z((i-1)*M+1:i*M,:)'*inv(sse_draw)*Z((i-1)*M+1:i*M,:);        
    end

    P_draw=inv(P_draw);
            
    if irep > nburn && mod((irep-nburn),nthin)==0 
        alpha_post(:,(irep-nburn)/nthin)=alpha_draw;    
    end       
end

sprintf('Average Acceptance rate in Metropolis Step sampling (A^T): %s percent',num2str(mean(100*acc_counter/(nrep + nburn))))

sprintf('Rate of out-of-bounds draws in Metropolis Step sampling (A^T): %s percent',num2str(mean(100*out_counter/(nrep + nburn))))

save([sprintf('%s_%d_%d_%d_%d_%d_%d_%d_%d_%d.mat','Example_constant',nrep+nburn,nburn,nthin,r,seed,alpha_true(1),t,prop,df)]);


%% Figures

    % Graph Multiple histograms
    
    LW=2;
    FS=16;

    temp1=alpha_post;
    
    bins=20;
    numa=size(temp1,1);
    
    figure
    set(gcf,'Color',[1 1 1])
    set(gcf,'defaultaxesfontsize',FS)    
    for i=1:numa
        subplot(1,numa,i)
        edges=linspace(min(temp1(i,:)),max(temp1(i,:)),bins);
        [n1, xout1]=hist(temp1(i,:),edges);
        plot(edges,n1,'-b','LineWidth',LW);
        hold on
        temp=alpha_true(i)*ones(1,length(xout1));
        [n3, xout3]=hist(temp,edges);
        bar(edges,n3);
        title(sprintf('\\alpha_{%d}',i),'FontSize',FS)    
    end
        AX=legend('Posterior','True');
        LEG = findobj(AX,'type','text');
        set(LEG,'FontSize',FS)
        
    clear temp1 temp2 temp
    
convcheck(alpha_post);

