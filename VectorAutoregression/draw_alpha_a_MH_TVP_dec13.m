%-------------------------------------------------------------------------------------------
%   STEP II: Draw A(t) from p(At|y,B,Sigma,V)
%-------------------------------------------------------------------------------------------

% Substract from the data y(t), the mean Z x B(t)
yhat = zeros(M,t);
for i = 1:t
    yhat(:,i) = y(:,i) - Z((i-1)*M+1:i*M,:)*Btdraw(:,i);
end
    
ytildeA=zeros(M,t);
Stemp = zeros(M,numa);
Zc = zeros(t*M,numa);
siga2temp=zeros(M,M,t);

for i = 1:t    
    ytildeA(:,i)=reshape(s_A,M,M)*yhat(:,i);
    Stemp=-kron(yhat(:,i)',eye(M))*S_A;
    Zc((i-1)*M+1:i*M,:) =  Stemp;
    sigatemp = sigt((i-1)*M+1:i*M,:);
    siga2temp(:,:,i) = sigatemp*sigatemp';        
end

%%  Extended Kalman Filter (EKF)

% Non linear State space model (see Harvey (1989) pp. 160 
% y(t)      =   Z_t(alpha(t))   +   F_t(alpha(t)) x epsilon(t)  
% alpha(t)  =   T_t(alpha(t-1)) +   R_t(alpha(t-1)) x eta(t)
%
% epsilon(t)~   N(0,Q)
% eta(t)    ~   N(0,S)

    alpha0_a_guess=A_OLS;
    alpha_a_ekf=zeros(numa,t);    
    Pt_t_1_a=zeros(numa,numa,t);
    P_ekf_a=zeros(numa,numa,t);
    P_ekf_sm_a=zeros(numa,numa,t);
    P_ekf_sm_star_a=zeros(numa,numa,t);
    
    Z1 = Zc(1:M,:);        % From measurement equation
    P0_a=A_0_prvar;       % Initial value of var(alpha) (unconditional = 1)    
    Q = siga2temp(:,:,1);% Variance of measurement equation error     
    if red_dim==1 && Xi_rand==0
        Q=reshape(S_A*alpha0_a_guess+s_A,M,M)*Q*reshape(S_A*alpha0_a_guess+s_A,M,M)';
    end
    T_t = eye(numa);        % Jacobian
    Pt_t_1_a(:,:,1) = T_t*P0_a*T_t'+Sadraw;    
    K1 = Pt_t_1_a(:,:,1)'*Z1'/(Z1*Pt_t_1_a(:,:,1)*Z1'+Q); % Kalman gain
    P_ekf_a(:,:,1) = Pt_t_1_a(:,:,1)-Pt_t_1_a(:,:,1)*Z1'*inv(Z1*Pt_t_1_a(:,:,1)*Z1'+Q)*Z1*Pt_t_1_a(:,:,1);
    alpha_a_ekf(:,1) = alpha0_a_guess + K1*(ytildeA(:,1)-Z1*alpha0_a_guess);

    % Filtering
    %----------
    for i=2:t
        % Prediction:        
        Pt_t_1_a(:,:,i) = T_t*P_ekf_a(:,:,i-1)*T_t'+Sadraw;    
        alphat_a_t_1    = T_t*alpha_a_ekf(:,i-1);        
        % Updating:
        Q = siga2temp(:,:,i); % Variance of measurement equation error
         if red_dim==1 && Xi_rand==0
            Q=reshape(S_A*alphat_a_t_1+s_A,M,M)*Q*reshape(S_A*alphat_a_t_1+s_A,M,M)';
         end
        Z1 = Zc((i-1)*M+1:i*M,:);
        K1 = Pt_t_1_a(:,:,i)'*Z1'/(Z1*Pt_t_1_a(:,:,i)*Z1'+Q);
        P_ekf_a(:,:,i)   = Pt_t_1_a(:,:,i)-Pt_t_1_a(:,:,i)*Z1'*inv(Z1*Pt_t_1_a(:,:,i)*Z1'+Q)*Z1*Pt_t_1_a(:,:,i);
        alpha_a_ekf(:,i) = alphat_a_t_1 + K1*(ytildeA(:,i)-Z1*alphat_a_t_1);    
    end

    % Smoothing
    %----------
    alpha_a_ekf_sm=zeros(numa,t);
    alpha_a_ekf_sm(:,t)=alpha_a_ekf(:,t);    
    P_ekf_sm_a(:,:,t)=P_ekf_a(:,:,t)-P_ekf_a(:,:,t)*T_t'*inv(T_t*P_ekf_a(:,:,t)*T_t'+Sadraw)*T_t*P_ekf_a(:,:,t);

    for i=t-1:-1:1
        alpha_a_ekf_sm(:,i)=alpha_a_ekf(:,i)+P_ekf_a(:,:,i)*T_t'*inv(Pt_t_1_a(:,:,i+1))*(alpha_a_ekf_sm(:,i+1)-T_t'*alpha_a_ekf(:,i));        
        P_ekf_sm_a(:,:,i)=P_ekf_a(:,:,i)-P_ekf_a(:,:,i)*T_t'*inv(Pt_t_1_a(:,:,i+1))*T_t*P_ekf_a(:,:,i);
    end
    
% Set proposal density variance
% -----------------------------
    for i=1:t
        P_ekf_sm_star_a(:,:,i)=P_ekf_sm_a(:,:,i);
    end
    
%%  Metropolis Hastings sampling starts here

%--------------------------------------------------------------------
% Draw alpha_a(T)
%--------------------------------------------------------------------

% Draw a candidate
alpha_a_can=zeros(numa,t);

for i=1:t
    % Generate a random draw from a Multivariate T-Student
    
    % (1) Multivariate normal draw
    P_temp=norm_rnd(c_a*P_ekf_sm_star_a(:,:,i)); 
    
    % (2) Chi-squared draw
    u_temp=chi2rnd(df);  
    
    % Combine (1) and (2)
    alpha_a_can(:,i) = Atdraw(:,i)+sqrt(df/u_temp)*P_temp;
end

% Reset densities
l_can_a=0;
l_draw_a=0;

for i=1:t   
    if i<t
    % Evaluate candidate log-density
         alpha_err=ytildeA(:,i)-Zc((i-1)*M+1:i*M,:)*alpha_a_can(:,i);
                
         l_can_a=l_can_a+0.5*M*log(2*pi)+log(det(reshape(S_A*alpha_a_can(:,i)+s_A,M,M)))...
                      -0.5*alpha_err'*inv(siga2temp(:,:,i))*alpha_err...
                   +log(mvnpdf(alpha_a_can(:,i+1),alpha_a_can(:,i),Sadraw));

    % Evaluate current state log-density                  
         alpha_err=ytildeA(:,i)-Zc((i-1)*M+1:i*M,:)*Atdraw(:,i);
         
         l_draw_a=l_draw_a+0.5*M*log(2*pi)+log(det(reshape(S_A*Atdraw(:,i)+s_A,M,M)))...
                  -0.5*alpha_err'*inv(siga2temp(:,:,i))*alpha_err+...
                   +log(mvnpdf(Atdraw(:,i+1),Atdraw(:,i),Sadraw));
    else
    % Evaluate candidate log-density
         alpha_err=ytildeA(:,i)-Zc((i-1)*M+1:i*M,:)*alpha_a_can(:,i);
                
         l_can_a=l_can_a+0.5*M*log(2*pi)+log(det(reshape(S_A*alpha_a_can(:,i)+s_A,M,M)))...
                      -0.5*alpha_err'*inv(siga2temp(:,:,i))*alpha_err;

    % Evaluate current state log-density                  
         alpha_err=ytildeA(:,i)-Zc((i-1)*M+1:i*M,:)*Atdraw(:,i);
         
         l_draw_a=l_draw_a+0.5*M*log(2*pi)+log(det(reshape(S_A*Atdraw(:,i)+s_A,M,M)))...
                  -0.5*alpha_err'*inv(siga2temp(:,:,i))*alpha_err;

    end
end

if imag(l_draw_a)==0
else
   l_draw_a=-Inf;
end

if imag(l_can_a)==0
else
   l_can_a=-Inf;
end

% Accept candidate with probability alpha
if l_can_a==-Inf
    l_alpha_a =-Inf;
else
   l_alpha_a = min(l_can_a-l_draw_a,0);
end

% Check bounds
if min(min(alpha_a_can > LB))==1
   if min(min(alpha_a_can < UB))==1
   else
      l_alpha_a=-9*1e200;
      out_counter=out_counter+1;
   end
else
    l_alpha_a=-9*1e200;
    out_counter=out_counter+1;
end

if log(rand)<l_alpha_a
   Atdraw=alpha_a_can;
   acc_counter_a=acc_counter_a+1;
end

 %--------------
 % Draw S_a
 %--------------
 
    %=====| Draw S, the covariance of A(t) (from iWishart)
    % Take the SSE in the state equation of A(t)
    alpha_a_temp = Atdraw(:,2:t)' - Atdraw(:,1:t-1)';
    sse_2 = zeros(numa,numa);
    for i = 1:t-1
        sse_2 = sse_2 + alpha_a_temp(i,:)'*alpha_a_temp(i,:);
    end
    % ...and subsequently draw S, the covariance matrix of A(t) 
   
    Sainv = inv(sse_2 + Sa_prmean);
    Sainvdraw = wish(Sainv,t+Sa_prvar);
    Sadraw = inv(Sainvdraw); % this is a draw from S  
