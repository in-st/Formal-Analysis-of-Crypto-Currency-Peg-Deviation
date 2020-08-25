%------------------------------------------------------------------------------------------
%   STEP III: Draw diagonal VAR covariance matrix log-SIGMA(t)
%------------------------------------------------------------------------------------------

    % Recover matrix A(t)
    
    capAt = zeros(M*t,M);
    for i = 1:t
        capAt((i-1)*M+1:i*M,:) = reshape(S_A*Atdraw(:,i)+s_A,M,M);
    end

    yhat2 = zeros(M,t);
    for i = 1:t
        yhat2(:,i) = capAt((i-1)*M+1:i*M,:)*yhat(:,i);
    end
    

%%  Extended Kalman Filter (EKF)

% Non linear State space model (see Harvey (1989) pp. 160 
% y(t)      =   Z_t(alpha(t))   +   F_t(alpha(t)) x epsilon(t)  
% alpha(t)  =   T_t(alpha(t-1)) +   R_t(alpha(t-1)) x eta(t)
%
% epsilon(t)~   N(0,Q)
% eta(t)    ~   N(0,S)

for m=1:M
    sigma0_a_guess=sigma_OLS(m,1);
    sigma_a_ekf=zeros(1,t);
    Pt_t_1_s=zeros(1,t);
    P_ekf_s=zeros(1,t);
    P_ekf_sm_s=zeros(1,t);
    P_ekf_sm_star_s=zeros(1,t);
    
    Z1 = 0;        % From measurement equation
    P0_s=sigma_prvar(m,m);       % Initial value of var(sigma)
    T_t = delta*sigma0_a_guess*(1-delta+delta*sigma0_a_guess^2+delta*0)^(-0.5);        % Jacobian
    Pt_t_1_s(1,1) = T_t*P0_s*T_t'+Wdraw(m,m);  
    sigmat_a_t_1    = sqrt((1-delta)+sigma0_a_guess^2+delta*0);
    Q = sigmat_a_t_1*sigmat_a_t_1';% Variance of measurement equation error     
    K1 = Pt_t_1_s(1,1)'*Z1'/(Z1*Pt_t_1_s(1,1)*Z1'+Q); % Kalman gain
    P_ekf_s(1,1) = Pt_t_1_s(1,1)-Pt_t_1_s(1,1)*Z1'*inv(Z1*Pt_t_1_s(1,1)*Z1'+Q)*Z1*Pt_t_1_s(1,1);
    y_for=yhat2(m,1)-Z1*sigmat_a_t_1;
    sigma_a_ekf(1,1) = sigmat_a_t_1 + K1*(y_for);

    % Filtering
    %----------
    for i=2:t
        % Prediction:        
        Pt_t_1_s(1,i) = T_t*P_ekf_s(1,i-1)*T_t'+Wdraw(m,m);    
        sigmat_a_t_1    = sqrt((1-delta)+delta*sigma_a_ekf(1,i-1)^2+delta*y_for^2);        
        % Updating:
        Q = sigmat_a_t_1*sigmat_a_t_1';% Variance of measurement equation error     
        Z1 = 0;
        K1 = Pt_t_1_a(1,i)'*Z1'/(Z1*Pt_t_1_a(1,i)*Z1'+Q);
        P_ekf_s(1,i)   = Pt_t_1_s(1,i)-Pt_t_1_s(1,i)*Z1'*inv(Z1*Pt_t_1_s(1,i)*Z1'+Q)*Z1*Pt_t_1_s(1,i);
        y_for=yhat2(1,i)-Z1*sigmat_a_t_1;
        sigma_a_ekf(1,i) = sigmat_a_t_1 + K1*(y_for);    
        T_t = delta*sigma_a_ekf(1,i)*(1-delta+delta*sigma_a_ekf(1,i)^2+delta*y_for^2)^(-0.5);        % Jacobian
    end

    % Smoothing
    %----------
    sigma_a_ekf_sm=zeros(1,t);
    sigma_a_ekf_sm(1,t)=sigma_a_ekf(1,t);    
    y_for=yhat2(1,t)-Z1*sigma_a_ekf(1,t);
    T_t = delta*sigma_a_ekf_sm(1,t)*(1-delta+delta*sigma_a_ekf_sm(1,t)^2+delta*y_for^2)^(-0.5);        % Jacobian
    P_ekf_sm_s(1,t)=P_ekf_s(1,t)-P_ekf_s(1,t)*T_t'*inv(T_t*P_ekf_s(1,t)*T_t'+Wdraw(m,m))*T_t*P_ekf_s(1,t);

    for i=t-1:-1:1
        Z1 = 0;
        y_for=yhat2(1,i)-Z1*sigma_a_ekf(1,i);
        sigma_a_ekf_sm(1,i)=sigma_a_ekf(1,i)+P_ekf_s(1,i)*T_t'*inv(Pt_t_1_s(1,i+1))*(sigma_a_ekf_sm(1,i+1)-sqrt((1-delta)+delta*sigma_a_ekf(1,i)^2+delta*y_for^2));
        T_t = delta*sigma_a_ekf_sm(1,i)*(1-delta+delta*sigma_a_ekf_sm(1,i)^2+delta*y_for^2)^(-0.5);        % Jacobian
        P_ekf_sm_s(1,i)=P_ekf_s(1,i)-P_ekf_s(1,i)*T_t'*inv(T_t*P_ekf_s(1,i)*T_t'+Wdraw(m,m))*T_t*P_ekf_s(1,i);
    end
    
% Set proposal density variance
% -----------------------------
    for i=1:t
        P_ekf_sm_star_s(1,i)=P_ekf_sm_s(1,i);
    end
    
%  Metropolis Hastings sampling starts here

%--------------------------------------------------------------------
% Draw sigma_a(T)
%--------------------------------------------------------------------

% Draw a candidate
sigma_a_can=zeros(1,t);

for i=1:t
    % Generate a random draw from a Multivariate T-Student
    
    % (1) Multivariate normal draw
    P_temp=norm_rnd(c_s*P_ekf_sm_star_s(1,i)); 
    
    % (2) Chi-squared draw
    u_temp=chi2rnd(df);  
    
    % Combine (1) and (2)
    sigma_a_can(1,i) = Sigtdraw(m,i)+sqrt(df/u_temp)*P_temp;
end

% Reset densities
l_can_s=0;
l_draw_s=0;

for i=1:t   
    if i<t
    % Evaluate candidate log-density
         sigma_err=yhat2(m,i)-0;
                
         l_can_s=l_can_s+0.5*log(2*pi)+...
                      -0.5*sigma_err'*inv(sigma_a_can(1,i)^2)*sigma_err...
                   +log(mvnpdf(sigma_a_can(1,i+1),sqrt(sigma_a_can(1,i)^2+yhat2(m,i)^2),Wdraw(m,m)));

    % Evaluate current state log-density                  
         sigma_err=yhat2(m,i)-0;
         
         l_draw_s=l_draw_s+0.5*log(2*pi)+...
                      -0.5*sigma_err'*inv(Sigtdraw(m,i)^2)*sigma_err...
                   +log(mvnpdf(Sigtdraw(m,i+1),sqrt(Sigtdraw(m,i)^2+yhat2(m,i)^2),Wdraw(m,m)));
    else
        
    % Evaluate candidate log-density
         sigma_err=yhat2(m,i)-0;                
         l_can_s=l_can_s+0.5*log(2*pi)-0.5*sigma_err'*inv(sigma_a_can(1,i)^2)*sigma_err;

    % Evaluate current state log-density                  
         sigma_err=yhat2(m,i)-0;         
         l_draw_s=l_draw_s+0.5*log(2*pi)-0.5*sigma_err'*inv(Sigtdraw(m,i)^2)*sigma_err;

    end
end

if imag(l_draw_s)==0
else
   l_draw_s=-Inf;
end

if imag(l_can_s)==0
else
   l_can_s=-Inf;
end

% Accept candidate with probability alpha
if l_can_s==-Inf
    l_sigma_a =-Inf;
else
   l_sigma_a = min(l_can_s-l_draw_s,0);
end

% Check bounds
if min(min(sigma_a_can > LB_s(m,:)))==1
   if min(min(sigma_a_can < UB_s(m,:)))==1
   else
      l_alpha_a=-9*1e200;
      out_counter_s(m,1)=out_counter_s(m,1)+1;
   end
else
    l_sigma_a=-9*1e200;
    out_counter_s(m,1)=out_counter_s(m,1)+1;
end

if log(rand)<l_sigma_a
   Sigtdraw(m,:)=sigma_a_can;
   acc_counter_s(m,1)=acc_counter_s(m,1)+1;
end
    
end
    
    %%
    
    % Draws in Sigtdraw are in logarithmic scale (log-volatilies). Create 
    % original standard deviations of the VAR covariance matrix
    sigtemp = eye(M);
    sigt = zeros(M*t,M);
    
    for i = 1:t
         if red_dim==1 && Xi_rand==0
            at=chol(reshape(S_A*Atdraw(:,i)+s_A,M,M)*Z((i-1)*M+1:i*M,:)*Z((i-1)*M+1:i*M,:)'*reshape(S_A*Atdraw(:,i)+s_A,M,M)')';            
             sigtemp = at*diag((Sigtdraw(:,i)));
             sigt((i-1)*M+1:i*M,:) = sigtemp;
         else
            for j = 1:M
                sigtemp(j,j) = Sigtdraw(j,i);
            end
            sigt((i-1)*M+1:i*M,:) = sigtemp;
        end
    end

    %=====| Draw W, the covariance of SIGMA(t) (from iWishart)
    % Get first differences of Sigtdraw to compute the SSE

     Wdraw=zeros(M,M);
    
     for k=1:M    
         Sigttemp = Sigtdraw(k,2:t).^2' - (1-delta)*ones(t-1,1)- delta*Sigtdraw(k,1:t-1).^2'- delta*yhat2(k,1:t-1).^2';
         sse_2 = 0;
         for i = 1:t-1
             sse_2 = sse_2 + Sigttemp(i,:)'*Sigttemp(i,:);
         end
         Winv = inv(sse_2 + W_prmean(k,k));
         Winvdraw = wish(Winv,t+W_prvar);
         Wdraw(k,k) = inv(Winvdraw);  % this is a draw from W    
     end
        
