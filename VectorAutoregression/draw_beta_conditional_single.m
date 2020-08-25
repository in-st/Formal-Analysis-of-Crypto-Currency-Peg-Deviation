% -----------------------------------------------------------------------------------------
%   STEP I: Sample B from p(B|y,A,Sigma,V)
%   Here we implement Koop & Potter (2011)'s Single-Move algorithm
% -----------------------------------------------------------------------------------------

%% Draw a candidate

Btdrawc=zeros(K,t);
Ind_beta=zeros(1,t);

Rcan=zeros(1,t); % This is the integrating constant
Rdraw=zeros(1,t); % This is the integrating constant

mu_t=zeros(K,t);
Omega_t=zeros(K,K,t);

for i=1:t
    R = Ht((i-1)*M+1:i*M,:);
    H=Z((i-1)*M+1:i*M,:);
    if i==1
        G_t=0.5*Qdraw*H'/(H*Qdraw*H'+R);
        Omega_t(:,:,i)=0.5*(eye(K)-G_t*H)*Qdraw;        
        mu_t(:,i)=0.5*(B_0_prmean+Btdraw(:,i+1))+...
            G_t*(y(:,i)-H*0.5*(B_0_prmean+Btdraw(:,i+1)));          
    elseif i<t                
        G_t=0.5*Qdraw*H'/(H*Qdraw*H'+R);
        Omega_t(:,:,i)=0.5*(eye(K)-G_t*H)*Qdraw;
        mu_t(:,i)=0.5*(Btdraw(:,i-1)+Btdraw(:,i+1))+...
            G_t*(y(:,i)-H*0.5*(Btdraw(:,i-1)+Btdraw(:,i+1)));        
    else
        G_t=Qdraw*H'/(H*Qdraw*H'+R);
        Omega_t(:,:,i)=(eye(K)-G_t*H)*Qdraw;
        mu_t(:,i)=Btdraw(:,i-1)+G_t*(y(:,i)-H*(Btdraw(:,i-1)));        
    end    
    Btdrawc(:,i) = mvnrnd(mu_t(:,i),c_b*Omega_t(:,:,i),1)';
    
    % Evaluate stability at period i
    ctemp1 = zeros(M*p,M*p);    
    for j = 1:p-1
        ctemp1(j*M+1:M*(j+1),M*(j-1)+1:j*M) = eye(M);
    end

    splace = 0;
    BBtempor = Btdrawc(M+1:end,i);
    for ii = 1:p
        for iii = 1:M
            ctemp1(iii,(ii-1)*M+1:ii*M) = BBtempor(splace+1:splace+M,1)';
            splace = splace + M;
        end
    end
    
    if max(abs(eig(ctemp1)))<1;
       Ind_beta(1,i)=1;
    end
    
    if i<t
        
        % Evaluate the integrating constant
    
        Ind_betaR=zeros(1,Rsimul);
        Ind_betaRd=zeros(1,Rsimul);
    
        for k=1:Rsimul            
            Bt_tempR=mvnrnd(Btdrawc(:,i),Qdraw,1)';        
            Bt_tempRd=mvnrnd(Btdraw(:,i),Qdraw,1)';        
            
            % Evaluate stability
    
            ctempR = zeros(M*p,M*p);            
            ctempRd = zeros(M*p,M*p);            
            for j = 1:p-1
                ctempR(j*M+1:M*(j+1),M*(j-1)+1:j*M) = eye(M);                
                ctempRd(j*M+1:M*(j+1),M*(j-1)+1:j*M) = eye(M);                
            end

            splace = 0;
            BBtempor = Bt_tempR(M+1:end,:);
            BBtempord = Bt_tempRd(M+1:end,:);
            
            for ii = 1:p
                for iii = 1:M
                    ctempR(iii,(ii-1)*M+1:ii*M) = BBtempor(splace+1:splace+M,1)';
                    ctempRd(iii,(ii-1)*M+1:ii*M) = BBtempord(splace+1:splace+M,1)';
                    splace = splace + M;
                end
            end
    
            if max(abs(eig(ctempR)))<0.9999  % Candidate B(t)
                Ind_betaR(1,k)=1;
            end    
            
            if max(abs(eig(ctempRd)))<0.9999 % Current B(t)
                Ind_betaRd(1,k)=1;
            end                
        end
    
        Rcan(1,i)=sum(Ind_betaR)/Rsimul;
        Rdraw(1,i)=sum(Ind_betaRd)/Rsimul;
    
        % Evaluate acceptance probability

        if Rcan(1,i)==0
           l_draw_b=0;
        else    
           l_draw_b=min(Ind_beta(1,i)*Rdraw(1,i)/Rcan(1,i),1);
        end

    else
        
        % Evaluate acceptance probability
        l_draw_b=Ind_beta(1,i);        

    end
    
    % Accept candidate with probability alpha
    
    if rand<l_draw_b
        Btdraw(:,i)=Btdrawc(:,i);
        Beta_counter(1,i)=Beta_counter(1,i)+1; 
    end    
end

%% Draw Q, the covariance of B(t) (from iWishart)

    % Take the SSE in the state equation of B(t)
    Btemp = Btdraw(:,2:t)' - Btdraw(:,1:t-1)';
    sse_2 = zeros(K,K);
    for i = 1:t-1
        sse_2 = sse_2 + Btemp(i,:)'*Btemp(i,:);
    end
    
    % Draw a candidate of Q, the covariance matrix of B(t)
    Qinv = inv(sse_2 + Q_prmean);
    Qinvdraw = wish(Qinv,t+Q_prvar);
    Qdrawc = inv(Qinvdraw);  
    
    RcanQ=zeros(1,t); % This is the integrating constant
    RdrawQ=zeros(1,t); % This is the integrating constant     
    
    for i=1:t
                
        % Evaluate the integrating constant
    
        Ind_QR=zeros(1,Rsimul);
        Ind_QRd=zeros(1,Rsimul);
    
        for k=1:Rsimul            
            
            if i>1            
                Bt_tempR=mvnrnd(Btdraw(:,i-1),Qdrawc,1)';        
                Bt_tempRd=mvnrnd(Btdraw(:,i-1),Qdraw,1)';                    
            else
                Bt_tempR=mvnrnd(B_0_prmean,Qdrawc,1)';        
                Bt_tempRd=mvnrnd(B_0_prmean,Qdraw,1)';                        
            end
            
            % Evaluate stability
    
            ctempR = zeros(M*p,M*p);            
            ctempRd = zeros(M*p,M*p);            
                
            for j = 1:p-1
                ctempR(j*M+1:M*(j+1),M*(j-1)+1:j*M) = eye(M);                
                ctempRd(j*M+1:M*(j+1),M*(j-1)+1:j*M) = eye(M);                
            end

            splace = 0;
            BBtempor = Bt_tempR(M+1:end,:);
            BBtempord = Bt_tempRd(M+1:end,:);
            
            for ii = 1:p
                for iii = 1:M
                    ctempR(iii,(ii-1)*M+1:ii*M) = BBtempor(splace+1:splace+M,1)';
                    ctempRd(iii,(ii-1)*M+1:ii*M) = BBtempord(splace+1:splace+M,1)';
                    splace = splace + M;
                end
            end
    
            if max(abs(eig(ctempR)))<0.9999
               Ind_QR(1,k)=1;
            end    
            
            if max(abs(eig(ctempRd)))<0.9999
               Ind_QRd(1,k)=1;
            end                
        end
    
        RcanQ(1,i)=sum(Ind_QR)/Rsimul;
        RdrawQ(1,i)=sum(Ind_QRd)/Rsimul;
    end
    
    % Evaluate acceptance probability
        
    lRcanQprod=0;
    lRdrawQprod=0;
    
    for i=1:t
        lRcanQprod=lRcanQprod+log(RcanQ(1,i));
        lRdrawQprod=lRdrawQprod+log(RdrawQ(1,i));
    end
    
    if lRcanQprod==-Inf
        l_draw_Q=0;
    else
        l_draw_Q=min(lRdrawQprod-lRcanQprod,0);
    end
    
    % Accept candidate with probability alpha
    if log(rand)<l_draw_Q
       Qdraw=Qdrawc;
       Q_counter=Q_counter+1; 
    end    
    