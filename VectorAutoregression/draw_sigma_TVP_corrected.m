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

    % yhat is the vector y(t) - Z x B(t) defined previously. Multiply yhat
    % with capAt, i.e the lower triangular matrix A(t). Then take squares
    % of the resulting quantity (saved in matrix y2)
    y2 = [];
    for i = 1:t
        if red_dim==1 && Xi_rand==0
            atinv=chol(reshape(S_A*Atdraw(:,i)+s_A,M,M)*Z((i-1)*M+1:i*M,:)*Z((i-1)*M+1:i*M,:)'*reshape(S_A*Atdraw(:,i)+s_A,M,M)')';            
            ytemps = atinv\yhat2(:,i);            
        else
            ytemps = yhat2(:,i);
        end
        
        y2 = [y2  (ytemps.^2)];        
    end   

    % Transform to the log-scale but also add the 'offset constant' to prevent
    % the case where y2 is zero (in this case, its log would be -Infinity) 
    yss = zeros(t,M);
    for i = 1:M
        yss(:,i) = log(y2(i,:)' + 1e-3);
    end
    
    % Draw s^T
    % --------
    % Draw statedraw (chi square approximation mixture component) conditional on previous Sigtdraw
    % This is used BEFORE sampling the log-volatilities Sigtdraw: Del Negro & Primiceri (2013)
    for jj = 1:M
        for i = 1:t
            for j = 1:numel(m_s)
                temp1= (1/sqrt(2*pi*u2_s(j)))*exp(-.5*(((yss(i,jj) - Sigtdraw(jj,i) - m_s(j) + 1.2704)^2)/u2_s(j)));
                prw(j,1) = q_s(j,1)*temp1;
            end
            prw = prw./sum(prw);
            cprw = cumsum(prw);
            trand = rand(1,1);
            if trand < cprw(1,1); imix=1;
            elseif trand < cprw(2,1), imix=2;
            elseif trand < cprw(3,1), imix=3;
            elseif trand < cprw(4,1), imix=4;
            elseif trand < cprw(5,1), imix=5;
            elseif trand < cprw(6,1), imix=6;
            else imix=7; 
            end
            statedraw(i,jj)=imix;  % this is a draw of the mixture component index
        end
    end    
    
    % Draw Sigma_(t)
    % --------------
    % In order to draw the log-volatilies, substract the mean and variance
    % of the 7-component mixture of Normal approximation to the measurement
    % error covariance
    vart = zeros(t*M,M);
    yss1 = zeros(t,M);
    for i = 1:t
        for j = 1:M
            imix = statedraw(i,j);
            vart((i-1)*M+j,j) = u2_s(imix);
            yss1(i,j) = yss(i,j) - m_s(imix) + 1.2704;
        end
    end
    
    % Sigtdraw is a draw of the diagonal log-volatilies, which will give us SIGMA(t)
    
    if red_dim==1 && Xi_rand==0
       Sigtdraw_old=Sigtdraw;
    end
    
    switch(diag_sigma)
        case 1
            Sigtdraw=zeros(M,t);
    
            for k=1:M
            temp_vart=vart(:,k);
            index_vart=find(temp_vart~=0);
            temp_vart=temp_vart(index_vart);
            [Sigtdraw_block,log_lik3] = carter_kohn1(yss1(:,k)',ones(t,1),temp_vart,Wdraw(k,k),1,1,t,sigma_prmean(k),sigma_prvar(k,k),ones(t,1));
            Sigtdraw(k,:)=Sigtdraw_block;
            end
        
        case 2
            [Sigtdraw,log_lik3] = carter_kohn1(yss1',Zs,vart,Wdraw,M,M,t,sigma_prmean,sigma_prvar,TVP_Sigma*ones(t,1));    
    end
    
    % Draws in Sigtdraw are in logarithmic scale (log-volatilies). Create 
    % original standard deviations of the VAR covariance matrix
    sigtemp = eye(M);
    sigt = zeros(M*t,M);
    
    if red_dim==1 && Xi_rand==0    
        Sigtdraw_a=zeros(M,t);
    end
    
    for i = 1:t
         if red_dim==1 && Xi_rand==0
            at=chol(reshape(S_A*Atdraw(:,i)+s_A,M,M)*Z((i-1)*M+1:i*M,:)*Z((i-1)*M+1:i*M,:)'*reshape(S_A*Atdraw(:,i)+s_A,M,M)')';            
             sigtemp = at*diag(exp(0.5*Sigtdraw(:,i)));
             sigt((i-1)*M+1:i*M,:) = sigtemp;
         else
            for j = 1:M
                sigtemp(j,j) = exp(0.5*Sigtdraw(j,i));
            end
            sigt((i-1)*M+1:i*M,:) = sigtemp;
        end
    end

    %=====| Draw W, the covariance of SIGMA(t) (from iWishart)
    % Get first differences of Sigtdraw to compute the SSE
    switch(diag_sigma)
        case 1    
            Wdraw=zeros(M,M);
    
            for k=1:M    
                Sigttemp = Sigtdraw(k,2:t)' - Sigtdraw(k,1:t-1)';
                sse_2 = 0;
                for i = 1:t-1
                    sse_2 = sse_2 + Sigttemp(i,:)'*Sigttemp(i,:);
                end
                Winv = inv(sse_2 + W_prmean);
                Winvdraw = wish(Winv,t+W_prvar);
                Wdraw(k,k) = inv(Winvdraw);  % this is a draw from W    
            end
        
        case 2
            Sigttemp = Sigtdraw(:,2:t)' - Sigtdraw(:,1:t-1)';

            sse_2 = zeros(M,M);
            for i = 1:t-1
                sse_2 = sse_2 + Sigttemp(i,:)'*Sigttemp(i,:);
            end
            Winv = inv(sse_2 + W_prmean);
            Winvdraw = wish(Winv,t+W_prvar);            
            Wdraw = inv(Winvdraw);  % this is a draw from W
    end