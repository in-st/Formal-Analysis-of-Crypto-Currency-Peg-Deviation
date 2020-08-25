    % ------------------------------
    % Draw Sigma
    % ------------------------------

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
    
    % In order to draw the log-volatilies, substract the mean and variance
    % of the 7-component mixture of Normal approximation to the measurement
    % error covariance
    vart = zeros(t*M,M);
    yss1 = zeros(t,M);
    for i = 1:t
        for j = 1:M
            imix = statedraw(i,j);
            vart((i-1)*M+j,j) = u2_s(imix);
            yss1(i,j) = (yss(i,j) - m_s(imix) + 1.2704);            
        end
    end

    sigma_mean=zeros(M,1);
    VSigma_OLS=zeros(M,1);        
    for k=1:M        
        Sigma_factor=0;
        for i = 1:t            
            VSigma_OLS(k,1) = VSigma_OLS(k,1) + inv(vart((i-1)*M+k,k));
            Sigma_factor= Sigma_factor+ inv(vart((i-1)*M+k,k))*yss1(i,k);
        end        
        VSigma_OLS(k,1)=inv(inv(sigma_prvar(k,k))+VSigma_OLS(k,1));
        sigma_mean(k,1)  = VSigma_OLS(k,1)*(inv(sigma_prvar(k,k))*sigma_prmean(k,1)+Sigma_factor);
    end
    
    sigma_draw=zeros(M,1);
    
    for k=1:M        
        sigma_draw(k,1)=mvnrnd(sigma_mean(k,1),VSigma_OLS(k,1),1)';
    end
    
    % Next draw statedraw (chi square approximation mixture component) conditional on Sigtdraw
    % This is used to update at the next step the log-volatilities Sigtdraw
    for jj = 1:M
        for i = 1:t
            for j = 1:numel(m_s)
                temp1= (1/sqrt(2*pi*u2_s(j)))*exp(-.5*(((yss(i,jj) - sigma_draw(jj,1) - m_s(j) + 1.2704)^2)/u2_s(j)));
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

    Sigtdraw=repmat(sigma_draw,1,t);        
    
    sigtemp = eye(M);
    sigt = zeros(M*t,M);    
    
    for i = 1:t
        for j = 1:M
            sigtemp(j,j) = exp(0.5*Sigtdraw(j,i));
        end
        sigt((i-1)*M+1:i*M,:) = sigtemp;
    end
        

    
    %sigma_draw=exp(0.5*sigma_draw);    
    