    % Btdrawc is a draw of the mean VAR coefficients, B(t)
    
    % Modified version of Hierarchical_TVP_VAR.m from Koop and Korobilis
    % (2009)
    
    %----------------------------------------------------------------------
    % Step 1: Sample beta_t ~ N(beta_post_mean,beta_post_var)
    %----------------------------------------------------------------------
    
    Btdrawc=zeros(K,t);
    
    Qinv=inv(Qdraw);
    
    for i=1:t        
        beta_post_var = inv(Qinv + Z((i-1)*M+1:i*M,:)'*inv((Ht((i-1)*M+1:i*M,:)))*Z((i-1)*M+1:i*M,:));
        beta_post_mean = beta_post_var*(Qinv*(Xi_draw*Thetatdraw(:,i)) + Z((i-1)*M+1:i*M,:)'*inv((Ht((i-1)*M+1:i*M,:)))*y(:,i));
        
        Btdrawc(:,i) = beta_post_mean + chol(beta_post_var)'*randn(K,1);    
    end
    
    switch(stat)
        case 1
        % Accept draw
        Btdraw = Btdrawc;
        case 2
        % or use the code below to check for stationarity
        %Now check for the polynomial roots to see if explosive

        counter = [];
        restviol=0;
        for i = 1:t;

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

            if max(abs(eig(ctemp1)))>0.9999;
                restviol=1;
                counter = [counter ; restviol]; 
            end
        end
        %if they have been rejected keep old draw, otherwise accept new draw 
        if sum(counter)==0
            Beta_counter=Beta_counter+1;
            Btdraw = Btdrawc;
            %disp('I found a keeper!');
        end
    end
    
    %----------------------------------------------------------------------
    % Step 2: Sample Q ~ iW(v_q,S_q)
    %----------------------------------------------------------------------
    for k=1:K
        sse_q = 0;
        for i=1:t
            sse_q =  sse_q + (Btdraw(k,i) - Xi_draw(k,:)*Thetatdraw(:,i))*(Btdraw(k,i) - Xi_draw(k,:)*Thetatdraw(:,i))';        
        end
        v_q = t + Q_prvar;
        S_q = inv(Q_prmean(k,k) + sse_q);
            Qinv = wish(S_q,v_q);
        Qdraw(k,k) = inv(Qinv);
    end

    %----------------------------------------------------------------------
    % Step 3: Sample theta_t using Carter and kohn
    %----------------------------------------------------------------------    
    [Thetatdrawc,log_lik] = carter_kohn_hom2(Btdraw,Xi_draw,Qdraw,Rdraw,dim_Theta,K,t,Theta_0_prmean,Theta_0_prvar);
    
    Thetatdraw = Thetatdrawc;
    %----------------------------------------------------------------------
    % Step 4: Sample R ~ iW(v_r,S_r)
    %----------------------------------------------------------------------
    
    for k=1:dim_Theta
        sse_r = 0;
        theta_temp = Thetatdraw(k,2:t) - Thetatdraw(k,1:t-1);
        for i=1:t-1
            sse_r =  sse_r + (theta_temp(:,i))*(theta_temp(:,i))';
        end
        v_r = t + R_prvar;
        S_r = inv(R_prmean(k,k) + sse_r);
        Rinv = wish(S_r,v_r);
        Rdraw(k,k) = inv(Rinv);
    end
    
    %----------------------------------------------------------------------
    % Step 5: Sample Xi ~ N(Xi_post_mean,Xi_post_var)
    %----------------------------------------------------------------------

    if Xi_rand==1
        for f = dim_Theta+1:K                
            VXipost = inv(inv(Xi_0_prvar) + Thetatdraw*inv(Qdraw(f,f))*Thetatdraw');
            Xi_mean  = VXipost*(inv(Xi_0_prvar)*Xi_0_prmean(f,:)'+inv(Qdraw(f,f))*Thetatdraw*Btdraw(f,:)');        
            Xi_draw(f,:)= mvnrnd(Xi_mean,VXipost,1)';        
        end    
    end
    