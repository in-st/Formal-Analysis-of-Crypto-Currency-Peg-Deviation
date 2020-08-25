    % ------------------------------
    % Draw B and Omega
    % ------------------------------    
        
    VB_OLS=zeros(K,K);
    B_factor=zeros(K,1);
    
    for i = 1:t    
        zhat1 = Z((i-1)*M+1:i*M,:);
        VB_OLS = VB_OLS + zhat1'*inv(Ht((i-1)*M+1:i*M,:))*zhat1;
        B_factor= B_factor+ zhat1'*inv(Ht((i-1)*M+1:i*M,:))*y(:,i);
    end
    
    VBpost=inv(inv(B_0_prvar)+VB_OLS);
    beta_mean  = VBpost*(inv(B_0_prvar)*B_0_prmean+B_factor);
    
    beta_can= mvnrnd(beta_mean,VBpost,1)';
    
    switch(stat)
        case 1
            beta_draw=beta_can;        
        case 2
            biga = zeros(M*p,M*p);    
            for j = 1:p-1
                biga(j*M+1:M*(j+1),M*(j-1)+1:j*M) = eye(M);
            end
            
            bbtemp = beta_can(M*constant+1:K,:);  % get the draw of B(t) at time i=1,...,T  (exclude intercept)
            splace = 0;
            for ii = 1:p
                for iii = 1:M
                    biga(iii,(ii-1)*M+1:ii*M) = bbtemp(splace+1:splace+M,1)';
                    splace = splace + M;
                end
            end
            
            if max(abs(eig(biga)))<0.9999
                Beta_counter=Beta_counter+1;
                beta_draw=beta_can;
            end
    end
    
    Btdraw=repmat(beta_draw,1,t);
