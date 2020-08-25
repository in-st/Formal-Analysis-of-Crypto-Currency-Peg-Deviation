    % Btdrawc is a draw of the mean VAR coefficients, B(t)
    [Thetatdrawc,log_lik] = carter_kohn1(y,Z*Xi_draw,Ht,Qdraw,dim_Theta,M,t,Theta_0_prmean,Theta_0_prvar,ones(t,1));    

    % Recover parameter vector
    
    Btdrawc = Xi_draw*Thetatdrawc;
    
    switch(stat)
        case 1
        % Accept draw
        Btdraw = Btdrawc;
        Thetatdraw=Thetatdrawc;
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
            Thetatdraw=Thetatdrawc;
            %disp('I found a keeper!');
        end
    end

    %=====| Draw Q, the covariance of Theta(t) (from iWishart)
    
    Qdraw=zeros(dim_Theta,dim_Theta);
    
    for k=1:dim_Theta    
        Thetatemp = Thetatdraw(k,2:t)' - Thetatdraw(k,1:t-1)';
        sse_2 = 0;
        for i = 1:t-1
            sse_2 = sse_2 + Thetatemp(i,:)'*Thetatemp(i,:);
        end
        Qinv = inv(sse_2 + Q_prmean(k,k));
        Qinvdraw = wish(Qinv,t+Q_prvar);
        Qdraw(k,k) = inv(Qinvdraw);  % this is a draw from Q_k    
    end
