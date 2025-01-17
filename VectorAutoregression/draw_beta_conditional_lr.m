% -----------------------------------------------------------------------------------------
%   STEP I: Sample B from p(B|y,A,Sigma,V)
% -----------------------------------------------------------------------------------------

% Btdrawc is a draw of the mean VAR coefficients, B(t)
if rho_b==1
    [Btdrawc,log_lik] = carter_kohn1(y,Z,Ht,Qdraw,K,M,t,B_0_prmean,B_0_prvar,ones(t,1));    
else
    [Btdrawc,log_lik] = carter_kohn2(y,Z,Ht,Qdraw,K,M,t,B_0_prmean,B_0_prvar,ones(t,1),rho_b);    
end
    
    switch(stat)
        case 1   % Accept draw
             Btdraw = Btdrawc;
        case 2   % Check for stationarity
        
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
                
                % Compute long run impact matrix

                A0_temp=reshape(S_A*Atdraw(:,i)+s_A,M,M);
                Sigmat_temp=diag(exp(0.5*Sigtdraw(:,i)));
        
                H_temp=A0_temp\Sigmat_temp;                
                
                %    H_temp=Htsd((i-1)*M+1:i*M,:);
        
                % Cumulative impacts matrix
                    temp_lr=bigj*inv(eye(size(ctemp1,1))-ctemp1)*bigj'*H_temp;
                
                    temp_lr(1,4)=0; % Money neutrality
        
                % Adjust draw of Bt
        
                    sumBt=ctemp1(1:M,1:M);
                                    
                    B2temp=eye(M)-sumBt-H_temp/temp_lr;
                            
                    Btdrawc(M+M*M*(p-1)+1:end,i)=reshape(B2temp',M*M,1);
               
                % Companion form again
                
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
                
                if max(abs(eig(ctemp1)))>0.9999
                    restviol=1;
                    counter = [counter ; restviol];
                end
            end
        
            if sum(counter)==0
                Beta_counter=Beta_counter+1;
                Btdraw = Btdrawc; % Accept new draw
            end
    end

    %=====| Draw Q, the covariance of B(t) (from iWishart)
    % Take the SSE in the state equation of B(t)
    Btemp = Btdraw(:,2:t)' - Btdraw(:,1:t-1)';
    sse_2 = zeros(K,K);
    for i = 1:t-1
        sse_2 = sse_2 + Btemp(i,:)'*Btemp(i,:);
    end
    
    % ...and subsequently draw Q, the covariance matrix of B(t)
    Qinv = inv(sse_2 + Q_prmean);
    Qinvdraw = wish(Qinv,t+Q_prvar);
    Qdraw = inv(Qinvdraw);  % this is a draw from Q