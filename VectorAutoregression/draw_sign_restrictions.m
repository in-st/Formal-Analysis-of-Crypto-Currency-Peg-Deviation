% Evaluate sign restrictions

if Beta_counter>0
    counter = zeros(t,horiz);
    for i = 1:t
        biga = zeros(M*p,M*p);
        for j = 1:p-1
            biga(j*M+1:M*(j+1),M*(j-1)+1:j*M) = eye(M);
        end

        splace = 0;
        BBtempor = Btdraw(M+1:end,i);
        for ii = 1:p
            for iii = 1:M
                biga(iii,(ii-1)*M+1:ii*M) = BBtempor(splace+1:splace+M,1)';
                splace = splace + M;
            end
        end

        A0_temp=reshape(S_A*Atdraw(:,i)+s_A,M,M);
        Sigmat_temp=diag(exp(0.5*Sigtdraw(:,i)));
        
        H_temp=A0_temp\Sigmat_temp;
        diagonal=diag(diag(H_temp));
        H_temp=diagonal\H_temp;
        bigai=biga;
        temp=bigj*bigai*bigj'*H_temp;
        
        if temp(4,4)>=0 && temp(5,4)<=0 % Liquidity effect
           counter(i,1)=1;
        end
        
        for k=1:horiz-1
            bigai=bigai*biga;
            temp=bigj*bigai*bigj'*H_temp;            
            if temp(4,4)>=0 && temp(5,4)<=0 % Liquidity effect
                counter(i,k+1)=1;
            end
        end
    end
    
    if sum(sum(counter))==t*horiz        
        % Accept draw
        Bt_old=Btdraw;
        At_old=Atdraw;
        Sigt_old=Sigtdraw; 
        sign_counter=sign_counter+1;
    else
        % Reject draw
        Btdraw=Bt_old;
        Atdraw=At_old;
        Sigtdraw=Sigt_old;
    end
else
    % Reject draw
    Btdraw=Bt_old;
    Atdraw=At_old;
    Sigtdraw=Sigt_old;    
end
