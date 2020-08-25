    % ------------------------------
    % Draw alpha and P(alpha)
    % ------------------------------    
    
    yhat = zeros(M,t);
    for i = 1:t
        yhat(:,i) = y(:,i) - Z((i-1)*M+1:i*M,:)*Btdraw(:,i);
    end

    ytildeA=zeros(M,t);
    Zc = zeros(t*M,numa);

    for i = 1:t    
        ytildeA(:,i)=reshape(s_A,M,M)*yhat(:,i);
        Stemp=-kron(yhat(:,i)',eye(M))*S_A;
        Zc((i-1)*M+1:i*M,:) =  Stemp;
        sigatemp = sigt((i-1)*M+1:i*M,:);
        siga2temp(:,:,i) = sigatemp*sigatemp';                
    end
    
    % Draw from proposal density alpha_{i+1} ~ N(alpha_{i},r*P(alpha_{i}))
    alpha_can = mvnrnd(alpha_draw,c_a*P_draw,1)';

    if min(alpha_can > LB)==1
        if min(alpha_can < UB)==1
            lpostcan=0;
            lpostdraw=0;
            for i=1:t
                alpha_err=ytildeA(:,i)-Zc((i-1)*M+1:i*M,:)*alpha_can;
                lpostcan=lpostcan-0.5*alpha_err'*inv(siga2temp(:,:,i))*alpha_err;
                alpha_err=ytildeA(:,i)-Zc((i-1)*M+1:i*M,:)*alpha_draw;
                lpostdraw=lpostdraw-0.5*alpha_err'*inv(siga2temp(:,:,i))*alpha_err;
            end
            
            % Add the Jacobian
            lpostcan=lpostcan + t*log(det(reshape(S_A*alpha_can+s_A,M,M)));            
            lpostdraw=lpostdraw + t*log(det(reshape(S_A*alpha_draw+s_A,M,M)));            

            laccprob = lpostcan-lpostdraw;
        else
            laccprob=-9e+200;
            out_counter=out_counter+1;
        end
    else
        laccprob=-9e+200;
        out_counter=out_counter+1;
    end
    
    %Accept candidate draw with log prob = laccprob, else keep old draw
    if log(rand)<laccprob
        lpostdraw=lpostcan;
        alpha_draw=alpha_can;
        acc_counter_a=acc_counter_a+1;     
    end

    sse_draw=zeros(M,M);
    for i=1:t
        err_alpha=ytildeA(:,i)-Zc((i-1)*M+1:i*M,:)*alpha_draw;    
        sse_draw=sse_draw+(err_alpha*err_alpha');
    end

    P_draw=zeros(numa,numa);
    for i=1:t
        P_draw=P_draw+Zc((i-1)*M+1:i*M,:)'*inv(sse_draw)*Zc((i-1)*M+1:i*M,:);
    end

    P_draw=inv(P_draw); % this is a draw from P  
    
    Atdraw=repmat(alpha_draw,1,t);
    