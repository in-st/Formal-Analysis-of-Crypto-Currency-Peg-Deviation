% This code produces the IRF plots.
% This code is based on "Estimating overidentified, non-recursive, time varying coefficients structural VARs"... 
% ...by Fabio Canova and Fernando J. Pérez Forero (Quantitative Economics, Volume 6, Issue 2 (July 2015)).

% Modified by
% Huachen Li

% -------------------------------------------------------------------------
    seed=10101010101010101;
    randn('state',seed);
    rand('twister',seed);
    
    shocks=1;     
    responses=3;
    nhor = 30; 
    Ndraws=nrep/nthin;    
    impresp = zeros(size(responses,2),size(shocks,2),t,nhor,Ndraws);
    impresp_temp= zeros(size(responses,2),size(shocks,2),t,nhor);    
    
for k=1:Ndraws
    
    R=round(unifrnd(1,nrep/nthin));    
    
    biga = zeros(M*p,M*p);
    
    for j = 1:p-1
        biga(j*M+1:M*(j+1),M*(j-1)+1:j*M) = eye(M);
    end
    
        At_posttemp=At_post{R};
        if TVP_Beta==1 && red_dim==1 && Xi_rand==0
            Bt_posttemp=Xi_draw*Thetat_post{R};
        else
            Bt_posttemp=Bt_post{R};
        end
        Sigt_temp=Sigt_post{R};
    
for i = 1:t 
    
    capatemp=reshape(S_A*At_posttemp(:,i)+s_A,M,M);          
    if TVP_Beta==1 && red_dim==1 && Xi_rand==0
        at=chol(capatemp*Z((i-1)*M+1:i*M,:)*Z((i-1)*M+1:i*M,:)'*capatemp')';            
        stem = at*diag(Sigt_temp(:,i));
        Ht = inv(capatemp)*stem*stem'*inv(capatemp)'+Z((i-1)*M+1:i*M,:)*Z((i-1)*M+1:i*M,:)';
        Hsd = chol(Ht)';
    else
        stem = diag(Sigt_temp(:,i));
        Ht = Hsd*Hsd';
        Hsd = capatemp\stem;            
    end
        
    bbtemp = Bt_posttemp(M*constant+1:K,i); 
    splace = 0;
    for ii = 1:p
        for iii = 1:M
            biga(iii,(ii-1)*M+1:ii*M) = bbtemp(splace+1:splace+M,1)';
            splace = splace + M;
        end
    end

    diagonal = diag(diag(Hsd));
    Hsd = Hsd*inv(diagonal);   
    temp_imp=Hsd;
    impresp_temp(:,:,i,1) = temp_imp(responses,shocks); 
    
    bigai = biga;
    for j = 1:nhor-1
        temp_imp=bigj*bigai*bigj'*Hsd;
       impresp_temp(:,:,i,j+1) = temp_imp(responses,shocks);  
       %impresp_temp(:,:,i,j+1) = temp_imp(responses,shocks)+impresp_temp(:,:,i,j);          
        bigai = bigai*biga;
    end
end 
    impresp(:,:,:,:,k)=impresp_temp;
end

clear impresp_temp
impresp=sort(impresp,5);
impresp=permute(impresp,[5 1 2 3 4]);
qus = 0.5;    
impresp=squeeze(quantile(impresp,qus));

    figure     
    set(gcf,'Color',[1 1 1])
    set(gcf,'defaultaxesfontsize',30)
    TT=2017+196/365:1/365:2020+182/365;
    HH=1:nhor; 
    surf(HH,TT,impresp(:,1:end)),axis tight
    clear At_posttemp Bt_posttemp Sigt_temp   