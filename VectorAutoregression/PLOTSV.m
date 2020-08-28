% This code produces the SV plots.
% This code is based on "Estimating overidentified, non-recursive, time varying coefficients structural VARs"... 
% ...by Fabio Canova and Fernando J. Pérez Forero (Quantitative Economics, Volume 6, Issue 2 (July 2015)).

% Modified by
% Huachen Li

% -------------------------------------------------------------------------
    seed=10101010101010101;
    randn('state',seed);
    rand('twister',seed);
    
    if TVP_Beta==1 && red_dim==1 && Xi_rand==0
        Sigt_temp=zeros(M,t,nrep/nthin);
        
        for k=1:nrep/nthin
            R=round(unifrnd(1,nrep/nthin)); 
            At_posttemp=At_post{R};                            
            Sigt_temp_k=Sigt_post{R};
            for i = 1:t
                at=chol(reshape(S_A*At_posttemp(:,i)+s_A,M,M)*Z((i-1)*M+1:i*M,:)*Z((i-1)*M+1:i*M,:)'*reshape(S_A*At_posttemp(:,i)+s_A,M,M)')';
                Sigt_temp_k(:,i)=at*Sigt_temp_k(:,i);                
            end            
            Sigt_temp(:,:,k)=Sigt_temp_k;
        end
    else
        Sigt_temp=cell2mat(Sigt_post);
        Sigt_temp=reshape(Sigt_temp,M,t,nrep/nthin);
    
    end
    Sigt_temp=sort(Sigt_temp,3);
    qus=[0.16 0.50 0.84];
    Sigt_temp=permute(Sigt_temp,[3 1 2]);
    Sigt_temp=quantile(Sigt_temp,qus);        
      
    FS=40;
    LW=3;
    variables3={'f','p','x','v'};
    resp=1:M;
    temp=variables3(resp);     
    
    gr_size=ceil(length(temp)/2);

    TT=2017+196/365:1/365:2020+182/365;    
    
    figure
    set(0,'DefaultAxesColorOrder',[1 0 0;0 0 1;1 0 0;],...
        'DefaultAxesLineStyleOrder','-|-|-|--')
    set(gcf,'Color',[1 1 1])
    set(gcf,'defaultaxesfontsize',FS-30)
    for i=1:length(temp)
        subplot(gr_size,2,i)        
        plot(TT,squeeze(Sigt_temp(:,resp(i),:))','LineWidth',LW-1.5)
        title(sprintf('\\sigma^{%s}_{t}',char(temp(i))),'FontSize',FS-25)    
    end    

   clear Sigt_temp
 