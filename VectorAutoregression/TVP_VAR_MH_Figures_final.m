% "A general algorithm for estimating SVARs"
% ------------------------------------------------------------------------------------
% This code reproduces main figures from Canova and Pérez Forero (2013)    

% This file calls the output estimation workspace:
% results_TVP_VAR_MH_(...).mat

% Cell evaluation usage is recommended(right click/Evaluate Current Cell)
  
%%  1. Plot a particular shock for a particular date with error bands
    
    % Shock
    shocks=4;     
    % Responses
    responses=1:M;
    
    %Insert date
    %------------
    yy=1990;   %year
    qq=1;      %month     
    date=yy+round2(qq/4,0.0001);
    
    nhor = 20;  % Impulse response horizon
    
    Ndraws=nrep/nthin;
        
    impresp = zeros(size(responses,2),size(shocks,2),size(date,1),nhor,Ndraws);
    impresp_temp= zeros(size(responses,2),size(shocks,2),size(date,1),nhor);    
    L=1;
    
for k=1:Ndraws
    
    impresp_tempa= zeros(size(responses,2),size(shocks,2),size(date,1),nhor);    
    
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
    
    for i = 1:length(date) %Get impulses recurssively for each time period
    
        for l=1:L

            capatemp=reshape(S_A*At_posttemp(:,find(yearlab==date(i)))+s_A,M,M);                
                        
            if TVP_Beta==1 && red_dim==1 && Xi_rand==0
                at=chol(capatemp*Z((find(yearlab==date(i))-1)*M+1:find(yearlab==date(i))*M,:)*Z((find(yearlab==date(i))-1)*M+1:find(yearlab==date(i))*M,:)'*capatemp')';            
                stem = at*diag(Sigt_temp(:,find(yearlab==date(i))));
                Ht = inv(capatemp)*stem*stem'*inv(capatemp)'+Z((find(yearlab==date(i))-1)*M+1:find(yearlab==date(i))*M,:)*Z((find(yearlab==date(i))-1)*M+1:find(yearlab==date(i))*M,:)';
                Hsd = chol(Ht)';
            else
                stem = diag(Sigt_temp(:,find(yearlab==date(i))));
                Ht = Hsd*Hsd';
                Hsd = capatemp\stem;            
            end
            bbtemp = Bt_posttemp(M*constant+1:K,find(yearlab==date(i)));  % get the draw of B(t) at time i=1,...,T  (exclude intercept)
            splace = 0;
            for ii = 1:p
                for iii = 1:M
                    biga(iii,(ii-1)*M+1:ii*M) = bbtemp(splace+1:splace+M,1)';
                    splace = splace + M;
                end
            end
        % ------------Identification code:                
        % St dev matrix for structural VAR
            diagonal = diag(diag(Hsd));
            Hsd = Hsd*inv(diagonal);    % Unit initial shock
     
        % Now get impulse responses for 1 through nhor future periods
            temp_imp=Hsd;
            impresp_temp(:,:,i,1) = temp_imp(responses,shocks); % First shock is the Cholesky of the VAR covariance
    
        bigai = biga;
        for j = 1:nhor-1
            temp_imp=bigj*bigai*bigj'*Hsd;
            impresp_temp(:,:,i,j+1) = temp_imp(responses,shocks);    
            bigai = bigai*biga;
        end
        
            impresp_tempa=impresp_tempa+impresp_temp;
        
        end
        
        impresp_tempa=impresp_tempa/L;
        
    end %END geting impulses for each time period
    impresp(:,:,:,:,k)=impresp_tempa;
end

clear impresp_temp

impresp=sort(impresp,5);
impresp=permute(impresp,[5 1 2 3 4]);
qus = [0.16 0.50 0.84];
impresp=squeeze(quantile(impresp,qus));

LW=3;
FS=30;

gr_size2=ceil(length(responses)/3);

     figure
     set(0,'DefaultAxesColorOrder',[1 0 0;0 0 1;1 0 0],...
        'DefaultAxesLineStyleOrder','-|--|-|--')
     set(gcf,'Color',[1 1 1])
     set(gcf,'defaultaxesfontsize',FS-10)
     for i=1:length(responses)
        subplot(gr_size2,3,i)
        plot(zeros(1,nhor),'k')
        hold on
        plot(squeeze(impresp(:,i,:))','LineWidth',LW)
        %title(sprintf('%s after %s, %.0f q%.0f',char(variables(i)),char(variables(shocks)),floor(date),4*(date-floor(date))),'Fontsize',FS)
        title(sprintf('%s',char(variables(i))),'Fontsize',FS)
        xlim([1 nhor])
        set(gca,'XTick',0:nhor/4:nhor)
     end    
   
    clear At_posttemp Bt_posttemp Sigt_temp
    
%% 2. Plot a particular shock for different dates without error bands
    
    % Shock
    shocks=4;     
    % Responses
    responses=1:M;
        
    %Insert dates
    %------------
    yy=1975;   %year
    qq=1;      %quarter     
    date(1,1)=yy+round2(qq/4,0.0001);
    
    yy=1981;   %year
    qq=1;      %quarter     
    date(2,1)=yy+round2(qq/4,0.0001);
    
    yy=1990;   %year
    qq=1;      %quarter     
    date(3,1)=yy+round2(qq/4,0.0001);
    
    yy=2005;   %year
    qq=1;      %quarter     
    date(4,1)=yy+round2(qq/4,0.0001);
    
    nhor = 20;  % Impulse response horizon
    
    Ndraws=nrep/nthin;    
    impresp = zeros(size(responses,2),size(shocks,2),size(date,1),nhor,Ndraws);
    impresp_temp= zeros(size(responses,2),size(shocks,2),size(date,1),nhor);    
    
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
        
for i = 1:length(date) %Get impulses recurssively for each time period
    
    capatemp=reshape(S_A*At_posttemp(:,find(yearlab==date(i)))+s_A,M,M);          
    if TVP_Beta==1 && red_dim==1 && Xi_rand==0
        at=chol(capatemp*Z((find(yearlab==date(i))-1)*M+1:find(yearlab==date(i))*M,:)*Z((find(yearlab==date(i))-1)*M+1:find(yearlab==date(i))*M,:)'*capatemp')';            
        stem = at*diag(Sigt_temp(:,find(yearlab==date(i))));
        Ht = inv(capatemp)*stem*stem'*inv(capatemp)'+Z((find(yearlab==date(i))-1)*M+1:find(yearlab==date(i))*M,:)*Z((find(yearlab==date(i))-1)*M+1:find(yearlab==date(i))*M,:)';
        Hsd = chol(Ht)';
    else
        stem = diag(Sigt_temp(:,find(yearlab==date(i))));
        Ht = Hsd*Hsd';
        Hsd = capatemp\stem;            
    end
        
    bbtemp = Bt_posttemp(M*constant+1:K,find(yearlab==date(i)));  % get the draw of B(t) at time i=1,...,T  (exclude intercept)
    splace = 0;
    for ii = 1:p
        for iii = 1:M
            biga(iii,(ii-1)*M+1:ii*M) = bbtemp(splace+1:splace+M,1)';
            splace = splace + M;
        end
    end
    % ------------Identification code:                
    % St dev matrix for structural VAR
     diagonal = diag(diag(Hsd));
     Hsd = Hsd*inv(diagonal);    % Unit initial shock
               
    % Now get impulse responses for 1 through nhor future periods
    temp_imp=Hsd;
    temp_imp=temp_imp(responses,shocks);
    impresp_temp(:,:,i,1) = temp_imp; % First shock is the Cholesky of the VAR covariance
    
    bigai = biga;
    for j = 1:nhor-1
        temp_imp=bigj*bigai*bigj'*Hsd;
        impresp_temp(:,:,i,j+1) = temp_imp(responses,shocks);    
        bigai = bigai*biga;
    end       
    
end %END geting impulses for each time period
    impresp(:,:,:,:,k)=impresp_temp;
end

clear impresp_temp

impresp=sort(impresp,5);
impresp=permute(impresp,[5 1 2 3 4]);

qus = 0.5;
impresp=squeeze(quantile(impresp,qus));
variables2=variables(responses);
LW=3;
FS=30;

gr_size2=ceil(length(responses)/3);

     figure
     set(0,'DefaultAxesColorOrder',[1 0 0;0 0 1;1 0 0],...
        'DefaultAxesLineStyleOrder','-|--|-|--')
     set(gcf,'Color',[1 1 1])
     set(gcf,'defaultaxesfontsize',FS-10)
     for i=1:length(responses)
        subplot(gr_size2,3,i)
        %plot(zeros(1,nhor),'k')
        grid on
        hold on
        plot(squeeze(impresp(i,1,:))','--b','LineWidth',LW)
        hold on
        plot(squeeze(impresp(i,2,:))','-*r','LineWidth',LW)
        hold on
        plot(squeeze(impresp(i,3,:))','-dg','LineWidth',LW)
        hold on
        plot(squeeze(impresp(i,4,:))','-.k','LineWidth',LW)
        title(sprintf('%s',char(variables2(i))),'Fontsize',FS)
        xlim([1 nhor])
        set(gca,'XTick',0:nhor/4:nhor)
     end  
     AX=legend(sprintf('%.0f q%.0f',floor(date(1)),4*(date(1)-floor(date(1)))),...
         sprintf('%.0f q%.0f',floor(date(2)),4*(date(2)-floor(date(2)))),...
         sprintf('%.0f q%.0f',floor(date(3)),4*(date(3)-floor(date(3)))),...
         sprintf('%.0f q%.0f',floor(date(4)),4*(date(4)-floor(date(4))))...
         ,'Location','best');
     LEG = findobj(AX,'type','text');
         
    clear At_posttemp Bt_posttemp Sigt_temp
    
%%  3. Plot differences of Impulse responses across dates (with error bands)
    
    % Shock
    shocks=3;     
    % Responses
    responses=2;
    
    %Insert dates
    %------------
    yy=1981;   %year
    qq=1;      %month     
    date(1,1)=yy+round2(qq/4,0.0001);
    
    yy=2005;   %year
    qq=1;      %month     
    date(2,1)=yy+round2(qq/4,0.0001);
           
    nhor = 20;  % Impulse response horizon
    
    Ndraws=nrep/nthin;    
    impresp = zeros(size(responses,2),size(shocks,2),size(date,1),nhor,Ndraws);
    impresp_temp= zeros(size(responses,2),size(shocks,2),size(date,1),nhor);    
    
    variables3={'y','p','u','mp','md','i'};
    
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
        
for i = 1:length(date) %Get impulses recurssively for each time period
    
    capatemp=reshape(S_A*At_posttemp(:,find(yearlab==date(i)))+s_A,M,M);          
    if TVP_Beta==1 && red_dim==1 && Xi_rand==0
        at=chol(capatemp*Z((find(yearlab==date(i))-1)*M+1:find(yearlab==date(i))*M,:)*Z((find(yearlab==date(i))-1)*M+1:find(yearlab==date(i))*M,:)'*capatemp')';            
        stem = at*diag(Sigt_temp(:,find(yearlab==date(i))));
        Ht = inv(capatemp)*stem*stem'*inv(capatemp)'+Z((find(yearlab==date(i))-1)*M+1:find(yearlab==date(i))*M,:)*Z((find(yearlab==date(i))-1)*M+1:find(yearlab==date(i))*M,:)';
        Hsd = chol(Ht)';
    else
        stem = diag(Sigt_temp(:,find(yearlab==date(i))));
        Ht = Hsd*Hsd';
        Hsd = capatemp\stem;            
    end
        
    bbtemp = Bt_posttemp(M*constant+1:K,find(yearlab==date(i)));  % get the draw of B(t) at time i=1,...,T  (exclude intercept)
    splace = 0;
    for ii = 1:p
        for iii = 1:M
            biga(iii,(ii-1)*M+1:ii*M) = bbtemp(splace+1:splace+M,1)';
            splace = splace + M;
        end
    end
    % ------------Identification code:                
    % St dev matrix for structural VAR
     diagonal = diag(diag(Hsd));
     Hsd = Hsd*inv(diagonal);    % Unit initial shock
     
    % Now get impulse responses for 1 through nhor future periods
    temp_imp=Hsd;
    impresp_temp(:,:,i,1) = temp_imp(responses,shocks); % First shock is the Cholesky of the VAR covariance
    
    bigai = biga;
    for j = 1:nhor-1
        temp_imp=bigj*bigai*bigj'*Hsd;
        impresp_temp(:,:,i,j+1) = temp_imp(responses,shocks);    
        bigai = bigai*biga;
    end
end %END geting impulses for each time period
    impresp(:,:,:,:,k)=impresp_temp;
end

clear impresp_temp

dimpresp=squeeze(impresp(:,:,2,:,:))-squeeze(impresp(:,:,1,:,:));

FS=20;
LW=3;

if length(responses)>1

dimpresp=sort(dimpresp,3);
dimpresp=permute(dimpresp,[3 1 2]);
qus = [0.16 0.5 0.84];
dimpresp=squeeze(quantile(dimpresp,qus));

gr_size2=ceil(length(responses)/3);

     figure
     set(0,'DefaultAxesColorOrder',[1 0 0;0 0 1;1 0 0],...
        'DefaultAxesLineStyleOrder','-|--|-|--')
     set(gcf,'Color',[1 1 1])
     set(gcf,'defaultaxesfontsize',FS)
     for i=1:length(responses)
        subplot(gr_size2,3,i)
        plot(zeros(1,nhor),'k')
        hold on
        plot(squeeze(dimpresp(:,i,:))','LineWidth',LW)
        %title(sprintf('%s after \\epsilon^{%s}_t',char(variables(i)),char(variables3(shocks))),'FontSize',FS)
        title(sprintf('%s',char(variables(i))),'FontSize',FS)
        xlim([1 nhor])
        set(gca,'XTick',0:5:nhor)
     end  
     ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,sprintf('%.0f q%.0f minus %.0f q%.0f',floor(date(2)),12*(date(2)-floor(date(2))),floor(date(1)),4*(date(1)-floor(date(1)))),'HorizontalAlignment','center','VerticalAlignment', 'top','Fontsize',16)
    
else
    
    dimpresp=sort(dimpresp,2);
    dimpresp=permute(dimpresp,[2 1]);
    qus = [0.16 0.5 0.84];
    dimpresp=squeeze(quantile(dimpresp,qus));

    figure
     set(0,'DefaultAxesColorOrder',[1 0 0;0 0 1;1 0 0],...
        'DefaultAxesLineStyleOrder','-|--|-|--')
     set(gcf,'Color',[1 1 1])
     set(gcf,'defaultaxesfontsize',FS)
    plot(zeros(1,nhor),'k')
    hold on
    plot(dimpresp','LineWidth',LW)
    title(sprintf('%s after \\epsilon^{%s}_t - %.0f q%.0f minus %.0f q%.0f',char(variables(responses)),char(variables3(shocks)),floor(date(2)),4*(date(2)-floor(date(2))),floor(date(1)),4*(date(1)-floor(date(1)))),'Fontsize',FS)
    xlim([1 nhor])
    set(gca,'XTick',0:5:nhor)    
end
    
    clear At_posttemp Bt_posttemp Sigt_temp 

    
%%  4. Plot the median response of a variable after particular shock for all dates
    
    % Shock
    shocks=4;     
    % Responses
    responses=2;
    
    yy=1971;   %year
    qq=1;      %quarter     
    date1=yy+round2(qq/4,0.0001);
    
    yy=2005;   %year
    qq=4;      %quarter     
    date2=yy+round2(qq/4,0.0001);
    
    nhor = 20;  % Impulse response horizon
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
    
for i = 1:t %Get impulses recurssively for each time period
    
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
        
    bbtemp = Bt_posttemp(M*constant+1:K,i);  % get the draw of B(t) at time i=1,...,T  (exclude intercept)
    splace = 0;
    for ii = 1:p
        for iii = 1:M
            biga(iii,(ii-1)*M+1:ii*M) = bbtemp(splace+1:splace+M,1)';
            splace = splace + M;
        end
    end
    % ------------Identification code:                
    % St dev matrix for structural VAR
     diagonal = diag(diag(Hsd));
     Hsd = Hsd*inv(diagonal);    % Unit initial shock
     
    % Now get impulse responses for 1 through nhor future periods
    temp_imp=Hsd;
    impresp_temp(:,:,i,1) = temp_imp(responses,shocks); % First shock is the Cholesky of the VAR covariance
    
    bigai = biga;
    for j = 1:nhor-1
        temp_imp=bigj*bigai*bigj'*Hsd;
        impresp_temp(:,:,i,j+1) = temp_imp(responses,shocks);    
        bigai = bigai*biga;
    end
end %END geting impulses for each time period
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
     surf(1:nhor,yearlab(find(yearlab==date1)+1:find(yearlab==date2)),impresp(find(yearlab==date1)+1:find(yearlab==date2),:))
     set(gca,'ylim',[yearlab(find(yearlab==date1)) yearlab(find(yearlab==date2))],'xlim',[1 nhor],'YDir','reverse')
     %title(sprintf('Response of %s to %s',char(variables(responses)),char(variables(shocks))))
     
     clear At_posttemp Bt_posttemp Sigt_temp         
     
%% 5. Variance decomposition (with bands)

    % Shock
    shocks=1:M;     
    % Responses
    responses=1:M;
        
    Ndraws=100;  % Memory!
        
    % Select horizon
    nhor=20;    
    
    % Select final date
    % -----------------
    yy1=2005;   %year
    qq1=4;      %quarter     
    date1=yy1+round2(qq1/4,0.0001);
    
    impresp = zeros(size(responses,2),size(shocks,2),t,nhor,Ndraws);
    impresp_temp= zeros(size(responses,2),size(shocks,2),t,nhor);
    vdecomp= zeros(size(responses,2),nhor,find(yearlab==date1),size(shocks,2),Ndraws);
    vdecomp_temp= zeros(size(responses,2),nhor,find(yearlab==date1),size(shocks,2));
    
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
    
for i = 1:t %Get impulses recurssively for each time period
    
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
        
    bbtemp = Bt_posttemp(M*constant+1:K,i);  % get the draw of B(t) at time i=1,...,T  (exclude intercept)
    splace = 0;
    for ii = 1:p
        for iii = 1:M
            biga(iii,(ii-1)*M+1:ii*M) = bbtemp(splace+1:splace+M,1)';
            splace = splace + M;
        end
    end
    % ------------Identification code:                
    % St dev matrix for structural VAR
     diagonal = diag(diag(Hsd));
     Hsd = Hsd*inv(diagonal);    % Unit initial shock
     
    % Now get impulse responses for 1 through nhor future periods
    temp_imp=Hsd;
    impresp_temp(:,:,i,1) = temp_imp(responses,shocks); % First shock is the Cholesky of the VAR covariance
    
    bigai = biga;
    for j = 1:nhor-1
        temp_imp=bigj*bigai*bigj'*Hsd;
        impresp_temp(:,:,i,j+1) = temp_imp(responses,shocks);    
        bigai = bigai*biga;
    end
end %END geting impulses for each time period
    impresp(:,:,:,:,k)=impresp_temp;
end

for k=1:Ndraws
    
    impresp_temp=impresp(:,:,:,:,k);
    
    t1=ones(nhor,nhor);
    t1=triu(t1);
    t2=kron(t1,ones(M,1));    
    
    for i = 1:t
        temp1=squeeze(impresp_temp(:,:,i,:));
        temp2=[];
        
        for r=1:nhor
            temp2=[temp2,squeeze(temp1(:,:,r))];
        end
        
        temp3=temp2.^2;
        temp4=temp3*t2;
       
        for j=1:M
            temp6=[];              
            for r=1:nhor
                if r==1    
                    temp7=squeeze(temp3(:,(r-1)*M+j));
                    temp6=[temp6,temp7];
                else
                    temp7=temp7+squeeze(temp3(:,(r-1)*M+j));
                    temp6=[temp6,temp7];
                end                              
            end
            vdecomp_temp(:,:,i,j)=temp6./temp4;
        end    
    end 
    vdecomp(:,:,:,:,k)=vdecomp_temp;
end

clear impresp impresp_temp vdecomp_temp temp1 temp2 temp3 temp4 temp6 temp7

vdecomp=sort(vdecomp,5);
vdecomp=permute(vdecomp,[5 1 2 3 4]);
qus = [0.16 0.5 0.84];
vdecomp=squeeze(quantile(vdecomp,qus));

shocks=3;
variables3_plot={'mp'};
variables_plot={'GDP','P','U'};
responses=[4 5 6];
FS=20;

    for i=1:length(variables3_plot)

     figure    
     set(gcf,'Color',[1 1 1])
     set(0,'DefaultAxesColorOrder',[1 0 0;0 0 1;1 0 0;0 0 0],...
        'DefaultAxesLineStyleOrder','-|--|-|--')
     set(gcf,'defaultaxesfontsize',FS)

        for j=1:length(variables_plot)
            subplot(ceil((length(variables_plot))/2)-1,length(variables_plot)-ceil((length(variables_plot)-1)/2)+1,j)            
            surf(yearlab,1:nhor,squeeze(vdecomp(2,responses(j),:,:,shocks(i))),'EdgeColor','none')
            index=find(yearlab==date1);
            set(gca,'xlim',[yearlab(1) yearlab(index)],'ylim',[1 nhor],'YDir','reverse','zlim',[0 1])                        
            title(sprintf('\\epsilon^{%s} to %s',char(variables3_plot(i)),char(variables_plot(j))),'FontSize',FS)
        end
    end

    nhor_v=[2 nhor/2 nhor];
    
    shocks=3;
    variables3_plot={'mp'};
    variables_plot={'GDP','P','U'};
    responses=[4 5 6];
    FS=20;
    LW=2;
    
    for k=1:length(nhor_v)

    for i=1:length(variables3_plot)

     figure    
     set(gcf,'Color',[1 1 1])
     set(0,'DefaultAxesColorOrder',[1 0 0;0 0 1;1 0 0;0 0 0],...
        'DefaultAxesLineStyleOrder','-|--|-|--')
     set(gcf,'defaultaxesfontsize',FS-5)

        for j=1:length(variables_plot)
            subplot(ceil((length(variables_plot))/2)-1,length(variables_plot)-ceil((length(variables_plot)-1)/2)+1,j)            
            plot(yearlab,squeeze(vdecomp(:,responses(j),nhor_v(k),:,shocks(i))),'LineWidth',LW)
            index=find(yearlab==date1);
            set(gca,'xlim',[yearlab(1) yearlab(index)],'ylim',[0 1])                        
            title(sprintf('\\epsilon^{%s} to %s',char(variables3_plot(i)),char(variables_plot(j))),'FontSize',FS-5)
        end
    end
    
    end

%%  6. Plot histograms of block A^T
    
    At_temp=cell2mat(At_post);
    At_temp=reshape(At_temp,numa,t,nrep/nthin);
    
    FS=14;
    % Select parameters from vector alpha_t
    coef=[2 6 11];
    
    periods=[1990.25 2004.25];
    bins=20;

    for i=1:length(coef)
        figure
        set(0,'DefaultAxesColorOrder',[1 0 0;0 0 1;1 0 0;0 0 0],...
            'DefaultAxesLineStyleOrder','-|-|-|-')
        set(gcf,'Color',[1 1 1])
        set(gcf,'defaultaxesfontsize',FS)
        for j=1:length(periods)
            subplot(1,length(periods),j)
            hist(squeeze(At_temp(coef(i),find(yearlab==periods(j)),:)),bins)
            title(sprintf('%.0fq%.0f',floor(periods(j)),4*(periods(j)-floor(periods(j)))),'Fontsize',FS)
        end
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    end
    clear At_temp
    
%%  7. Plot alpha_t
    
    At_temp=cell2mat(At_post);
    At_temp=reshape(At_temp,numa,t,nrep/nthin);
    At_temp=sort(At_temp,3);
    qus=[0.16 0.50 0.84];
    %qus=0.5;
    At_temp=permute(At_temp,[3 1 2]);
    At_temp=quantile(At_temp,qus);
    
    FS=20;
    LW=2;
    gr_size=ceil(numa/3); 
    
    figure
    set(0,'DefaultAxesColorOrder',[1,0,0;0 0 1;1,0,0],...
        'DefaultAxesLineStyleOrder','-|-|-|-')
    set(gcf,'Color',[1 1 1])    
    set(gcf,'defaultaxesfontsize',FS-10)

    for i=1:numa
        subplot(gr_size,3,i)
        plot(yearlab(1:end,:),zeros(1,t),'-k','Linewidth',LW)
        hold on
        plot(yearlab,squeeze(At_temp(:,i,:))','LineWidth',LW-0.5)
        title(sprintf('\\alpha_{%s}',num2str(i)),'FontSize',FS-5)
        xlim([yearlab(1) yearlab(end)])
    end
    
    % Select parameters from vector alpha_t
        
    coef=[4,7,8,10,12]; % Pcom equation
    %coef=[3,6,9];    % Money demand equation
    %coef=11;           % Money supply equation
    %coef=[1,2,5];    % Non-policy block
            
    for i=1:length(coef)
        figure
        set(0,'DefaultAxesColorOrder',[1,0,0;0 0 1;1,0,0],...
        'DefaultAxesLineStyleOrder','-|-|-|-')
        set(gcf,'Color',[1 1 1])        
        set(gcf,'defaultaxesfontsize',FS)
        plot(yearlab,zeros(1,t),'-k')
        hold on
        plot(yearlab,squeeze(At_temp(:,coef(i),:))','LineWidth',LW)
        title(sprintf('\\alpha_{%s,t}',num2str(coef(i))),'FontSize',FS)
        xlim([yearlab(1) yearlab(end)])
        set(gca,'XTick',yearlab(2):5:yearlab(end))    
    end

    gr_size=ceil(length(coef)/2);        
    
    figure
    set(0,'DefaultAxesColorOrder',[1,0,0;0 0 1;1,0,0],...
    'DefaultAxesLineStyleOrder','-|-|-|-')
    set(gcf,'Color',[1 1 1])        
    set(gcf,'defaultaxesfontsize',FS)    
    
    for i=1:length(coef)
        subplot(gr_size,2,i)
        plot(yearlab,zeros(1,t),'-k')
        hold on
        plot(yearlab,squeeze(At_temp(:,coef(i),:))','LineWidth',LW)
        title(sprintf('\\alpha_{%s,t}',num2str(coef(i))),'FontSize',FS)
        xlim([yearlab(1) yearlab(end)])
        set(gca,'XTick',yearlab(2):10:yearlab(end))    
    end    
    
   clear At_temp
        
%% 8. Plot beta
    
    Bt_temp=cell2mat(Bt_post);
    Bt_temp=reshape(Bt_temp,K,t,nrep/nthin);
    Bt_temp=sort(Bt_temp,3);
    qus=[0.16 0.50 0.84];
    %qus=0.5;
    Bt_temp=permute(Bt_temp,[3 1 2]);
    Bt_temp=quantile(Bt_temp,qus);
    FS=20;
    LW=2;    
    % Select specific parameters
    entries=1:constant*M;
    %entries=[entries,M+1:2*M];

    figure
    set(0,'DefaultAxesColorOrder',[1,0,0;0 0 1;1,0,0],...
        'DefaultAxesLineStyleOrder','-|-|-|-')
    set(gcf,'Color',[1 1 1])    
    set(gcf,'defaultaxesfontsize',FS)
    for i=1:length(entries)
        %subplot(constant,M,i)
        subplot(3,round(length(entries)/3),i)
        plot(yearlab,squeeze(Bt_temp(:,entries(i),:))','LineWidth',1.5)
        title(sprintf('B_{%s,t}',num2str(entries(i))),'FontSize',FS)
        xlim([yearlab(1) yearlab(end)])
    end
    
    clear Bt_temp
    
    if red_dim==1
        Thetat_temp=cell2mat(Thetat_post);
        Thetat_temp=reshape(Thetat_temp,dim_Theta,t,nrep/nthin);
        Thetat_temp=sort(Thetat_temp,3);
        qus=[0.16 0.50 0.84];
        %qus=0.5;
        Thetat_temp=permute(Thetat_temp,[3 1 2]);
        Thetat_temp=quantile(Thetat_temp,qus);
        FS=20;
        LW=2;    
        
        % Select specific parameters
        entries=1:dim_Theta;

        figure
        set(0,'DefaultAxesColorOrder',[1,0,0;0 0 1;1,0,0],...
            'DefaultAxesLineStyleOrder','-|-|-|-')
        set(gcf,'Color',[1 1 1])    
        set(gcf,'defaultaxesfontsize',FS-10)
        for i=1:length(entries)
            %subplot(constant,M,i)
            subplot(5,round(length(entries)/5)+1,i)
            plot(yearlab,squeeze(Thetat_temp(:,entries(i),:))','LineWidth',1.5)
            title(sprintf('\\theta_{%s,t}',num2str(entries(i))),'FontSize',FS)
            xlim([yearlab(1) yearlab(end)])
        end
    
        clear Bt_temp       
        
    end
    
    if single_beta==1
        FS=20;
        bins=20;
    
        figure
         set(0,'DefaultAxesColorOrder',[1 0 0;0 0 1;1 0 0],...
            'DefaultAxesLineStyleOrder','-|--|-|--')
         set(gcf,'Color',[1 1 1])
         set(gcf,'defaultaxesfontsize',FS-10)
         hist(Beta_counter/irep*100,bins)
         title('Acceptance rates of B_t','Fontsize',FS)        
    end
    
%%  9. Stochastic volatility
    
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
    
    % Initial date
    % ------------
    yy1=1970;   %year
    qq1=3;      %quarter     
    date1=yy1+round2(qq1/4,0.0001);

    % Final date
    % ------------
    yy2=2005;   %year
    qq2=4;      %quarter     
    date2=yy2+round2(qq2/4,0.0001);    
    
    FS=40;
    LW=3;
    
    variables3={'y','p','u','mp','md','i'};
    
    resp=1:M;
    
    temp=variables3(resp);     
    
    gr_size=ceil(length(temp)/2);    
    
    figure
    set(0,'DefaultAxesColorOrder',[1 0 0;0 0 1;1 0 0],...
        'DefaultAxesLineStyleOrder','-|-|-|--')
    set(gcf,'Color',[1 1 1])
    set(gcf,'defaultaxesfontsize',FS-30)
    for i=1:length(temp)
        subplot(gr_size,2,i)        
        plot(yearlab,squeeze(Sigt_temp(:,resp(i),:))','LineWidth',LW-1.5)
        title(sprintf('\\sigma^{%s}_{t}',char(temp(i))),'FontSize',FS-25)
        index=find(yearlab==date1);
        index2=find(yearlab==date2);
        xlim([yearlab(index) yearlab(index2)])        
    end    

    for i=1:length(temp)
        figure
        set(0,'DefaultAxesColorOrder',[1 0 0;0 0 1;1 0 0],...
        'DefaultAxesLineStyleOrder','-|-|-|--')
        set(gcf,'Color',[1 1 1])
        set(gcf,'defaultaxesfontsize',FS-10)
        plot(yearlab,squeeze(Sigt_temp(:,resp(i),:))','LineWidth',LW)
        title(sprintf('\\sigma^{%s}_{t}',char(temp(i))),'FontSize',FS)
        index=find(yearlab==date1);
        index2=find(yearlab==date2);
        xlim([yearlab(index) yearlab(index2)])        
    end

   clear Sigt_temp
   
%%  10. Plot Data
    
   
    % Initial date
    % ------------
    yy1=1970;   %year
    qq1=3;      %quarter     
    date1=yy1+round2(qq1/4,0.0001);

    % Final date
    % ------------
    yy2=2005;   %year
    qq2=4;      %quarter     
    date2=yy2+round2(qq2/4,0.0001);    
    
    FS=40;
    LW=3;
    
    resp=1:M;
    
    temp=variables(resp);     
    
    gr_size=ceil(length(temp)/2);    
    
    figure
    set(0,'DefaultAxesColorOrder',[1 0 0;0 0 1;1 0 0],...
        'DefaultAxesLineStyleOrder','-|-|-|--')
    set(gcf,'Color',[1 1 1])
    set(gcf,'defaultaxesfontsize',FS-30)
    for i=1:length(temp)
        subplot(gr_size,2,i)        
        plot(yearlab,y(i,:),'-b','LineWidth',LW-1.5)
        title(sprintf('%s',char(temp(i))),'FontSize',FS-25)
        index=find(yearlab==date1);
        index2=find(yearlab==date2);
        xlim([yearlab(index) yearlab(index2)])        
    end    
    
    figure
    set(0,'DefaultAxesColorOrder',[1 0 0;0 0 1;1 0 0],...
        'DefaultAxesLineStyleOrder','-|-|-|--')
    set(gcf,'Color',[1 1 1])
    set(gcf,'defaultaxesfontsize',FS-30)
    plot(yearlab,y(1,:)','-b','LineWidth',LW-1.5)
    hold on
    plot(yearlab,y(2,:)','-db','LineWidth',LW-1.5)
    hold on    
    plot(yearlab,y(3,:)','-*b','LineWidth',LW-1.5)
    hold on        
    plot(yearlab,y(4,:)','--b','LineWidth',LW-1.5)
    hold on            
    plot(yearlab,y(5,:)','-.b','LineWidth',LW-1.5)
    hold on         
    plot(yearlab,y(6,:)','-ob','LineWidth',LW-1.5)
    hold on         
    legend(variables(1:6));
    index=find(yearlab==date1);
    index2=find(yearlab==date2);
    xlim([yearlab(index) yearlab(index2)])            
   
%%  10. Diagnosis of convergence A^T
       
    At_temp=cell2mat(At_post);
    At_temp=reshape(At_temp,numa,t,nrep/nthin);

    period=25; 

    MCMC=squeeze(At_temp(:,period,:));
    n=size(MCMC,1);
    m=size(MCMC,2);
    x=[];

    for j=11:m
        X=diag(cov(MCMC(:,10:j)'));
        x=[x X];
    end

    sqrn=n^.5;
    
    LW=2;
    FS=30;
    
    % Plot rolling variances
    
    figure
    set(gcf,'Color',[1 1 1])
    set(gcf,'defaultaxesfontsize',FS-15)
    for j=1:n;
        subplot(ceil(sqrn),ceil(sqrn)-1,j);
        plot(x(j,:),'Linewidth',2*LW);
        title(sprintf('Var of \\alpha_{%s}',num2str(j)),'FontSize',FS)
        xlim([1 size(x,2)])
    end
    
    % Plot draws for the entire sample
    
    param=6;

    figure
    set(gcf,'Color',[1 1 1])
    set(gcf,'defaultaxesfontsize',FS-15)
    surf(yearlab,1:nrep/nthin,squeeze(At_temp(param,:,:))','EdgeColor','none')
    set(gca,'xlim',[yearlab(1) yearlab(end)],'ylim',[1 nrep/nthin],'YDir','reverse')                        
    
    clear At_temp
    
    Bt_temp=cell2mat(Bt_post);
    Bt_temp=reshape(Bt_temp,K,t,nrep/nthin);

    period=25; 

    MCMC=squeeze(Bt_temp(:,period,:));
    n=size(MCMC,1);
    m=size(MCMC,2);
    x=[];

    for j=11:m
        X=diag(cov(MCMC(:,10:j)'));
        x=[x X];
    end

    sqrn=n^.5;
    
    LW=2;
    FS=30;
    
    % Plot rolling variances
    
    figure
    set(gcf,'Color',[1 1 1])
    set(gcf,'defaultaxesfontsize',FS-23)
    for j=1:n;
        subplot(ceil(sqrn),ceil(sqrn),j);
        plot(x(j,:),'Linewidth',LW);
        title(sprintf('Var B_{%s}',num2str(j)),'FontSize',FS-23)
        xlim([1 size(x,2)])
    end
    
    % Plot draws for the entire sample
    
    param=6;

    figure
    set(gcf,'Color',[1 1 1])
    set(gcf,'defaultaxesfontsize',FS-15)
    surf(yearlab,1:nrep/nthin,squeeze(Bt_temp(param,:,:))','EdgeColor','none')
    set(gca,'xlim',[yearlab(1) yearlab(end)],'ylim',[1 nrep/nthin],'YDir','reverse')                        
    
    clear Bt_temp       
      
%% 11. Convergence of Beta

info=struct('q',0.025,'r',0.025,'s',0.95);

num_param=K;
vnames=cell(num_param,1);
for i=1:num_param
    vnames{i}=num2str(i);
end
Bt_temp=cell2mat(Bt_post);
Bt_temp=reshape(Bt_temp,K,t,nrep/nthin);

IF_beta=zeros(num_param,t);

for i=1:t
    diagn_beta=coda(squeeze(Bt_temp(:,i,:))',vnames,info);
    for j=1:num_param    
        IF_beta(j,i)=1/diagn_beta(j).rne1;
    end
end

clear Bt_temp

plot(reshape(IF_beta,size(IF_beta,1)*size(IF_beta,2),1))

%% 12. Convergence of alpha

info=struct('q',0.025,'r',0.025,'s',0.95);
num_param=numa;
vnames=cell(num_param,1);
for i=1:num_param
    vnames{i}=num2str(i);
end
At_temp=cell2mat(At_post);
At_temp=reshape(At_temp,numa,t,nrep/nthin);

IF_alpha=zeros(num_param,t);

for i=1:t
    diagn_alpha=coda(squeeze(At_temp(:,i,:))',vnames,info);
    for j=1:num_param    
        IF_alpha(j,i)=1/diagn_alpha(j).rne1;        
    end
end

%prt_coda(diagn_alpha)

clear At_temp

plot(reshape(IF_alpha,size(IF_alpha,1)*size(IF_alpha,2),1))

%% 13. Convergence of Sigma

info=struct('q',0.025,'r',0.025,'s',0.95);
num_param=M;
vnames=cell(num_param,1);
for i=1:num_param
    vnames{i}=num2str(i);
end
Sigt_temp=cell2mat(Sigt_post);
Sigt_temp=reshape(Sigt_temp,M,t,nrep/nthin);

IF_sigma=zeros(M,t);

for i=1:t
    diagn_sigma=coda(squeeze(Sigt_temp(:,i,:))',vnames,info);
    for j=1:M    
        IF_sigma(j,i)=1/diagn_sigma(j).rne1;
    end
end

clear Sigt_temp

plot(reshape(IF_sigma,size(IF_sigma,1)*size(IF_sigma,2),1))


%% 14. Convergence of Q

info=struct('q',0.025,'r',0.025,'s',0.95);
num_param=K^2;
vnames=cell(num_param,1);
for i=1:num_param
    vnames{i}=num2str(i);
end

diagn_Q=coda(Q_post',vnames,info);
prt_coda(diagn_Q)

IF_Q=zeros(num_param,1);
tic
for j=1:num_param    
    IF_Q(j,1)=1/diagn_Q(j).rne1;
end
toc
plot(IF_Q)

%% 15. Convergence of Sa

info=struct('q',0.025,'r',0.025,'s',0.95);
num_param=numa^2;
vnames=cell(num_param,1);
for i=1:num_param
    vnames{i}=num2str(i);
end

diagn_Sa=coda(Sa_post',vnames,info);
%prt_coda(diagn_Sa)

IF_Sa=zeros(num_param,1);

for j=1:num_param    
    IF_Sa(j,1)=1/diagn_Sa(j).rne1;
end

plot(IF_Sa)

%% 16. Convergence of W

info=struct('q',0.025,'r',0.025,'s',0.95);
num_param=M;
vnames=cell(num_param,1);
for i=1:num_param
    vnames{i}=num2str(i);
end

diagn_W=coda(W_post',vnames,info);
%prt_coda(diagn_W)

IF_W=zeros(num_param,1);

for j=1:num_param    
    IF_W(j,1)=1/diagn_W(j).rne1;
end

plot(IF_W)
%% 17. Plot RL

temp_RL=[reshape(IF_alpha,size(IF_alpha,1)*size(IF_alpha,2),1);...
    IF_Q;...
    IF_Sa;...    
    IF_W;...
    reshape(IF_sigma,size(IF_sigma,1)*size(IF_sigma,2),1);...
    reshape(IF_beta,size(IF_beta,1)*size(IF_beta,2),1)];

FS=30;

    figure
     set(0,'DefaultAxesColorOrder',[1 0 0;0 0 1;1 0 0],...
        'DefaultAxesLineStyleOrder','-|--|-|--')
     set(gcf,'Color',[1 1 1])
     set(gcf,'defaultaxesfontsize',FS-10)
     plot(temp_RL)
     xlim([1 size(temp_RL,1)])

    % saveas(gcf,'RL_plot', 'emf')
clear temp_RL