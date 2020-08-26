%% This code executes the Carter and Kohn algorithm.
% This code is based on "Estimating overidentified, non-recursive, time varying coefficients structural VARs"... 
% ...by Fabio Canova and Fernando J. P�rez Forero (Quantitative Economics, Volume 6, Issue 2 (July 2015)).
% Modified by
% Huachen Li

function [bdraw,log_lik] = carter_kohn1beta(y,Z,Ht,Qt,m,p,t,B0,V0,kdraw,cc)
% Carter and Kohn (1994), On Gibbs sampling for state space models.

% Kalman Filter
bp = B0;
Vp = V0;
bt = zeros(t,m);
Vt = zeros(m^2,t);
log_lik = 0;
cc = cc;
for i=1:t
    R = Ht((i-1)*p+1:i*p,:);
    H = Z((i-1)*p+1:i*p,:);
    % F = eye(m);
    cfe = y(:,i) - H*bp;   % conditional forecast error
    f = H*Vp*H' + R;    % variance of the conditional forecast error
    %inv_f = inv(f);
    inv_f = H'/f;
    %log_lik = log_lik + log(det(f)) + cfe'*inv_f*cfe;
    btt = bp + Vp*inv_f*cfe;
    %btt = bp + Vp*H'*(f\cfe);
    Vtt = Vp - Vp*inv_f*H*Vp;
%    Vtt = Vp - Vp*H'*(f\H)*Vp;
    if i < t
        bp = btt;
        Vp = Vtt + kdraw(i,:)*Qt;
    end
    bt(i,:) = btt';
    Vt(:,i) = reshape(Vtt,m^2,1);
end

% draw Sdraw(T|T) ~ N(S(T|T),P(T|T))
bdraw = zeros(t,m);
bdraw(t,:) = mvnrnd(btt,cc*Vtt,1);

% Backward recurssions
for i=1:t-1
    bf = bdraw(t-i+1,:)';
    btt = bt(t-i,:)';
    Vtt = reshape(Vt(:,t-i),m,m);
    f = Vtt + kdraw(t-i,:)*Qt;
    %inv_f = inv(f);
    inv_f = Vtt/f;
    cfe = bf - btt;
    bmean = btt + inv_f*cfe;
    %bmean = btt + Vtt*f\cfe;
    bvar = Vtt - inv_f*Vtt;
    %bvar = Vtt - Vtt*f\Vtt;
    bdraw(t-i,:) = mvnrnd(bmean,cc*bvar,1); 
end
bdraw = bdraw';