function [LogL,cost_flag]=loglike_svar_AB(theta,varargin)

% Constructs the Log-likelihood function for a Structural VAR
% Fernando Pérez Forero (September-2011)

M=varargin{1};
numa = varargin{2};
omega1=varargin{3};
t=varargin{4};
s_A=varargin{5};
S_A=varargin{6};
numa_c = varargin{7};

theta_a=theta(1:numa)';
A=reshape(S_A*theta_a+s_A,M,M);

if numa_c==0
    B=eye(M);
else
    s_B=varargin{8};
    S_B=varargin{9};
    theta_b=theta(numa+1:numa+numa_c)';
    B=reshape(S_B*theta_b+s_B,M,M);
end

C0=B\A;

sigma_index= find(eye(M)==1);
D=zeros(M,M);
temp2=length(sigma_index);
for i=1:temp2 
    D(sigma_index(i))=exp(theta(numa+numa_c+i));
end

LogL=-0.5*(t*M)*log(2*pi)+0.5*t*log(det(C0'*inv(D*D')*C0))-0.5*t*trace(C0'*inv(D*D')*C0*omega1);

 LogL=-LogL;
 
 if imag(LogL)==0
     cost_flag=1;
 else
     cost_flag=[];
 end

end




