%%  CS_SBIL1  L1_SplitBregmanIteration
% ÊäÈë£ºf---²âÁ¿ÐÅºÅ M X 1
%          A---»Ö¸´¾ØÕó M X N
%          mu---²½³¤
%          lambda---¹æÔò»¯Òò×Ó
%          Niter---×î´óµü´ú´ÎÊý
% Êä³ö£ºu---»Ö¸´µÄÐÅºÅ N X 1
%    
% 
%  minimize ||x||_1
%  subject to Ax-y=0
%  ±à³ÌÈË£º ºÎÁõ                                    Email: aresmiki@163.com
%  ±à³ÌÊ±¼ä£º2017Äê04ÔÂ30ÈÕ  Î÷ÄÏ½»Í¨´óÑ§Ç£Òý¶¯Á¦¹ú¼ÒÖØµãÊµÑéÊÒ
%                                        SWJTU  TPL
% ²Î¿¼ÎÄÏ×£ºYin W, Osher S, Goldfarb D, et al. 
% Bregman Iterative Algorithms for L1 Minimization with Applications to Compressed Sensing[J].
% Siam Journal on Imaging Sciences, 2008, 1(1):143-168.
%---------------------------------------------------------------------------------------------------------------------%


function u=CS_SBIL1(f,A,mu,lambda,Niter)

N=max(size(A));

d=zeros(N,1);
b=zeros(N,1);
u=zeros(N,1);

Z=zeros(N,1);
Ft=mu*A'*f;
IV=inv(mu*(A'*A)+lambda*eye(N));

err=norm(f,2);
tol=1e-3*err;
K=0;
while ((err>tol) && (K<Niter)),
    K=K+1;
    up=u;
    u=IV*(Ft+lambda*(d-b));
    tmp=u+b;
    d=sign(tmp).*max(Z,abs(tmp)-1/lambda);
    b=tmp-d;
    err=norm(u-up,2);
end