%% CS_IHT  Iterative Hard Thresholding algorithms for compressive sensing
% ÊäÈë£ºy---²âÁ¿ÐÅºÅ M X 1
%          A---»Ö¸´¾ØÕó M X N
%          K---ÐÅºÅµÄÏ¡Êè¶È
% Êä³ö£ºx0---»Ö¸´µÄÐÅºÅ N X 1
%    
% 
%  minimize ||x||_1
%  subject to ||Ax-y||_2<eps
%  ±à³ÌÈË£º ºÎÁõ                                    Email: aresmiki@163.com
%  ±à³ÌÊ±¼ä£º2017Äê04ÔÂ30ÈÕ  Î÷ÄÏ½»Í¨´óÑ§Ç£Òý¶¯Á¦¹ú¼ÒÖØµãÊµÑéÊÒ
%                                        SWJTU  TPL
%  ²Î¿¼ÎÄÏ×£ºBlumensath T, Davies M E. 
% Iterative hard thresholding for compressed sensing[J]. 
% Applied & Computational Harmonic Analysis, 2009, 27(3):265-274.
% Blumensath T, Davies M E.
% Iterative Thresholding for Sparse Approximations[J].
% Journal of Fourier Analysis and Applications, 2008, 14(5):629-654.
%---------------------------------------------------------------------------------------------------------------------%

function x0=CS_IHT(y,A,K)

M=min(size(A));
N=max(size(A));
if nargin<3;
    K=floor(M/4);        %×îÉÙµü´ú´ÎÊý£¬Ò»°ãµÈÓÚÏ¡Êè¶ÈÊý£¬¾­Ñé¹«Ê½³ýÒÔ4£¬Îª±£Ö¤ÖØ¹¹¾«¶È£¬iter¿ÉÒÔÑ¡´óÒ»µã
end
x0=zeros(N,1);         % ³õÊ¼»¯½â¿Õ¼äÏòÁ¿
u=0.5;                       % Ó°ÏìÏµÊý

for times=1:M
    
    x1=A'*(y-A*x0);
    
    x2=x0+u*x1;
    
    [val,pos]=sort(abs(x2),'descend');  % ½µÐòÅÅÁÐ
    
    x2(pos(K+1:end))=0;   % ±£Áô×î´óµÄÇ°iters¸öÊýµÄÊý¾Ý£¬itersÎªÏ¡Êè¶È

    x0=x2;          % ¸üÐÂÖµµ½ÏÂÒ»²½Ñ­»·

end
