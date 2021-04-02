%% CS_NSRAL0  Null-Space Reweigthted Approximate l0-Pseudonorm Algorithm
% ÊäÈë£ºy---²âÁ¿ÐÅºÅ M X 1
%          A---»Ö¸´¾ØÕó M X N
%          deltaT---×îÐ¡µÄdelta ÏÞÖÆdelta²»ÄÜÌ«Ð¡,deltaÔ½Ð¡Ô½½Ó½üÓÚl0·¶Êý
%          r---deltaµÄÊÕËõ²½³¤
%          t---ºÜÐ¡µÄÒ»¸öÕý²ÎÊý£¬È·±£¿ªÊ¼µÄdeltaÔÚ×î´óÖµÉÏÉÔÎ¢´óÒ»µã
%          eps---È¨ÖØ¸üÐÂµÄÒ»¸ö½ÏÐ¡µÄÁ¿eps£¬È·±£·ÖÄ¸²»Îª0
% Êä³ö£ºyk---»Ö¸´µÄÐÅºÅ N X 1
%          valL0---Ã¿´Îµü´úµÄÏ¡Êè¶È×·×ÙÇé¿ö
% 
%  minimize ||x||_0
%  subject to Ax=y
%  ±à³ÌÈË£º ºÎÁõ                                    Email: aresmiki@163.com
%  ±à³ÌÊ±¼ä£º2017Äê04ÔÂ30ÈÕ  Î÷ÄÏ½»Í¨´óÑ§Ç£Òý¶¯Á¦¹ú¼ÒÖØµãÊµÑéÊÒ
%                                        SWJTU  TPL
%  ²Î¿¼ÎÄÏ×£ºJ. K. Pant, W.-S. Lu, and A. Antoniou, 
% ¡°Reconstruction of sparse signals by minimizing a re-weighted approximate l0-norm in the
% null space of the measurement matrix,¡± IEEE Inter. Midwest Symp. on Circuits-Syst, pp. 430¨C433, 2010.
%---------------------------------------------------------------------------------------------------------------------%
%%
function [yk,valL0]=CS_NSRAL0(y,A,deltaT,r,t,eps);
% rÎªdeltaËõÐ¡²½³¤£¬<1¡£
y=y(:);
N=max(size(A));
M=min(size(A));
% AÎ¬¶ÈÎªM X N£»ÆäÁã¿Õ¼ä VÎ¬¶ÈÎª N X (N-M)   Áã¿Õ¼äµÄÒ»¸öÏòÁ¿epsigÎª (N-M) X 1
%% ³õÊ¼»¯Null Space Áã¿Õ¼äÖÐÒ»¸öÏòÁ¿£¬ÎªÍ¨½â×ö×¼±¸
epsig=zeros((N-M),1);   %Áã¿Õ¼äÏòÁ¿  
ys=A'*inv(A*A')*y; %Ò»¸ö×îÐ¡¶þ³Ë½âÊÇÌØ½â
w=ones(N,1);     %È¨ÖØ¾ØÕó,¿ªÊ¼Ã¿¸öÈ¨ÖØ¶¼Îª1
delta=max(y)+t;  %a reasonable initial value of theta È·±£¿ªÊ¼º¯ÊýÊÇÍ¹µÄ
k=0;
% V=null(A);   %matlab ¿ÉÒÔÖ±½ÓÃüÁîÇó½â
[Q,R]=qr(A');   %ÇóQR·Ö½â
V=Q(:,M+1:N);  % QµÄ×îºóN-MÊÇ¾ØÕóAµÄÁã¿Õ¼ä¾ØÕó
valL0(1)=0;
 %% ¹¹ÔìÓÅ»¯ÎÊÌâ
while (delta>deltaT)
    k=k+1;  %¼ÇÂ¼µü´ú´ÎÊý
%  epsig=fminunc(@(epsigx) w'*(1-exp(-(ys+V*epsigx).^2./(2*delta.*delta)')),epsig);   %Çó×îÓÅ»¯ÎÊÌâ
 [epsig,val,iters]=bfgs(@(epsigx) w'*(1-exp(-(ys+V*epsigx).^2./(2*delta.*delta)')),...
                   @(epsigx) (V'*(w.*((ys+V*epsigx).*exp(-(ys+V*epsigx).^2./(2*delta.*delta)')))),epsig);
 yk=ys+V*epsig;
 w=1./(abs(yk)+eps);  %È¨ÖØ¸üÐÂ
 delta=r*delta;
 [valL0(k+1)]= ones(1,N)*(1-exp(-(yk).^2./(2*deltaT.*deltaT)'));
 if (abs(valL0(k+1)-valL0(k))<1e-4);
     valL0=valL0(2:end);
     break;
 end
end
valL0=valL0(2:end);
end




