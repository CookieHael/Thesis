%% CS_UALP  Minimization of Approximate Lp Pseudonorm Using a Quasi-Newton Algorithm
% ÊäÈë£ºy---²âÁ¿ÐÅºÅ M X 1
%          A---»Ö¸´¾ØÕó M X N
%          deltaT---×îÐ¡µÄdelta ÏÞÖÆdelta²»ÄÜÌ«Ð¡,deltaÔ½Ð¡Ô½½Ó½üÓÚl0·¶Êý
%          r---deltaµÄÊÕËõ²½³¤
%          mu---ºÜÐ¡µÄ¸üÐÂ²ÎÊý£¬ÕýÊý 2-4È¡Öµ ¼ûÎÄÕÂ
%          L---Ñ­»·´ÎÊý Ä¬ÈÏÈ¡Öµ3
% Êä³ö£ºxs---»Ö¸´µÄÐÅºÅ N X 1
%          valL0---Ã¿´Îµü´úµÄÏ¡Êè¶È×·×ÙÇé¿ö
% 
%  minimize ||x||_0
%  subject to Ax=y
%  ±à³ÌÈË£º ºÎÁõ                                    Email: aresmiki@163.com
%  ±à³ÌÊ±¼ä£º2017Äê04ÔÂ30ÈÕ  Î÷ÄÏ½»Í¨´óÑ§Ç£Òý¶¯Á¦¹ú¼ÒÖØµãÊµÑéÊÒ
%                                        SWJTU  TPL
%  ²Î¿¼ÎÄÏ×£ºPant J K, Lu W S, Antoniou A. 
% Unconstrained regularized Lp -norm based algorithm for the reconstruction of sparse signals[C]
% IEEE International Symposium on Circuits and Systems. IEEE, 2011:1740-1743.
%---------------------------------------------------------------------------------------------------------------------%
%% UALPÊÇPant J KÔÚ²©Ê¿ÂÛÎÄÖÐµÄ¶¨Òå£¬URLPÊÇIEEEÂÛÎÄÖÐµÄ¶¨Òå£¬Ëã·¨Ò»Ñù
function yk=CS_UALP(y,A,p)
y=y(:);
N=max(size(A));
M=min(size(A));
% AÎ¬¶ÈÎªM X N£»ÆäÁã¿Õ¼ä VÎ¬¶ÈÎª N X (N-M)   Áã¿Õ¼äµÄÒ»¸öÏòÁ¿epsigÎª (N-M) X 1
%% ³õÊ¼»¯Null Space Áã¿Õ¼äÖÐÒ»¸öÏòÁ¿£¬ÎªÍ¨½â×ö×¼±¸
epsig=zeros((N-M),1);   %Áã¿Õ¼äÏòÁ¿  
ys=A'*inv(A*A')*y; %Ò»¸ö×îÐ¡¶þ³Ë½âÊÇÌØ½â
% k=0;
% V=null(A);   %matlab ¿ÉÒÔÖ±½ÓÃüÁîÇó½â
[Q,R]=qr(A');   %ÇóQR·Ö½â
V=Q(:,M+1:N);  % QµÄ×îºóN-MÊÇ¾ØÕóAµÄÁã¿Õ¼ä¾ØÕó
eps1=sqrt(1-p)*max(abs(ys));
epsT=1e-5;
T=9;
belta=log(eps1/epsT)/(T-1);
for i=2:(T-1)
    eps(i)=exp(-belta*i);
end
eps=[eps1,eps,epsT];
w=ones(N,1);     %È¨ÖØ¾ØÕó,Ã¿¸öÈ¨ÖØ¶¼Îª1,ÎªÁËÇóºÍ
for k=1:T
      epsig=fminunc(@(epsigx) w'*((ys+V*epsigx).^2+eps(k).^2).^(p/2),epsig);   %Çó×îÓÅ»¯ÎÊÌâ
%     [epsig,val,iters]=bfgs(@(epsigx) w'*((ys+V*epsigx).^2+eps(k).^2).^(p/2),...
%                    @(epsigx) p*V'*((((ys+V*epsigx).^2+eps(k).^2).^(p/2-1)).*(ys+V*epsigx)),epsig);             
end
yk=ys+V*epsig;



%%  BFGSÖÐµÄÏßÐÔËÑË÷
%     function [alpha1]=LSBFP(epsigk,xs,V,p,dk,deltaT,epsk)
%         lalpha1=0;deltaA=deltaT+1;
%         while(deltaA>deltaT)
%             gy=(xs+V*epsigk+alpha1*V*dk).^2+epsk*epsk;
%             Geps=-sum(((xs+V*epsigk).*(V*dk).*gy.^(p/2-1)))/sum(V*dk.*gy.^(p/2-1));
%             deltaA=Geps-alpha1;
%             alpha1=Geps;
%         end
%         
%     end

end

