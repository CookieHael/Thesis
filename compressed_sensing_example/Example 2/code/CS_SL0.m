%% CS_SL0  Smoothed l0-Pseudonorm Algorithm
% ÊäÈë£ºy---²âÁ¿ÐÅºÅ M X 1
%          A---»Ö¸´¾ØÕó M X N
%          deltaT---×îÐ¡µÄdelta ÏÞÖÆdelta²»ÄÜÌ«Ð¡,deltaÔ½Ð¡Ô½½Ó½üÓÚl0·¶Êý
%          r---deltaµÄÊÕËõ²½³¤
%          mu---ºÜÐ¡µÄ¸üÐÂ²ÎÊý£¬¼ûÎÄÕÂ
%          L---Ñ­»·´ÎÊý Ä¬ÈÏÈ¡Öµ3
% Êä³ö£ºxs---»Ö¸´µÄÐÅºÅ N X 1
%          valL0---Ã¿´Îµü´úµÄÏ¡Êè¶È×·×ÙÇé¿ö
% 
%  minimize ||x||_0
%  subject to Ax=y
%  ±à³ÌÈË£º ºÎÁõ                                    Email: aresmiki@163.com
%  ±à³ÌÊ±¼ä£º2017Äê04ÔÂ30ÈÕ  Î÷ÄÏ½»Í¨´óÑ§Ç£Òý¶¯Á¦¹ú¼ÒÖØµãÊµÑéÊÒ
%                                        SWJTU  TPL
%  ²Î¿¼ÎÄÏ×£ºH. Mohimani, M. Babie-Zadeh, and C. Jutten,
% ¡°A fast approach for overcomplete sparse decomposition based on smoothed l0-norm,¡±
% IEEE Trans. Signal Process., vol. 57, no. 1, pp. 289-301, Jan. 2009.
%---------------------------------------------------------------------------------------------------------------------%
%%
function [xs,valL0]=CS_SL0(y,A,deltaT,r,mu,L);
% rÎªdeltaËõÐ¡²½³¤£¬<1¡£
y=y(:);
N=max(size(A));
M=min(size(A));
%% ³õÊ¼»¯Ò»Ð©²ÎÊý
pinv_A=pinv(A);
% pinv_A=A'*inv(A*A');
xs=pinv_A*y;  %µÃµ½Ò»¸ö×îÐ¡¶þ³Ë½â×÷Îª³õÊ¼Öµ
delta=2*max(abs(xs));  %ÎÄÕÂÌá³ö2-4±¶£¬Èç¹ûdelta>4*max(abs(ys)),exp(-s.^2/(2*delta.^2))=1
k=0; %¼ÇÂ¼Ñ­»·´ÎÊý
valL0(1)=0;
%  maximizing F_delta(x)=sum(f_delta(x(i)))=sum(exp(-x(i).^2)/(2*delta.^2))
%  subject Ax=y;
%  ÓÃLagrangian ÍÆµ¼£ºL(x,lambda)= F_delta(x)-lambda'*(Ax-y)
%  ¶ÔxºÍlambda·Ö±ðÇóµ¼ÊýµÃµ½KKTÌõ¼þ
%  xÆ«µ¼Êý£º[x(1)exp(-x(1).^2)/(2*delta.^2)),...,x(i)exp(-x(i).^2)/(2*delta.^2))-A'*lambda                (1)
%  lambdaÆ«µ¼Êý£ºAx-y=0
% Ax=yµÄ×îÐ¡¶þ³Ë½âÎª£ºx=pinv(A)*y,Ïàµ±ÓÚ¶ÔÓÅ»¯·½³Ì
% maximizing 0.5*x*x'
% subject Ax=y
% ÓÃLagrangian ÍÆµ¼ L(x,lambda)= 0.5*x*x'-lambda'*(Ax-y)
%  xÆ«µ¼Êý£º[exp(-x(1).^2)/(2*delta.^2)),...,exp(-x(i).^2)/(2*delta.^2))-A'*lambda                (2)
%  lambdaÆ«µ¼Êý£ºAx-y=0
%  (1)ºÍ(2)¶Ô±È£¬·¢ÏÖdeltaÇ÷½üÎÞÇî´ó£¬exp((-x.^2)/(2*delta.^2))=1£¬Á½¸öÓÅ»¯ÎÊÌâµÄ½â¼¸ºõÒ»Ñù
%  delta>>max(xs),xsÊÇ·½³ÌµÄ×îÐ¡¶þ³Ë½â
while delta>deltaT
    k=k+1;
for i=1:L   %L´Î×îËÙÉÏÉýËã·¨
    t_delta=xs.*exp(-abs(xs).^2/(2*delta^2));  %¸üÐÂÌÝ¶ÈÖµ
    xs=xs-mu*t_delta;  %¹Ì¶¨²½³¤µÄÌÝ¶ÈÉÏÉý£¬»ñÈ¡ÉÏÉýºóµÄº¯ÊýÖµ
    xs=xs-pinv_A*(A*xs-y); %Í¶Ó°µ½¿ÉÐÐ¼¯ÉÏ
end
   valL0(k+1)=N-sum(exp(-abs(xs).^2/(2*deltaT^2)));  %ÒªÓÃ´Ë¹«Ê½µÃµ½ºÃµÄÖµdeltaT²»ÄÜÌ«Ð¡
   delta = delta * r;  %ÊÕËõ²½³¤,Öð½¥ËõÐ¡delta£¬Ô½Ð¡Ô½ºÃ£¬Ö»Òª²»Ð¡ÓÚdeltaT
  if (abs(valL0(k+1)-valL0(k))<1e-4);
     valL0=valL0(2:end);
     break;
 end
end
valL0=valL0(2:end);




