
function [x,val,k]=bfgs(fun,gfun,x0)
%¹¦ÄÜ£ºÓÃBFGSËã·¨Çó½âÎÞÔ¼ÊøÎÊÌâ£ºmin f(x)
% ÊäÈë£ºx0ÊÇ³õÊ¼µã£¬fun,gfun·Ö±ðÊÇÄ¿±êº¯Êý¼°ÆäÌÝ¶È£»
%vararginÊÇÊäÈë¿É±ä²ÎÊý±äÁ¿£¬¼òµ¥µ÷ÓÃbfgsÊ±¿ÉÒÔºöÂÔËü£¬
% µ«ÊÇÆäËû³ÌÐòÑ­»·µ÷ÓÃÊ±½«»á·¢»ÓÖØÒª×÷ÓÃ
%Êä³ö£ºx,val·Ö±ðÊÇ½üËÆ×îÓÅµãºÍ×îÓÅÖµ£¬kÊÇµü´ú´ÎÊý¡£
% syms x1 x2;
maxk=500;       %¸ø³ö×î´óµü´ú´ÎÊý
rho=0.55; sigma=0.4; epsilon=1e-6;      %¸ø³öÒ»Ð©³£Êý²ÎÊý¼°¾«¶ÈÎó²î
k=0; n=length(x0);
Bk=eye(n);      %Bk=feval('Hesse',x0);
x=x0;
%%
while(k<maxk)
    gk=feval(gfun,x0);      %¼ÆËã¾«¶È µ¼ÊýÖµÌ«Ð¡£¬Ð¡ÓÚepsilon¾Í²»¼ÆËãÁË
    if(norm(gk)<epsilon)
        break;
    end                         %¼ìÑéÖÕÖ¹×¼Ôò
    dk=-Bk\gk;      %½â·½³Ì×é£¬¼ÆËãËÑË÷·½Ïò
    m=0;mk=0;
    while(m<20)     %ÓÃArmijoËÑË÷Çó²½³¤
        newf=feval(fun,x0+rho^m*dk);
        oldf=feval(fun,x0);
        if(newf<oldf+sigma*rho^m*gk'*dk)
            mk=m;
            break;
       end
    m=m+1;
    end
%BFGS½ÃÕý
x=x0+rho^mk*dk;
sk=x-x0;  yk=feval(gfun,x)-gk;
if(yk'*sk>0)
    Bk=Bk-(Bk*sk*sk'*Bk)/(sk'*Bk*sk)+(yk*yk')/(yk'*sk);
end
k=k+1;      x0=x;
% val=feval(fun,x0);
end
%%
val=feval(fun,x0);