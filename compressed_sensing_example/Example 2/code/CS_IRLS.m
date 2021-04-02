%% CS_IRLS  Iteratively reweighted algorithms for compressive sensing
% ÊäÈë£ºy---²âÁ¿ÐÅºÅ M X 1
%          A---»Ö¸´¾ØÕó M X N
%          p---Î±·¶Êý 1--0Ö®¼ä
%          deltaT---×îÐ¡µÄdelta ÏÞÖÆdelta²»ÄÜÌ«Ð¡,deltaÔ½Ð¡Ô½½Ó½üÓÚl0·¶Êý
%          r---deltaµÄÊÕËõ²½³¤
%          mu---ºÜÐ¡µÄ¸üÐÂ²ÎÊý£¬ÕýÊý 2-4È¡Öµ ¼ûÎÄÕÂ
%          L---Ñ­»·´ÎÊý Ä¬ÈÏÈ¡Öµ3
% Êä³ö£ºys---»Ö¸´µÄÐÅºÅ N X 1
%    
% 
%  minimize ||x||_p
%  subject to Ax=y
%  ±à³ÌÈË£º ºÎÁõ                                    Email: aresmiki@163.com
%  ±à³ÌÊ±¼ä£º2017Äê04ÔÂ30ÈÕ  Î÷ÄÏ½»Í¨´óÑ§Ç£Òý¶¯Á¦¹ú¼ÒÖØµãÊµÑéÊÒ
%                                        SWJTU  TPL
%  ²Î¿¼ÎÄÏ×£º Chartrand and W. Yin,
% ¡°Iteratively Reweighted Algorithms for Compressed Sensing,¡± 2008.
%---------------------------------------------------------------------------------------------------------------------%
function ys=CS_IRLS(y,A,p,epsT,r)
N=max(size(A));
M=min(size(A));
y=y(:);

eps=1;
ys=inv(A'*A)*A'*y;  %¹À¼ÆÒ»¸ö×îÐ¡¶þ³Ë½â
% ys=pinv(A)*y;

while (eps>epsT)
    
    w=(ys.^2+eps).^(p/2-1);  %¸üÐÂÈ¨ÖØ
%     w=(abs(ys)+eps).^(p-1);  %¸üÐÂÈ¨ÖØ
    Q=diag(1./w);
    yx=Q*A'*inv(A*Q*A')*y;  %¸üÐÂ
    if(norm(yx-ys,2) < sqrt(eps)*r.^2)
        eps=eps*r;  %¸üÐÂeps  r±¶¸üÐÂ
    end
    ys=yx;
    
end

end
