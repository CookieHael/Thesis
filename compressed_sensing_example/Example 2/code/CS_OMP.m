%% CS_OMP  Algorithm
%-------------------------------------------------------------------------------------%
%  CS_OMP  Algorithm (Õý½»Æ¥Åä×·×Ù·¨ Orthogonal Matching Pursuit)   
%  ÊäÈë£ºy---²âÁ¿ÐÅºÅ  M X 1
%           A---»Ö¸´¾ØÕó  M X N
%           K---µü´ú´ÎÊý
% Êä³ö £ºtheta---¹À¼ÆµÄÏ¡ÊèÏòÁ¿ N X 1
%            erro_rn---Ã¿´Îµü´úµÄÎó²î
%  ±à³ÌÈË£º ºÎÁõ                                    Email: aresmiki@163.com
%  ±à³ÌÊ±¼ä£º2017Äê04ÔÂ26ÈÕ  Î÷ÄÏ½»Í¨´óÑ§Ç£Òý¶¯Á¦¹ú¼ÒÖØµãÊµÑéÊÒ
%                                        SWJTU  TPL
%  ²Î¿¼ÎÄÏ×£ºJoel A. Tropp and Anna C. Gilbert 
%  Signal Recovery From Random Measurements Via Orthogonal Matching
%  Pursuit£¬IEEE TRANSACTIONS ON INFORMATION THEORY, VOL. 53, NO. 12,
%------------------------------------------------------------------------------------------%
%%   
function [ theta,erro_rn ] = CS_OMP( y,A,K )
N=max(size(A));
M=min(size(A));
theta=zeros(1,N);   %  ´ýÖØ¹¹µÄÏòÁ¿    
Base_t=[];              %  ¼ÇÂ¼»ùÏòÁ¿µÄ¾ØÕó
r_n=y;                  %  ²Ð²îÖµ
for times=1:K;                                    %  µü´ú´ÎÊý(ÓÐÔëÉùµÄÇé¿öÏÂ,¸Ãµü´ú´ÎÊýÎªK)
    for col=1:N;                                  %  »Ö¸´¾ØÕóµÄËùÓÐÁÐÏòÁ¿
        product(col)=abs(A(:,col)'*r_n);          %  »Ö¸´¾ØÕóµÄÁÐÏòÁ¿ºÍ²Ð²îµÄÍ¶Ó°ÏµÊý(ÄÚ»ýÖµ) 
    end
    [val,pos]=max(product);                       %  ×î´óÍ¶Ó°ÏµÊý¶ÔÓ¦µÄÎ»ÖÃ£¬valÖµ£¬posÎ»ÖÃ
    Base_t=[Base_t,A(:,pos)];                       %  ¾ØÕóÀ©³ä£¬¼ÇÂ¼×î´óÍ¶Ó°µÄ»ùÏòÁ¿
    A(:,pos)=zeros(M,1);                          %  Ñ¡ÖÐµÄÁÐÖÃÁã£¨ÊµÖÊÉÏÓ¦¸ÃÈ¥µô£¬ÎªÁË¼òµ¥ÎÒ°ÑËüÖÃÁã£©
    aug_y=(Base_t'*Base_t)^(-1)*Base_t'*y;   %  ×îÐ¡¶þ³Ë,Ê¹²Ð²î×îÐ¡
    r_n=y-Base_t*aug_y;                            %  ²Ð²î
    erro_rn(times)=norm(r_n,2);      %µü´úÎó²î
    pos_array(times)=pos;                         %  ¼ÍÂ¼×î´óÍ¶Ó°ÏµÊýµÄÎ»ÖÃ
    if erro_rn(times)<1e-6 %
            break; %Ìø³öforÑ­»·
    end
end
theta(pos_array)=aug_y;                           %  ÖØ¹¹µÄÏòÁ¿
end