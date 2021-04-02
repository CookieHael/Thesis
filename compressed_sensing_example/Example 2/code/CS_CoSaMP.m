%% CS_CoSaMP  Algorithm 
%-----------------------------------------------------------------------------------------%
%  CS_CoSaMP  Algorithm (Ñ¹Ëõ²ÉÑùÕý½»Æ¥Åä×·×Ù·¨ Orthogonal Matching Pursuit)   
%  ÊäÈë£ºy---²âÁ¿ÐÅºÅ  M X 1
%           A---»Ö¸´¾ØÕó  M X N
%           K---µü´ú´ÎÊý
% Êä³ö £ºtheta---¹À¼ÆµÄÏ¡ÊèÏòÁ¿ N X 1
%            erro_res---Ã¿´Îµü´úµÄÎó²î
%  ±à³ÌÈË£º ºÎÁõ                                    Email: aresmiki@163.com
%  ±à³ÌÊ±¼ä£º2017Äê04ÔÂ26ÈÕ  Î÷ÄÏ½»Í¨´óÑ§Ç£Òý¶¯Á¦¹ú¼ÒÖØµãÊµÑéÊÒ
%                                        SWJTU  TPL
%  ²Î¿¼ÎÄÏ×1£ºNeedell D£¬Tropp J A 
%  CoSaMP£ºIterative signal recovery from incomplete and inaccurate samples[J]£®
% Applied and Computation Harmonic Analysis£¬2009£¬26£º301-321.
% ²Î¿¼ÎÄÏ×2£ºD.Needell, J.A. Tropp£®
% CoSaMP: Iterative signal recoveryfrom incomplete and inaccurate samples[J]. 
% Communications of theACM£¬2010£¬53(12)£º93-100.
%------------------------------------------------------------------------------------------%

%%
function [ theta,erro_res ] = CS_CoSaMP( y,A,K )
    [m,n] = size(y);
    if m<n
        y = y'; 
    end
    [M,N] = size(A); %´«¸Ð¾ØÕóAÎªM*N¾ØÕó
    theta = zeros(N,1); %ÓÃÀ´´æ´¢»Ö¸´µÄtheta(ÁÐÏòÁ¿)
    pos_num = []; %ÓÃÀ´µü´ú¹ý³ÌÖÐ´æ´¢A±»Ñ¡ÔñµÄÁÐÐòºÅ
    res = y; %³õÊ¼»¯²Ð²î(residual)Îªy
    for kk=1:K %×î¶àµü´úK´Î
        %(1) Identification
        product = A'*res; %´«¸Ð¾ØÕóA¸÷ÁÐÓë²Ð²îµÄÄÚ»ý
        [val,pos]=sort(abs(product),'descend');
        Js = pos(1:2*K); %Ñ¡³öÄÚ»ýÖµ×î´óµÄ2KÁÐ
        %(2) Support Merger
        Is = union(pos_num,Js); %Pos_thetaÓëJs²¢¼¯
        %(3) Estimation
        %AtµÄÐÐÊýÒª´óÓÚÁÐÊý£¬´ËÎª×îÐ¡¶þ³ËµÄ»ù´¡(ÁÐÏßÐÔÎÞ¹Ø)
        if length(Is)<=M
            At = A(:,Is); %½«AµÄÕâ¼¸ÁÐ×é³É¾ØÕóAt
        else %AtµÄÁÐÊý´óÓÚÐÐÊý£¬ÁÐ±ØÎªÏßÐÔÏà¹ØµÄ,At'*At½«²»¿ÉÄæ
            if kk == 1
                theta_ls = 0;
            end
            break; %Ìø³öforÑ­»·
        end
        %y=At*theta£¬ÒÔÏÂÇóthetaµÄ×îÐ¡¶þ³Ë½â(Least Square)
        theta_ls = (At'*At)^(-1)*At'*y; %×îÐ¡¶þ³Ë½â
        %(4) Pruning
        [val,pos]=sort(abs(theta_ls),'descend');
        %(5) Sample Update
        pos_num = Is(pos(1:K));
        theta_ls = theta_ls(pos(1:K));
        %At(:,pos(1:K))*theta_lsÊÇyÔÚAt(:,pos(1:K))ÁÐ¿Õ¼äÉÏµÄÕý½»Í¶Ó°
        res = y - At(:,pos(1:K))*theta_ls; %¸üÐÂ²Ð²î 
        erro_res(kk)=norm(res,2);
        if norm(res)<1e-6 %Repeat the steps until r=0
            break; %Ìø³öforÑ­»·
        end
    end
    theta(pos_num)=theta_ls; %»Ö¸´³öµÄtheta
end