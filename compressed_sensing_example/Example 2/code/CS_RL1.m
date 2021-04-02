%% CS_RL1 Reweighted L1 Minimization Algorithm
% ÊäÈë£ºy---²âÁ¿ÐÅºÅ M X 1
%          A---»Ö¸´¾ØÕó M X N
%          iter---×î´óµü´ú´ÎÊý  ÖÁÉÙ´óÓÚ2£¬Ð¡ÓÚ2¾ÍÏàµ±ÓÚL1×îÐ¡»¯£¨Ã»ÓÐ¼ÓÈ¨£©

% Êä³ö£ºys---»Ö¸´µÄÐÅºÅ N X 1
%    
% 
%  minimize W||x||_1
%  subject to Ax=y
%  ±à³ÌÈË£º ºÎÁõ                                    Email: aresmiki@163.com
%  ±à³ÌÊ±¼ä£º2017Äê04ÔÂ30ÈÕ  Î÷ÄÏ½»Í¨´óÑ§Ç£Òý¶¯Á¦¹ú¼ÒÖØµãÊµÑéÊÒ
%                                        SWJTU  TPL
%  ²Î¿¼ÎÄÏ×£º Cand¨¨s E J, Wakin M B, Boyd S P. 
% Enhancing sparsity by reweighted L1 minimization.[J]. 
% Journal of Fourier Analysis & Applications, 2007, 14(5):877-905.
%---------------------------------------------------------------------------------------------------------------------%
function xh=CS_RL1(y,A,iter)
N=max(size(A));
M=min(size(A));
Li=(M/(4*log(N/M)));
y=y(:);
W=ones(N,1);   %³õÊ¼»¯È¨ÖØÏòÁ¿
QW=diag(W);    %³õÊ¼»¯È¨ÖØ¾ØÕó
% delta=0.01;
for i=1:iter
    QWt=inv(QW);
    At=A*QWt;
    x0=At'*y;  %×îÐ¡¶þ³Ë½â¹À¼ÆÒ»¸ö³õÊ¼Öµ
    xh=l1eq_pd(x0,At,[],y,1e-3);
    delta=max(norm(xh,Li),1e-3) ;%¶¯Ì¬¸üÐÂ ×ÔÊÊÓ¦¸üÐÂdelta
%     delta should be set slightly smaller than the expected nonzero magnitudes of x0. 
%     In general, the recovery process tends to be reasonably robust to the choice of delta.--Ô­ÎÄÖÐµÄ»°
    xh=QWt*xh;
    QW=diag(1./(abs(xh)+delta));
end

end
