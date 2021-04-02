%----------------------------------------------------------------------------------%
%  1-DÐÅºÅÑ¹Ëõ´«¸ÐµÄÊµÏÖ(Õý½»Æ¥Åä×·×Ù·¨Orthogonal Matching Pursuit)   
%  ²âÁ¿ÊýM>=K*log(N/K),KÊÇÏ¡Êè¶È,NÐÅºÅ³¤¶È,¿ÉÒÔ½üºõÍêÈ«ÖØ¹¹
%  ±à³ÌÈË£º ºÎÁõ  Email: aresmiki@163.com
%  ±à³ÌÊ±¼ä£º2017Äê04ÔÂ26ÈÕ  Î÷ÄÏ½»Í¨´óÑ§Ç£Òý¶¯Á¦¹ú¼ÒÖØµãÊµÑéÊÒ
%                                        SWJTU  TPL
%---------------------------------------------------------------------------------%
clc
clear all
close all
%% 1. Éú³ÉÔ­Ê¼ÐÅºÅ
fs=400;     %²ÉÑùÆµÂÊ
f1=25;         %µÚÒ»¸öÐÅºÅÆµÂÊ
f2=50;      %µÚ¶þ¸öÐÅºÅÆµÂÊ
f3=100;     %µÚÈý¸öÐÅºÅÆµÂÊ
f4=200;    %µÚËÄ¸öÐÅºÅÆµÂÊ
N=1024;    %ÐÅºÅ³¤¶È
t=0:1/fs:(N-1)/fs;   
% x=0.3*cos(2*pi*f1*t)+0.6*cos(2*pi*f2*t)+0.1*cos(2*pi*f3*t)+0.9*cos(2*pi*f4*t);  %¹¹ÔìÐÅºÅ
x=cos(2*pi*f1*t)+cos(2*pi*f2*t)+cos(2*pi*f3*t)+cos(2*pi*f4*t);  %¹¹ÔìÐÅºÅ

%% 1.1²é¿´Ê±ÓòºÍ¸µÀïÒ¶Æ×
% fx=abs(fftshift(fft(x)))*2/N;
% fsf=(fs/N)*((1:N)-N/2-1);
% figure
% plot(fsf,fx)
%% ¿ÉÒÔÓÃÏÂÃæ´úÂëÖ±½Ó½«ºÜÐ¡µÄÖµÉèÖÃÎª0
% fft_x=fft(x);
% fft_x(find(abs(fft_x)*2/N<0.1))=0;
% figure
% plot(fsf,fx,fsf,fftshift(fft_x*2/N),'--')
% xx=real(ifft(fft_x));
% figure
% plot(t,x,t,xx,'--')
% x=xx;
%% 2. Ê±ÓòÐÅºÅÑ¹Ëõ´«¸Ð£¬»ñÈ¡²âÁ¿Öµ
K=8;   %ÐÅºÅÏ¡Êè¶È£¬¸µÀïÒ¶Æ×ÖÐ¿´³öÀ´
M=ceil(K*log(N/K));  %²âÁ¿Êý,²âÁ¿¾ØÕóÑ¹Ëõ³Ì¶È£¬¾­Ñé¹«Ê½
randn('state',2)
Phi=randn(M,N);  %  ²âÁ¿¾ØÕó(¸ßË¹·Ö²¼°×ÔëÉù)
Phi=orth(Phi')';    %Õý½»»¯
y=Phi*x';     %  »ñµÃÏßÐÔ²âÁ¿ 

%% 3. L_1·¶Êý×îÓÅ»¯ÖØ¹¹ÐÅºÅ£¨ÓÐ²âÁ¿ÖµyÖØ¹¹x£©
Psi=fft(eye(N,N))/sqrt(N);    %  ¸µÀïÒ¶Õý±ä»»¾ØÕó,ÈÏÎªÐÅºÅxÔÚ¸µÀïÒ¶×ÖµäÉÏÏ¡Êè£ºtheta=Psi*x.  Ôò£ºx=Psi'*theta.
% ×îÐ¡»¯ÎÊÌâ minimize:       ||theta||_0;
%                  subject to:     y=Phi*Psi'*theta;     ==   Áî A=Phi*Psi'.   
A=Phi*Psi';                         %  »Ö¸´¾ØÕó(²âÁ¿¾ØÕó*Õý½»·´±ä»»¾ØÕó);   x=Psi'*theta.

%%  4. Õý½»Æ¥Åä×·×ÙÖØ¹¹ÐÅºÅ
[ fft_y,erro_rn ] = CS_OMP( y,A,2*K );
figure
plot(erro_rn,'-*')
legend('OMPÕý½»Æ¥Åä×·×ÙÎó²î')
r_x=real(Psi'*fft_y');                         %  ×öÄæ¸µÀïÒ¶±ä»»ÖØ¹¹µÃµ½Ê±ÓòÐÅºÅ

%% 5. »Ö¸´ÐÅºÅºÍÔ­Ê¼ÐÅºÅ¶Ô±È

figure;
hold on;
plot(t,r_x,'k.-')                                 %  ÖØ½¨ÐÅºÅ
plot(t,x,'r')                                       %  Ô­Ê¼ÐÅºÅ
xlim([0,t(end)])
legend('OMP»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')


%% 6. CoSaMP ²ÉÓÃÑ¹Ëõ²ÉÑùÆ¥ÅäµÄ·½·¨½âl1×îÐ¡»¯ÎÊÌâ
A=Phi*Psi';                         %  »Ö¸´¾ØÕó(²âÁ¿¾ØÕó*Õý½»·´±ä»»¾ØÕó);   x=Psi'*theta.
[theta,erro_rnn]=CS_CoSaMP( y,A,K );
figure
plot(erro_rnn,'-*')
legend('CoSaMPÑ¹Ëõ²ÉÑù×·×ÙÎó²î')
%% ÖØ¹¹²¢¶Ô±È
cor_x=real(Psi'*theta);                         %  ×öÄæ¸µÀïÒ¶±ä»»ÖØ¹¹µÃµ½Ê±ÓòÐÅºÅ
figure;
hold on;
plot(t,cor_x,'k.-')                                 %  ÖØ½¨ÐÅºÅ
plot(t,x,'r')                                       %  Ô­Ê¼ÐÅºÅ
xlim([0,t(end)])
legend('CoSaMP»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')


%% ISTAºÍFISTAÁ½ÖÖËã·¨½âµÄÄ£ÐÍÎª£º
% minimize   ||x||_1
% subject to:||Ax-y||_2<=eps
% second-order cone program(SOCP)£¬epsÊÇÔëÉùÇ¿¶È  
% ²âÁ¿¾ØÕóyº¬ÓÐÔëÉù£¬Òª¹À¼Æx¾ÍÊÇ½âÉÏÃæµÄÍ¹ÓÅ»¯ÎÊÌâ
% ÎÞÔëÉùÇé¿öÊÇÏÂÃæµÄl1ÎÊÌâ
% minimize ||x||_1
% subject to: Ax=y;
%% 7. ISTA·½·¨½âl1×îÐ¡»¯ÎÊÌâ---ÓÐ×îÓÅ»¯·½·¨×ö
A=Phi*Psi';                         %  »Ö¸´¾ØÕó(²âÁ¿¾ØÕó*Õý½»·´±ä»»¾ØÕó);   x=Psi'*theta.
[theta,erro_rnn]=CS_ISTA( y,A,0.1,2000);
figure
plot(erro_rnn,'-*')
legend('ISTAÎó²î')
%% 8.»Ö¸´ÐÅºÅºÍÔ­Ê¼ÐÅºÅ±È½Ï
figure
plot(t,real(Psi'*theta),'k.-',t,x,'r-')
xlim([0,t(end)])
legend('ISTA»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')

%% 9. ISTA·½·¨½âl1×îÐ¡»¯ÎÊÌâ   ÓÐ×îÓÅ»¯·½·¨×ö
A=Phi*Psi';                         %  »Ö¸´¾ØÕó(²âÁ¿¾ØÕó*Õý½»·´±ä»»¾ØÕó);   x=Psi'*theta.
[theta,erro_rnn]=CS_FISTA( y,A,0.1,2000);
figure
plot(erro_rnn,'-*')
legend('FISTAÎó²î')
%% 10.»Ö¸´ÐÅºÅºÍÔ­Ê¼ÐÅºÅ±È½Ï
figure
plot(t,real(Psi'*theta),'k.-',t,x,'r-')
xlim([0,t(end)])
legend('FISTA»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')