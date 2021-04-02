

%%
%----------------------------------------------------------------------------------%
%  1-DÐÅºÅÑ¹Ëõ´«¸ÐµÄÊµÏÖ(l1-MAGICºÍl1_ls½âl1ÎÊÌâ)     ÐÅºÅ±¾Éí¾ÍÊÇÏ¡ÊèµÄ£¬
%  ²»ÐèÒªÏ¡Êè¾ØÕó£¬»Ö¸´¾ØÕóAÊÇµ¥Î»Õý½»¾ØÕó£¬ÓÃOMP·½·¨Çó½âl1ÎÊÌâ.
%  ²âÁ¿ÊýM>=K*log(N/K),KÊÇÏ¡Êè¶È,NÐÅºÅ³¤¶È,¿ÉÒÔ½üºõÍêÈ«ÖØ¹¹
%  ±à³ÌÈË--Î÷ÄÏ½»Í¨´óÑ§Ç£Òý¶¯Á¦¹ú¼ÒÖØµãÊµÑéÊÒ ºÎÁõ  Email: aresmiki@163.com
%  ±à³ÌÊ±¼ä£º2017Äê04ÔÂ27ÈÕ
%---------------------------------------------------------------------------------%
clc;clear all;close all;
%% 1. ²úÉúÏ¡ÊèµÄÐÅºÅ
N=1024;
K=50;
x=zeros(N,1);
rand('state',8)
q=randperm(N); %Ëæ»úÅÅÁÐ1µ½NµÄÕûÊý
randn('state',10)
x(q(1:K))=randn(K,1); %½«K¸öËæ»úÊýËæ»ú·Åµ½xÖÐ
t=0:N-1;
%% 2. ¹¹Ôì¸ÐÖª¾ØÕó
M=2*ceil(K*log(N/K));
Phi=randn(M,N);  %¸ßË¹¾ØÕó×÷Îª¸ÐÖª¾ØÕó
Phi=orth(Phi')';  %Õý½»»¯

%% 3. Á¿²âÐÅºÅ
y=Phi*x;
A=Phi;%»Ö¸´¾ØÕó,Ï¡Êè»¯¾ØÕóÎªµ¥Î»¾ØÕó£¬ÒòÎªÐÅºÅ±¾Éí¾ÍÊÇÏ¡ÊèµÄ£¬²»ÐèÒª×öÈÎºÎÏ¡Êè±ä»»
%% 4. ÖØ¹¹ÐÅºÅ l1×îÐ¡»¯   Using  l1-MAGIC  ¸Ã·½·¨»Ö¸´½ÏºÃ Ä£ÐÍÆ¥Åämin_x ||x||_1  s.t.  Ax = b
x0=A'*y;  %×îÐ¡¶þ³Ë½â¹À¼ÆÒ»¸ö³õÊ¼Öµ
xh1=l1eq_pd(x0,A,[],y,1e-3);

%% 5. ÖØ¹¹ÐÅºÅl1×îÐ¡»¯   Using l1_ls ¸Ã·½·¨ÉÔÎ¢²îµã Ä£ÐÍÎª minimize ||A*x-y||^2 + lambda*sum|x_i|,
lambda  = 0.01; % ÕýÔò»¯²ÎÊý
rel_tol = 1e-3; % Ä¿±êÏà¶Ô¶ÔÅ¼¼äÏ¶
quiet=1;   %²»Êä³öÖÐ¼ä½á¹û
[xh2,status]=l1_ls(A,y,lambda,rel_tol,quiet);

%% 6. »Ö¸´ÐÅºÅºÍÔ­Ê¼ÐÅºÅ±È½Ï 
figure
plot(t,xh1,'ko',t,x,'r.')
xlim([0,t(end)])
legend('l1-MAGIC»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')

figure
plot(t,xh2,'ko',t,x,'r.')
xlim([0,t(end)])
legend('l1-ls»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')

%% 7. ÓÃÕý½»Æ¥Åä×·×ÙµÄ·½·¨¼ÆËãl1×îÓÅ»¯ÎÊÌâ
[ xh,erro_rn ] = CS_OMP( y,A,2*K );
figure
plot(erro_rn,'-*')
legend('OMPÕý½»Æ¥Åä×·×ÙÎó²î')
figure
plot(t,xh,'ko',t,x,'r.')
xlim([0,t(end)])
legend('OMP»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')

%% 8. ÓÃÑ¹Ëõ²ÉÑùÆ¥Åä×·×Ù(CoSaMP)µÄ·½·¨¼ÆËãl1×îÓÅ»¯ÎÊÌâ
%Needell D, Tropp J A. CoSaMP: Iterative signal recovery from incomplete 
% and inaccurate samples [J]. Applied & Computational Harmonic Analysis, 2008, 26(3):301-321.
% Ò»´ÎÐÔÑ¡Ôñ2*K¸ö½Ï´óµÄ»ù£¬Ã¿´ÎÑ­»·²»¶ÏÉ¾³ýºÍ²¹½øÒ»¶¨ÊýÄ¿µÄ»ù£¬×î¿ì´ïµ½¸ø¶¨¾«¶È

A=Phi;    %»Ö¸´¾ØÕó

%% CoSaMP
 [theta,erro_rnn]=CS_CoSaMP( y,A,K );
 figure
plot(erro_rnn,'-*')
legend('CoSaMPÑ¹Ëõ²ÉÑù×·×ÙÎó²î')
 %%    
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('CoSaMP»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')

%% CS_NSRAL0
A=Phi;    %»Ö¸´¾ØÕó 
deltaT=1e-3;
r=1/3;
te=0.01;
eps=0.09;
[theta,Spare_L0]=CS_NSRAL0(y,A,deltaT,r,te,eps);
figure
plot(1:length(Spare_L0),Spare_L0,'b*-',1:length(Spare_L0),ones(length(Spare_L0),1)*K,'r^-')
legend('½üËÆL0·¶Êý×·×Ù½á¹û','ÕæÊµÏ¡Êè¶È')
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('NSRAL0»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')

%% CS_SL0
A=Phi;    %»Ö¸´¾ØÕó 
[theta,Spare_L0]=CS_SL0( y,A,0.001,0.9,2,3);
figure
plot(1:length(Spare_L0),Spare_L0,'b*-',1:length(Spare_L0),ones(length(Spare_L0),1)*K,'r^-')
legend('½üËÆL0·¶Êý×·×Ù½á¹û','ÕæÊµÏ¡Êè¶È')
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('SL0»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')

%% CS_UALP
A=Phi;    %»Ö¸´¾ØÕó 
[theta]=CS_UALP( y,A,0.1);
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('UALP»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')

%% CS_RSL0
A=Phi;    %»Ö¸´¾ØÕó 
[theta,Spare_L0]=CS_RSL0( y,A,0.001,0.9,2,3,0.001);
figure
plot(1:length(Spare_L0),Spare_L0,'b*-',1:length(Spare_L0),ones(length(Spare_L0),1)*K,'r^-')
legend('½üËÆL0·¶Êý×·×Ù½á¹û','ÕæÊµÏ¡Êè¶È')
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('RSL0»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')


%% CS_IRLS
A=Phi;    %»Ö¸´¾ØÕó 
[theta]=CS_IRLS( y,A,0,1e-8,0.1);
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('IRLS»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')

%% CS_RL1    ¼ÓÈ¨L1×îÐ¡»¯£¬Ïàµ±ÓÚL0×îÐ¡»¯
A=Phi;    %»Ö¸´¾ØÕó 

[theta]=CS_RL1( y,A,3);
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('RL1»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')

%% CS_IHT    ¼ÓÈ¨L1×îÐ¡»¯£¬Ïàµ±ÓÚL0×îÐ¡»¯
A=Phi;    %»Ö¸´¾ØÕó 

% [theta]=CS_IHT( y,A,K);
[theta]=CS_IHT( y,A);
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('IHT»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')


%% CS_SBI    
A=Phi;    %»Ö¸´¾ØÕó 

theta=CS_SBIL1(y,A,2000,100,1e4)
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('SBIL1»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')

%% ISTA
A=Phi;    %»Ö¸´¾ØÕó  
[theta,erro_rnn]=CS_ISTA( y,A,0.00819); %0.00819
figure
plot(erro_rnn,'-*')
legend('ISTAÎó²î')
%% »Ö¸´ÐÅºÅºÍÔ­Ê¼ÐÅºÅ±È½Ï
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('ISTA»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')

%% FISTA
A=Phi;    %»Ö¸´¾ØÕó  
[theta,erro_rnn]=CS_ISTA( y,A,0.00819); %0.00819
figure
plot(erro_rnn,'-*')
legend('FISTAÎó²î')
%% »Ö¸´ÐÅºÅºÍÔ­Ê¼ÐÅºÅ±È½Ï
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('FISTA»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')





