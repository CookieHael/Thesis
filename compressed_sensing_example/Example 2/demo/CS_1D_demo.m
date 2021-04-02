%% º¬ÓÐÔëÉù
% minimize ||x||_1
% subject to: (||Ax-y||_2)^2<=eps;
% minimize :  (||Ax-y||_2)^2+lambda*||x||_1
% y´«ÊäÖÐ¿ÉÄÜº¬Ôë y=y+w
%
%%
clc;clear all;
%% 1.¹¹ÔìÒ»¸öÁ½¸öÐ³²¨ÐÅºÅ
lam=0.37;
itrs=400;
m=380;
sig=0.5;
n=1024;
dt=1/2000;
T=1023*dt;
t=0:dt:T;
t=t(:);
x=sin(697*pi*t)+sin(1975*pi*t);
Dn=dctmtx(n);

%% 2.¹¹Ôì²âÁ¿¾ØÕó 
rand('state',15);
q=randperm(n);
q=q(:);
y=x(q(1:m));
randn('state',7)
w=sig*randn(m,1);  %²úÉúÔëÉù
yn=y+w;  %Ñ¹Ëõ¾ØÕóÓÐÔëÉù
Psi1=Dn';
%% 3. ÖØ¹¹ÐÅºÅ  ISTA
A=Psi1(q(1:m),:);
% [xh,err]=CS_ISTA(yn,A,lam,itrs);  %ISTA
[xh,err]=CS_ISTA(yn,A,lam);  %ISTA
xx=Psi1*xh;
figure
plot(err,'*-')
legend('ISTAÎó²î')

figure
plot(t,x,'b',t,xx,'r');
legend('DCT-Ï¡ÊèÐÅºÅ','ISTAÖØ¹¹ÐÅºÅ')

figure
t1=50*dt:dt:100*dt;
plot(t1,x(50:100),'b',t1,xx(50:100),'r','linewidth',1.5)
legend('DCT-Ï¡ÊèÐÅºÅ','ISTAÖØ¹¹ÐÅºÅ')

%% 3. ÖØ¹¹ÐÅºÅ  FISTA
A=Psi1(q(1:m),:);
% [xh,err]=CS_FISTA(yn,A,lam,itrs);  %FISTA
[xh,errr]=CS_FISTA(yn,A,lam);  %FISTA
xx=Psi1*xh;
figure
plot(errr,'*-')
legend('FISTAÎó²î')

figure
plot(t,x,'b',t,xx,'r');
legend('DCT-Ï¡ÊèÐÅºÅ','FISTAÖØ¹¹ÐÅºÅ')

figure
t1=50*dt:dt:100*dt;
plot(t1,x(50:100),'b',t1,xx(50:100),'r','linewidth',1.5)
legend('DCT-Ï¡ÊèÐÅºÅ','FISTAÖØ¹¹ÐÅºÅ')

%% 4. ÖØ¹¹ÐÅºÅ  OMP
A=Psi1(q(1:m),:);
[xh,errr]=CS_OMP(yn,A,100);  %OMP
xx=Psi1*xh';
figure
plot(errr,'*-')
legend('OMPÎó²î')

figure
plot(t,x,'b',t,xx,'r');
legend('DCT-Ï¡ÊèÐÅºÅ','OMPÖØ¹¹ÐÅºÅ')

figure
t1=50*dt:dt:100*dt;
plot(t1,x(50:100),'b',t1,xx(50:100),'r','linewidth',1.5)
legend('DCT-Ï¡ÊèÐÅºÅ','OMPÖØ¹¹ÐÅºÅ')

%% 4. ÖØ¹¹ÐÅºÅ  CoSaMP
A=Psi1(q(1:m),:);
[xh,errr]=CS_CoSaMP(yn,A,100);  %CoSaMP
xx=Psi1*xh;
figure
plot(errr,'*-')
legend('CoSaMPÎó²î')

figure
plot(t,x,'b',t,xx,'r');
legend('DCT-Ï¡ÊèÐÅºÅ','CoSaMPÖØ¹¹ÐÅºÅ')

figure
t1=50*dt:dt:100*dt;
plot(t1,x(50:100),'b',t1,xx(50:100),'r','linewidth',1.5)
legend('DCT-Ï¡ÊèÐÅºÅ','CoSaMPÖØ¹¹ÐÅºÅ')

%% 5. l1ÎÊÌâ×îÐ¡»¯ l1-Magic¹¤¾ßÏä  Õâ¸öÉÔ²î£¬Ä£ÐÍÎªmin_x ||x||_1  s.t.  Ax = b£¬±¾Éí²»ÊÊºÏ¸ÃÄ£ÐÍ
A=Psi1(q(1:m),:); 
% x0=A'*y;   %×îÐ¡¶þ³Ë×÷Îªl1×îÐ¡»¯µÄ³õÊ¼Öµ¹À¼Æ
% ÓÃl1-MAGICµÄMATLAB¹¤¾ßÏä½âl1×îÐ¡»¯ÎÊÌâ
% xh1=l1eq_pd(x0,A,[],y,1e-3);
xh=l1eq_pd(zeros(n,1),A,[],y,1e-3);  %¿ÉÒÔ²»¸ø³õÊ¼µÄ¹À¼Æ
xx=Psi1*xh;
figure
plot(t,x,'b',t,xx,'r');
legend('DCT-Ï¡ÊèÐÅºÅ','l1-MagicÖØ¹¹ÐÅºÅ')

figure
t1=50*dt:dt:100*dt;
plot(t1,x(50:100),'b',t1,xx(50:100),'r','linewidth',1.5)
legend('DCT-Ï¡ÊèÐÅºÅ','l1-MagicÖØ¹¹ÐÅºÅ')

%% 5. l1ÎÊÌâ×îÐ¡»¯ l1-ls¹¤¾ßÏä   Õâ¸öÄ£ÐÍminimize ||A*x-y||^2 + lambda*sum|x_i|,·Ç³£·ûºÏ¸ÃÎÊÌâ
A=Psi1(q(1:m),:); 
lambda  = 0.01; % ÕýÔò»¯²ÎÊý
rel_tol = 1e-3; % Ä¿±êÏà¶Ô¶ÔÅ¼¼äÏ¶
quiet=1;   %²»Êä³öÖÐ¼ä½á¹û
[xh,status]=l1_ls(A,y,lambda,rel_tol,quiet);
xx=Psi1*xh;
figure
plot(t,x,'b',t,xx,'r');
legend('DCT-Ï¡ÊèÐÅºÅ','l1-lsÖØ¹¹ÐÅºÅ')

figure
t1=50*dt:dt:100*dt;
plot(t1,x(50:100),'b',t1,xx(50:100),'r','linewidth',1.5)
legend('DCT-Ï¡ÊèÐÅºÅ','l1-lsÖØ¹¹ÐÅºÅ')



%% 6. SL0   Õâ¸öÄ£ÐÍminimize £¨A*x-y£© + lambda*sum|x_i|
A=Psi1(q(1:m),:); 

[xh,Spare_L0]=CS_SL0( y,A,0.001,0.9,2,3);

xx=Psi1*xh;
figure
plot(t,x,'b',t,xx,'r');
legend('DCT-Ï¡ÊèÐÅºÅ','SL0ÖØ¹¹ÐÅºÅ')

figure
t1=50*dt:dt:100*dt;
plot(t1,x(50:100),'b',t1,xx(50:100),'r','linewidth',1.5)
legend('DCT-Ï¡ÊèÐÅºÅ','SL0ÖØ¹¹ÐÅºÅ')

%% 6. SL0   Õâ¸öÄ£ÐÍminimize ||A*x-y||_2 + lambda*sum|x_i|
A=Psi1(q(1:m),:); 

[xh,Spare_L0]=CS_RSL0( y,A,0.001,0.9,2,3,0.01);

xx2=Psi1*xh;
figure
plot(t,x,'b',t,xx2,'r');
legend('DCT-Ï¡ÊèÐÅºÅ','RSL0ÖØ¹¹ÐÅºÅ')

figure
t1=50*dt:dt:100*dt;
plot(t1,x(50:100),'b',t1,xx2(50:100),'r','linewidth',1.5)
legend('DCT-Ï¡ÊèÐÅºÅ','RSL0ÖØ¹¹ÐÅºÅ')

figure
plot(t,x,'b',t,xx,'r',t,xx2,'k');
legend('DCT-Ï¡ÊèÐÅºÅ','SL0ÖØ¹¹ÐÅºÅ','RSL0ÖØ¹¹ÐÅºÅ')

figure
bar([norm(x-xx),norm(x-xx2)],'r')



