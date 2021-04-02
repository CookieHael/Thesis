%% CS_RL1_test
%------------------------------------------------------------------------------------------%
% The weighted L1 minimization can be viewed as a relaxation of a weighted L0 minimization problem
%  minimize W||x||_0
%  subject to Ax=y
%  ±à³ÌÈË£º ºÎÁõ                                    Email: aresmiki@163.com
%  ±à³ÌÊ±¼ä£º2017Äê05ÔÂ3ÈÕ  Î÷ÄÏ½»Í¨´óÑ§Ç£Òý¶¯Á¦¹ú¼ÒÖØµãÊµÑéÊÒ
%                                        SWJTU  TPL
%------------------------------------------------------------------------------------------%
clc;clear all;close all
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
% M=2*ceil(K*log(N/K));
M=152;
Phi=randn(M,N);  %¸ßË¹¾ØÕó×÷Îª¸ÐÖª¾ØÕó
Phi=orth(Phi')';  %Õý½»»¯

%% 3. Á¿²âÐÅºÅ
y=Phi*x;
A=Phi;%»Ö¸´¾ØÕó,Ï¡Êè»¯¾ØÕóÎªµ¥Î»¾ØÕó£¬ÒòÎªÐÅºÅ±¾Éí¾ÍÊÇÏ¡ÊèµÄ£¬²»ÐèÒª×öÈÎºÎÏ¡Êè±ä»»
%% CS_RL1    ¼ÓÈ¨L1×îÐ¡»¯£¬Ïàµ±ÓÚL0×îÐ¡»¯
[theta]=CS_RL1( y,A,1);
figure
subplot(3,1,1)
scatter(x,theta)
hold on
plot([min(x):0.01:max(x)],[min(x):0.01:max(x)],'r')
hold off
[theta]=CS_RL1( y,A,2);
subplot(3,1,2)
scatter(x,theta)
hold on
plot([min(x):0.01:max(x)],[min(x):0.01:max(x)],'r')
hold off

[theta]=CS_RL1( y,A,5);
subplot(3,1,3)
scatter(x,theta)
hold on
plot([min(x):0.01:max(x)],[min(x):0.01:max(x)],'r')
hold off
figure
plot(t,theta,'ko',t,x,'r.')
xlim([0,t(end)])
legend('RL1»Ö¸´ÐÅºÅ','Ô­Ê¼ÐÅºÅ')