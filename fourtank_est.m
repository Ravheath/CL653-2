clear all
phi=[0.9233,0,0.1813,0;
    0,0.9462,0,0.1493;
    0,0,0.8112,0;
    0,0,0,0.8465];
tau=[0.4001,0.02276;
    0.01209,0.3055;
    0,0.2159;
    0.1438,0];
C=[0.5,0,0,0;
    0,0.5,0,0];
xo= [0.4;0.4;0.1;0.2];
u=[0.6;0.6]; %constant input
R=[0.01,0;
    0,0.02];
randn('state',0);
N=50;
V= mvnrnd([0,0],R,N);
V=V';
Y=[];
X_true=xo;
for i=1:N
    xk= phi*(X_true(:,end))+tau*u;
    yk=C*X_true(:,end)+V(:,i);
    X_true =[X_true,xk];
    Y=[Y,yk];
end
Z=Y(:,1);
psi=C;
N1=10; % no. of measurements used in batch estimation
for i=1:(N1-1)
    psi=[psi;C*(phi^i)];
end
for i=1:N1-1
    zk=Y(:,i+1);
    for j=1:i
        zk=zk-C*(phi^(j-1))*tau*u;
    end
    Z=[Z;zk];
end
Wr=[];
Rinv=R^-1;
for i=1:N1
    Wr=blkdiag(Wr,Rinv);
end
%% initial state estimation
xo_est= ((psi'*Wr*psi)^-1)*psi'*Wr*Z;
Xest=xo_est;
%% estimation of successive states
for i=1:N1
    xk= phi*(Xest(:,end))+tau*u;
    Xest=[Xest,xk];
end
% plot(Xest(1,:),'r');
% hold on;
% plot(X_true(1,1:N1+1),'b');
% figure();
% plot(Xest(2,:),'r');
% hold on;
% plot(X_true(2,1:N1+1),'b');
% figure();
% plot(Xest(3,:),'r');
% hold on;
% plot(X_true(3,1:N1+1),'b');
% figure();
% plot(Xest(4,:),'r');
% hold on;
% plot(X_true(4,1:N1+1),'b');

diff=xo-xo_est;

%% recursive estimation



        
    
    
    
