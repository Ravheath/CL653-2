clear all;
N=100; %sample size
Q=diag([0.01,0.01,0.001,0.001]); %covariance matrix for process noise 
rand('state',0);
W=mvnrnd([0,0,0,0],Q,N-1);
W=W';
T=5; %sampling time
Xs = [12.4,12.7,1.8,1.4]'; %Steady state
U=zeros(2,N-1); %Input matrix
pert=[1,-1,1,-1]'; %initial perturbation
sol = Xs+pert; %initial solution
 C=[0.5,0,0,0;
    0,0.5,0,0]; % C from y=C*x+v
R=diag([0.01,0.01]); %covariance matrix for measurement noise
Y=zeros(2,N); 
rand('state',0);
V=mvnrnd([0,0],R,N);
V=V';
Y(:,1)=C*sol+V(:,1);
% Input generation
for i=1:N-1
    U(:,i)=[3;3]+[0.1*sin(0.025*(i-1))+0.15*cos(0.02*(i-1));
              0.15*sin(0.02*(i-1))-0.1*cos(0.025*(i-1))];
end
% True state and measurement generation by integrating non-linear differential equation
for i=1:N-1 
    [tsol,hsol] = ode45(@(t,h) fourtank1(t,h,U(:,i),W(:,i)),[(i-1)*T, i*T], sol(:,end));
    sol= [sol,hsol(end,:)'];
    Y(:,i+1)=C*(sol(:,end))+V(:,i+1);
end
%% EKF impementation
X0=[100,5,6.9,5]'; % initial guess
P0=diag([1000,2000,3000,1000]); % covariance of initial guess; arbitrarily chosen very high values
%using the measurement at t=0 to improve the initial guess of X0
Lk=P0*(C')*((C*P0*(C')+R)^-1); 
ek=Y(:,1)-C*X0;
X0=X0+Lk*ek;
P0=[eye(4)-Lk*C]*P0;
% estimating further states
X_est=zeros(4,N);
X_est(:,1)=X0;
P=zeros(4,4*N);
P(:,1:4)=P0;
for i=1:N-1 
    % Prediction
    [time,hieght]=ode45(@(t,h) fourtank1(t,h,U(:,i),[0;0;0;0]),[0, T], X_est(:,i));
    xk=hieght(end,:)';
    Pk= P(:,4*i-3:4*i); %extracting the current covariance matrix 
    %covariance prediction
    f=@fourtank;
    [A,B,H]=compute_jacobian(f,X_est(:,i),U(:,i),[0,0,0,0]');
    sys_c=ss(A,H,[],[]);
    sys_d=c2d(sys_c,T);
    phi=sys_d.a;
    tau=sys_d.b;
    Pk_pred=phi*Pk*(phi')+tau*Q*(tau');
    %kalman gain computation
    Lk=Pk_pred*(C')*((C*Pk_pred*(C')+R)^-1);
    %update step
    ek=Y(:,i+1)-C*xk;
    xk_upd=xk+Lk*ek;
    Pk_upd=(eye(4)-Lk*C)*Pk_pred;
    if (xk_upd(1))<0
       xk_upd(1)=0;
    end   
    if (xk_upd(2))<0
       xk_upd(2)=0;
    end
    if (xk_upd(3))<0
       xk_upd(3)=0;
    end
    if (xk_upd(4))<0
       xk_upd(4)=0;
    end     
    %storage
    P(:,4*i+1:4*(i+1))=Pk_upd;
    X_est(:,i+1)=xk_upd;
end
figure();
plot(sol(1,:),'r');
hold on;
plot(X_est(1,:),'b');
figure();
plot(sol(2,:),'r');
hold on;
plot(X_est(2,:),'b');
figure();
plot(sol(3,:),'r');
hold on;
plot(X_est(3,:),'b');
figure();
plot(sol(4,:),'r');
hold on;
plot(X_est(4,:),'b');