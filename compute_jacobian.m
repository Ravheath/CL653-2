%% dX/dt = F(X,U,D)=[f1;f2;f3], the following function returns the jacobian matrices `     
function [A,B,H] = compute_jacobian(F,X,U,D)
n=size(X);
n=n(1);
m=size(U);
m=m(1);
p=size(D);
p=p(1);
dh=0.000005;
I=eye(n);
for j=1:n  
e=I(j,:)';
A(:,j) = (F(X+dh*e,U,D)-F(X,U,D))/(dh);
end
I=eye(m);
for j=1:m  
e=I(j,:)';
B(:,j) = (F(X,U+dh*e,D)-F(X,U,D))/(dh);
end
I=eye(p);
for j=1:p 
e=I(j,:)';
H(:,j) = (F(X,U,D+dh*e)-F(X,U,D))/(dh);
end

