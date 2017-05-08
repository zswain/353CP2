function [m,be] = SOR_1(a,b,tol)
%where a is external sparse, rhs b, and tolerance; m iterations, backward
%error
n=size(a,1);
k=1;
x=zeros(n,1);
xc=zeros(n,1);
w=1.5;
c=zeros([n,1]);
e=zeros([n,1]);
n1=0;
err=0;
relerr=1;
D=diag(diag(a));
L=tril(-a,-1);
U=triu(-a,1);
tw=inv(D-w*L)*((1-w)*D+w*U);
cw=w*inv(D-w*L)*b;
while (relerr>tol)
    x(:,k+1)=tw*x(:,k)+cw;
    err=abs(norm(x(:,k+1)-x(:,k)));
    n1=n1+1;
    k=k+1;
end
xc=x(:,k);
m=n1;
for i=1:n
    for j=1:n
    if (mod(j,2)==0)
        xa(j)=-1;
    else
        xa(j)=1;
    end
    c(i)=c(i)+a(i,j)*xa(j);
    e(i)=e(i)+a(i,j)*xc(j);
    end
end
for i=1:n
    dif(i)=abs(xa(i)-xc(i));
    dif2(i)=abs(e(i)-c(i));
end
fe=max(dif,[],2);
be=max(dif2,[],2);
end

  

