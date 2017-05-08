function [m,fe,be] = jacobi_2_que(a,b,tol)
%where a is sparse from other fn, rhs b and tolerance; m max iterations,
%forward error, backward error
n=length(b);
d=diag(a);
r=a-diag(d);
xc=zeros(n,1);
p=zeros(n,1);
c=zeros([n,1]);
e=zeros([n,1]);
n1=0;
err=0;
relerr=1;
while (relerr>tol)
    xc=(b-r*xc)./d;
    err=abs(norm(xc-p,inf));
    relerr=err/(norm(xc)+eps);
    p=xc;
    n1=n1+1;
end
x=xc;
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