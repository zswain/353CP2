function [m,be] = gauss_seidel(a,b,tol)
%where a is external sparse, rhs b, and tolerance; m iterations, forward
%error, backward error
n=length (b);
p=zeros(n,1);
c=zeros([n,1]);
e=zeros([n,1]);
n1=0;
err=0;
relerr=1;
while (relerr>tol)
    for j=1:n
        if j==n
            x(1)=(b(1)-a(1,2:n)*p(2:n))/a(1,1);
        elseif j==n
            x(n)=(b(n)-a(n,1:n-1)*(x(1:n,1))')/a(n,n);
        else
            x(j)=(b(j)-a(j,1:j-1)*(x(1:j-1))'-a(j,j+1:n)*p(j+1:n,1))/a(j,j);
        end
    end
    err=abs(norm(x'-p));
    relerr=err/(norm(x)+eps);
    p=x';
    n1=n1+1;
end
x=x';
m=n1;
for i=1:n
    for j=1:n
        if (mod(j,2)==0)
            xa(j)=-1;
        else
            xa(j)=1;
        end
        c(i)=c(i)+a(i,j)*xa(j);
        e(i)=e(i)+a(i,j)*x(j);
    end
end
for i=1:n
    dif(i)=abs(xa(i)-x(i));
    dif2(i)=abs(e(i)-c(i));
end
fe=max(dif,[],2);
be=max(dif2,[],2);
end