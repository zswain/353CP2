function [a,b] = sparsesetup1(n)
%where n is system size, a and b are sparse matrices
e=ones(n,1);
a=spdiags([e 2*e e],-1:1 ,n, n);
b=zeros(n,1);
b(1)=1;
b(n)=-1;
b(2:n-1)=0;
end