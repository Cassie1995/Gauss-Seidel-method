function y= Gauss_Seidel(A,F,x0,ep )
D=diag(diag(A));
L=tril(A,-1);
U=triu(A,1);
BS=-(D+L)\U;        
g=(D+L)\F;
y=BS*x0+g;           
n=1;
while norm(y-x0,inf)>=ep
    x0=y;
    y=BS*x0+g;
    n=n+1;
end

disp('Number of Iterations');
disp(n);
end

