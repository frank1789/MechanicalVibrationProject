function [xx]=mult(d,x,n)
%Mult please refefer "Mechanical Vibrations - fith edition SI" -Singiresu
%S.Rao / chapter 7.5 Matrix Iteration Method
xx = zeros(length(n));
for i=1:n
    xx(i)=0;
    for j=1:n
        xx(i)=xx(i)+d(i,j)*x(j);
    end
end
