function [freq,eig,b,c,xx]=miter(d,xs,n,nvec,xm,eps)
%Mult please refefer "Mechanical Vibrations - fith edition SI" -Singiresu
%S.Rao / chapter 7.5 Matrix Iteration Method
con=xs(1);
for i=1:n
    x(i)=xs(i)/con;
end
icon=0;
while icon==0
    [xx]=mult(d,x,n);
    alam=xx(1);
    for i=1:n
        xx(i)=xx(i)/alam;
    end
    for i=1:n
        if (abs((xx(i)-x(i))/x(i))>eps)
            icon=0;
        else
            icon=1;
        end
    end
    for i=1:n
        x(i)=xx(i);
    end
end
icon=0;
freq(1)=sqrt(1/alam);
for i=1:n
    eig(i,1)=x(i);
end
ii=2;
while ii<=nvec
    %ii=ii+1;
    sum=0.0;
    for i=1:n
        sum=sum+x(i)^2;
    end
    alp=sqrt(1/sum);
    for i=1:n
        x(i)=x(i)*alp;
    end
    for i=1:n
        for j=1:n
            c(i,j)=x(i)*x(j);
        end
    end
    b=c*xm;
    for i=1:n
        for j=1:n
            d(i,j)=d(i,j)-alam*b(i,j);
        end
    end
    con=xs(1);
    for i=1:n
        x(i)=xs(i)/con;
    end
    icon=0;
    while icon==0
        [xx]=mult(d,x,n);
        alam=xx(1);
        for i=1:n
            xx(i)=xx(i)/alam;
        end
        for i=1:n
            if (abs((xx(i)-x(i))/x(i))>eps)
                icon=0;
            else
                icon=1;
            end
        end
        for i=1:n
            x(i)=xx(i);
        end
    end
    icon=0;
    freq(ii)=sqrt(1/alam);
    for i=1:n
        eig(i,ii)=x(i);
    end
    ii=ii+1;
end