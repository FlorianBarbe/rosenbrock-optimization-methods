function [f,g,H,r,J]=rosenbrock(x,params)
%ROSENBROCK Fonction de Rosenbrock, derivees et forme moindres carres.
%   La convention utilisee est f(x)=r(x)'*r(x), donc g(x)=2*J(x)'*r(x).

    x=x(:);
    n=length(x);
    b=params.b;

    if n < 2
        error('rosenbrock:Dimension','Le vecteur x doit avoir au moins deux composantes.');
    end
    if ~isscalar(b) || ~isfinite(b) || b <= 0
        error('rosenbrock:Parameter','Le parametre b doit etre un scalaire strictement positif.');
    end

    %% Calcul de f(x)
    f=0;
    for i=1:n-1
        f=f+b*(x(i+1)-x(i)^2)^2+(1-x(i))^2;
    end

    %% Calcul de g(x)
    g=zeros(n,1);
    g(1)=-4*b*x(1)*(x(2)-x(1)^2)-2*(1-x(1));
    g(n)=2*b*(x(n)-x(n-1)^2);

    for i=2:n-1
        g(i)=2*b*(x(i)-x(i-1)^2) ...
            -4*b*x(i)*(x(i+1)-x(i)^2)-2*(1-x(i));
    end

    %% Calcul de H(x)
    H=zeros(n,n);
    H(1,1)=2-4*b*(x(2)-3*x(1)^2);
    H(n,n)=2*b;

    for i=1:n-1
        H(i,i+1)=-4*b*x(i);
        H(i+1,i)=H(i,i+1);
    end

    for i=2:n-1
        H(i,i)=2+2*b-4*b*(x(i+1)-3*x(i)^2);
    end

    %% Calcul de r(x) et J(x)
    r=zeros(2*n-2,1);
    J=zeros(2*n-2,n);

    for i=1:n-1
        r(2*i-1)=sqrt(b)*(x(i+1)-x(i)^2);
        r(2*i)=1-x(i);

        J(2*i-1,i)=-2*sqrt(b)*x(i);
        J(2*i-1,i+1)=sqrt(b);
        J(2*i,i)=-1;
    end
end
