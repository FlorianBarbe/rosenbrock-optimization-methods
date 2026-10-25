clear;
clc;

n = 10;

%% Variables symboliques
syms b real
x = sym('x', [n 1], 'real');

f=0;
for i=1:(length(x)-1)
    f=f+b*(x(i+1)-x(i)^2)^2+(1-x(i))^2;
end

g=gradient(f,x);
H=hessian(f,x);


disp('Gradient')
disp(g)
disp('Hessienne')
disp(H)

%%Graphes des fonctions pour n=2
b=100; %arbitrairement
x=linspace(-100,100,200);
[X1,X2]=meshgrid(x,x);
F=b*(X2-X1^2)^2 +(1-X1)^2;

figure;

surf(X1, X2, F);

xlabel('x_1');
ylabel('x_2');
zlabel('f(x_1,x_2)');

title('Fonction de Rosenbrock pour n=2');

grid on;
shading interp;

%% Courbes de niveau
figure;

contour(X1, X2, F, 50);

xlabel('x_1');
ylabel('x_2');

title('Courbes de niveau');

grid on;

