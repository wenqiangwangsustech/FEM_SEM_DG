function[I,x,y,w]=QuadLTri(func,n)
% function[I,x,y,w]=QuadLTri(func,n)
% Call [I]=QuadLTri(n) : displays nodes on triangle
% Call [I,x,y,z]=QuadLTri(n) : returns nodes (x,y) and weights (w)
% Call [I]=QuadLTri(func,n) : evaluates numerical integral over triangle
% The triangle is defined by the nodes (0,0), (1,0), and (0,1)
% Inputs:
% func: function handle of two arguments f(x,y)
% n: integer (n>0), the order of the Gauss-Legendre quadrature
% Outputs:
% I: the integral reuslt or the quadrature order n
% x: nodes on alpha axis
% y: nodes on beta axis
% w: weights
%% Main Function
switch nargin
    case 1
        if func<1
            error("Order n must be greater than zero.");
        end
        n       =   func;
        %showNodes(round(n));
        [x,w]	=	GaussLegendre(round(n));
        [zeta,eta]=meshgrid(x,x);
        x       =   (3+3*zeta-eta-zeta.*eta)/8;
        y       =   (3-zeta+3*eta-zeta.*eta)/8;
        I       =   func;
        w    	=   w.*w'.*(2-zeta-eta)/16;
        x       =   reshape(x,1,[]);
        y       =   reshape(y,1,[]);
        w       =   reshape(w,1,[]);
    case 2
        if n<1
            error("Order must be greater than zero.");
        end
        [x,w]	=	GaussLegendre(round(n));
        [zeta,eta]=meshgrid(x,x);
        x       =   (3+3*zeta-eta-zeta.*eta)/8;
        y       =   (3-zeta+3*eta-zeta.*eta)/8;
        J     	=   (2-zeta-eta)/16;
        I       =   w'*(func(x,y).*J)*w;
        w    	=   w.*w'.*(2-zeta-eta)/16;
        x       =   reshape(x,1,[]);
        y       =   reshape(y,1,[]);
        w       =   reshape(w,1,[]);
    otherwise
        error("Not enough input arguments.");
end
end
%%
function[]=showNodes(n)
%% Gauss-Legendre
[x,~]	=	GaussLegendre(n);
[zeta,eta] = meshgrid(x,x);
x       =   (3+3*zeta-eta-zeta.*eta)/8;
y       =   (3-zeta+3*eta-zeta.*eta)/8;
%% Show Nodes
figure()
hold on
fill([0 1 0.5 0],[0 0 0.5 1],[0.9 0.9 0.9],'LineWidth',1)
plot(0,0,'.b','MarkerSize',15)
plot(1,0,'.b','MarkerSize',15)
plot(0.5,0.5,'.b','MarkerSize',15)
plot(0,1,'.b','MarkerSize',15)
plot(x,y,'.k','MarkerSize',10)
hold off
xlabel('$\alpha$','Interpret','Latex','FontSize',15)
ylabel('$\beta$','Interpret','Latex','FontSize',15)
set(gca,'TickLabel','Latex','FontSize',15)
axis equal
xlim([0 1])
ylim([0 1])
str     =   sprintf('$n=%d$',n);
title(str,'Interpret','Latex','FontSize',15)
%% 
end
%% Gauss-Legendre Nodes and Weights
function [x,w]=GaussLegendre(n)
% © Geert Van Damme
% geert@vandamme-iliano.be
% February 21, 2010  
i       =   1:n-1;
a       =   i./sqrt(4*i.^2-1);
CM      =   diag(a,1) + diag(a,-1);
[V,L]   =   eig(CM);
[x,ind] =   sort(diag(L));
V       =   V(:,ind)';
w       =   2*V(:,1).^2;
end