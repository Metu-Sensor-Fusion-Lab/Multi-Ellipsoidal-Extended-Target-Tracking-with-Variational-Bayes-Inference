function h=drawEllipse(x,P,c,ax,linewidth,color)

[nx,N]=size(x);
if nx~=2
    error('nx must be 2');
end
phi=linspace(0,2*pi,100);
z=[cos(phi);sin(phi)];

h=zeros(1,N);
htext=zeros(1,N);
if (N==1)
    [U,S,V]=svd(P);
    Sroot=sqrt(S);
    Ph=U*Sroot;
    h=plot(x(1)+c*Ph(1,:)*z,x(2)+c*Ph(2,:)*z,['-'],'Color',color,'parent',ax,'linewidth',linewidth);
end