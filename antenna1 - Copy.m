function beem= antenna1 (R,numb,P)
%cercle
theta = linspace(0, 2*pi);
rc =R *cos(theta);
xc=rc.*cos(theta);
yc=rc.*sin(theta);  
figure(2)
plot(xc,yc);
point=length (theta)/numb;
beem=point*2*(P/100+1)*3.6;
end