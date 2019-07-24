clc
clf
clear all


n=0.5;
B0=10;
num=18;
P=0;
R=16;
resolution = 10000;


% %cercle
 theta = linspace(-pi/2, pi/2);
rc =R *cos(theta);
xc=rc.*cos(theta);
yc=rc.*sin(theta);  
theta = linspace(-pi/2, pi/2);

%plot the vaiation of beamwdth vs n
for i=1:20
%radiation intensity polar
U =B0 *cos(theta).^i;
figure(1)
polarplot(theta,U)
hold on
grid on

%radiation inensity cartesian
x=U.*cos(theta);
y=U.*sin(theta); 
figure(2)
plot(x,y)
grid on
hold on
end

%rotation of the radiation pattern

U =B0 *cos(theta).^n;
x=U.*cos(theta);
y=U.*sin(theta);
figure(7)
polarplot(theta,U)

%beamwidth calculation
x3Db=max(x)/2;
[value,row]= min(abs(x-x3Db));
y3Db=abs(y(row));
beamwidth=2*atan(y3Db/x3Db)*180/pi
% for i=1:15
%   U =B0 *cos(theta).^i;
% x=U.*cos(theta);
% y=U.*sin(theta);
% x3Db1=max(x)/2
% [value,row]= min(abs(x-x3Db));
% y3Db1=abs(y(row))
% beamwidth(i)=2*atan(y3Db1/x3Db1)*180/pi;
% end
    


for i=1:num
% Given data vectors X and Y.  
% Want to rotate the data  through angle "ang" about rotation center Xc, Yc
X = x;
Y = y;
ang = 2*i*pi/num-pi*P/num*100;
% Specify the coordinates of the center of rotation
Xc = R/2 ;  % Rotate about the 1/4 chord point
Yc = 0 ;
% The data is roated in a three-step process
% Step 1) Shift the data to the rotation center
Xs = X - Xc;  % shifted data
Ys = Y - Yc;
% Step 2) Rotate the data
Xsr =  Xs*cos(ang) + Ys*sin(ang);    % shifted and rotated data
Ysr = -Xs*sin(ang) + Ys*cos(ang);    %
% Step 3) Un-shift the data (back to the original coordinate system)
Xr = Xsr + Xc;  % Rotated data
Yr = Ysr + Yc;
Xr1(i,:) = Xsr + Xc;  % Rotated data
Yr1(i,:)= Ysr + Yc;
if i<7
hold on
figure(3)
plot(xc,yc,Xr,Yr,x,y)
end
end

%create a grid
% given values
pos = [0 0 ;    % startpoint
       0 R/2 ] ;  % endpoint
nturns = 100 ;    % number of turns (integer value)
% engine
dp = diff(pos,1,1) ;
R = hypot(dp(1), dp(2)) ;
phi0 = atan2(dp(2), dp(1)) ;
phi = linspace(0, nturns*2*pi, resolution) ; % 10000 = resolution
r = linspace(0, R, numel(phi)) ;
xg = pos(1,1) + r .* cos(phi + phi0)+R ;
yg = pos(1,2) + r  .* sin(phi + phi0) ;
figure(4)
plot(xg,yg) ; % nturns crossings, including end point

%calculate the intersection 
L=[];

for l=1:resolution
        for i=1:6
    j(i)=inpolygon(xg(l),yg(l),Xr1(i,:),Yr1(i,:));
        end
        k(l)=sum(j);
end
for l=1:resolution
       if k(l)==6
          L=[L ; xg(l) yg(l)];
       end
end

figure(5)
plot(xc,yc,L(:,1),L(:,2))

for i=1:num
    %rotating the intersection dats
% Given data vectors X and Y.  
% Want to rotate the data  through angle "ang" about rotation center Xc, Yc
X = L(:,1);
Y = L(:,2);
ang = 2*i*pi/num;
% Specify the coordinates of the center of rotation
Xc = R ;  % Rotate about the 1/4 chord point
Yc = 0 ;
% The data is roated in a three-step process
% Step 1) Shift the data to the rotation center
Xs = X - Xc;  % shifted data
Ys = Y - Yc;
% Step 2) Rotate the data
Xsr =  Xs*cos(ang) + Ys*sin(ang);    % shifted and rotated data
Ysr = -Xs*sin(ang) + Ys*cos(ang);    %
% Step 3) Un-shift the data (back to the original coordinate system)
Xr = Xsr + Xc;  % Rotated data
Yr = Ysr + Yc;
Xr2(i,:) = Xsr + Xc;  % Rotated data
Yr2(i,:)= Ysr + Yc;
hold on
figure(8)
plot(xc,yc,Xr,Yr)
end

%calculating beamwidth of intersection 
x3Db1=max(Xr)/2;
[value,row]= min(abs(Xr-x3Db1));
y3Db1=sqrt(Yr(row)^2+(Xr(row)-R)^2);
beamwidth_intersection=2*atan(y3Db1/x3Db1)*180/pi

%calculating the ring
in=inpolygon(Xr2(1,:),Yr2(1,:),Xr2(2,:),Yr2(2,:));
for i=length(in):-1:1
if in(i)==1
    row1=i;
    break
end
end
ring=abs(Yr2(1,row1)-R)




