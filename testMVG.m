
% function A=illumination(R,W,L,P,x1,y1);
R=10;W=3;L=10;P=0;
x1=7;
y1=0;
l=0;
% if x1>=0 && x1<=R && y1<=R/2 && y1>=-R/2 && L<=R
  %cercle
theta = linspace(-pi/2, pi/2);
rc =R *cos(theta);
xc=rc.*cos(theta);
yc=rc.*sin(theta);  

beam=atand(W/L) %beem
num=ceil(360/(beam-beam*(2*P/100)))+1 %number of antennas

% Given data vectors X and Y.  
% Want to rotate the data  through angle "ang" about rotation center Xc, Yc
X = x;
Y = y;
ang = pi/2;
% Specify the coordinates of the center of rotation
Xc = 0 ;  % Rotate about the 1/4 chord point
Yc = 0 ;
% The data is roated in a three-step process
% Step 1) Shift the data to the rotation center
Xs = X - Xc;  % shifted data
Ys = Y - Yc;
% Step 2) Rotate the data
Xsr =  Xs*cos(ang) + Ys*sin(ang);    % shifted and rotated data
Ysr = -Xs*sin(ang) + Ys*cos(ang);    %
% Step 3) Un-shift the data (back to the original coordinate system)
x = Xsr + Xc;  % Rotated data
y = Ysr + Yc;
figure(1)
plot(xc,yc,x,y)
title('visual representation of the scanned circle and the beam of the radiation pattern')


%rotation of the radiation pattern
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
hold on
figure(2)
plot(xc,yc,Xr,Yr,x,y,x1,y1,'c*')
end

for i=1:num
    j=inpolygon(x1,y1,Xr1(i,:),Yr1(i,:));
if j==1     
    l=l+1;
   end

end
disp('number of antennas illuminating the point is')
    disp(l)
% else
% 
%  error('Error. \L Input must be  changed if not the radiation pattern is not covring the centre of the circle or thhe point is not inside the circle.')
%     return 
% end