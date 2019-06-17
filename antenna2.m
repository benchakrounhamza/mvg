function num= antenna2(R,W,L,P)
clf
clc
%clear all

%R=10;W=3;L=6;P=10;

if L<R/2
error('Error. \L Input must be  changed if not the radiation pattern is not covring the centre of the circle.')
    return
else
%cercle
theta = linspace(-pi/2, pi/2);
rc =R *cos(theta);
xc=rc.*cos(theta);
yc=rc.*sin(theta);  

beem=atand(W/L) %beem
num=ceil(360/(beem-beem*(2*P/100)))+1 %number of antennas

%ellipse
a=L; % horizontal radius
b=L*tan(pi/num); % vertical radius
x0=0; % x0,y0 ellipse centre coordinates
y0=0;
x=x0+a*cos(theta);
y=y0+b*sin(theta);
figure(1)
plot(x,y)
h=plot(x,y,xc,yc);
title('visual representation of the scaned circel and the beem of the radiation pattern')


%rotation of the radiation pattern
for i=1:num
% Given data vectors X and Y.  
% Want to rotate the data  through angle "ang" about rotation center Xc, Yc
X = x;
Y = y;
ang = 2*i*pi/num;
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
hold on
figure(2)
plot(xc,yc,Xr,Yr,x,y)

end
hold off
 %overlap space
S=pi*a*b/2; % surface of half am elips
Sp=S*2*(P/100); %surface over laped sides
p=abs((R-L)/L); %purcentage of overlapping biger than R
St=pi*abs(R-L)*((1+p)*W)/2; %surface overllapng biger than R
So=(S+Sp+St)*num; %surfcae coverred in mm
Sc=pi*R^2; %surface of the circle en mm
p1=So*100/Sc %pourcentage of coveragein %
end

%representation of coverage area with the number of antennas
for i=1:num

S=pi*a*b/2; % surface of half am elips
Sp=S*2*(P/100); %surface over laped sides
p=abs((R-L)/L); %purcentage of overlapping biger than R
St=pi*abs(R-L)*((1+p)*W)/2; %surface overllapng biger than R
So=(S+Sp+St)*i; %surfcae coverred in mm
Sc=pi*R^2; %surface of the circle en mm
p2(1,i)=i;
p2(2,i)=So*100/Sc; %pourcentage of coveragein %

end
figure (3)
plot (p2(1,:),p2(2,:))
title('representation of coverage area with the number of antennas')
xlabel('number of antennas')
ylabel('% of the surface covered')

end
