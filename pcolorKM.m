 %Front end to pcolor to get it to plot everything
%
%   x:  x-axis coordinate (vector,constant spacing);
%   y:  y-axis coordinate (vector, constant spacing);
%   z:  z-axis coardinate (lengh(x) by length(y) matrix);
%
%   x and y should referance the middle of each grid-point.

function [h] = pcolorKM(x,y,z);


if min(size(x))==1,
  y=y(:);
  x=x(:);
  dy=diff(y(1:2));
  dx=diff(x(1:2));
  y=[y; y(end)+dy];
  x=[x; x(end)+dx];
  y=y-dy/2;
  x=x-dx/2;
  [ly lx]=size(z);
  %if ly~=length(y); z=z';  end;
else,
  dx=max([diff(x(1,1:2)) diff(x(1:2,1))]);
  dy=max([diff(y(1:2,1)) diff(y(1,1:2))]);
  x=[x; x(end,:)];
  x=[x x(:,end)]-dx/2;
  y=[y; y(end,:)];
  y=[y y(:,end)];%-dy/2;
end;

%Duplicate last row and column so that pcolor will plot them
z=[z; z(end,:)];
z=[z z(:,end)];

h = pcolor(x,y,z);

shading flat;