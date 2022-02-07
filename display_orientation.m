function display_orientation(imin,orientim,spacing)
imin=imin/max(max(imin));
imin(imin<0)=0;
imin=1-imin;
[rows, cols] = size(orientim);
len = 0.8*spacing;  % length of orientation lines
s_orient = orientim(1:spacing:rows,1:spacing:cols);
xoff = len/2*cos(s_orient);
yoff = len/2*sin(s_orient);
% Determine placement of orientation vectors
[x,y] = meshgrid(1:spacing:cols,1:spacing:rows);
x = x-xoff;
y = y-yoff;
u = xoff*2;
v = yoff*2;
figure;%contour(imin(end:-1:1,:));
imshow(imin,[])
hold on 
% u1=u;
% v1=v;
% u1(u>0)=u(u>0)+pi/2;
% v1(v>0)=v(v>0)+90;
% 
% u1(u<0)=u(u<0)-90;
% v1(v<0)=v(v<0)-90;

%u=mod(u,360);
%v=mod(v,360);
%quiver(x,y,u,v,0,'.');
%quiver(x,y,u1,v1);
%quiver(x,y,u+3*pi/2,v);
quiver(x,y,u+3*pi/2,v);
%quiver(x,y,u,v);

hold off