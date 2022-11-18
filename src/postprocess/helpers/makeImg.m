function [I, J] = makeImg(sizex, sizey, helem)
xx = 1:sizex/helem;
yy = 1:sizey/helem;
[JJ, II] = meshgrid(xx, yy);
I = reshape(II, [], 1);
J = reshape(JJ, [], 1);
end