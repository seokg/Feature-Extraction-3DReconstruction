function output = deformImg(H,img)
%% make deformed image
% input
%   coordinate: coordinate of the new points
%   img: original image
% output 
%   output: new deformed image
[row,col,~]=size(img);
[X,Y]=meshgrid(1:col,1:row);

x = reshape(X,[1,row*col]);
y = reshape(Y,[1,row*col]);
coordinate = [x;y;ones(1,row*col)];
newCoordinate = H * coordinate;
newCoordinate = newCoordinate./repmat(newCoordinate(3,:),[3,1]);

clear x y;
x = reshape(newCoordinate(1,:),[row,col]);
y = reshape(newCoordinate(2,:),[row,col]);

red = img(:,:,1);
green = img(:,:,2);
blue = img(:,:,3);
newred   = uint8(griddata(x,y,double(red)  ,X,Y,'cubic'));
newgreen = uint8(griddata(x,y,double(green),X,Y,'cubic'));
newblue  = uint8(griddata(x,y,double(blue) ,X,Y,'cubic'));

newimg = uint8(zeros(row,col,3));
newimg(:,:,1) = newred;
newimg(:,:,2) = newgreen;
newimg(:,:,3) = newblue;

output = newimg;
end