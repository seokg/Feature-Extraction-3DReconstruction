function output = funcSobel(image,type)
% sobel operator of the image
% image: original image
% output: sobel img

switch type
    case 'x'
      sobelX = [1 0 -1;...
          2 0 -2;...
          1 0 -1];
      
      output = funcConv2d(image,sobelX);

    case 'y'
        sobelY = [-1 -2 -1 ;...
          0 0 0 ;...
          1 2 1];
      output = funcConv2d(image,sobelY);

    otherwise
        disp('wrong input type')
        return
        
end



end