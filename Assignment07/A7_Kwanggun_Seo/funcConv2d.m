function output = funcConv2d(img,kernel)
%% function 2D convolution of image
% img: original image
% kernel: filter/kernel
% output: output image of the convolved image
[row, col , ~] = size(img);
[rowK,colK, ~] = size(kernel);
output = zeros(row,col);
tempImg = zeros(row+2,col+2);
tempImg(2:end-1,2:end-1) = img;
for i = 1:row
   for j = 1:col
       
       for a = 1:rowK
           for b = 1:colK
               output(i,j) = output(i,j) + tempImg(i+a-1,j+b-1)*kernel(a,b);
           end
       end
       
   end
end
output = uint8(output);

end
