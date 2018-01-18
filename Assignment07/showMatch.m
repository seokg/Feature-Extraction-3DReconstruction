function showMatch(pt1,pt2,linetype)
x1 = pt1(:,1);
y1 = pt1(:,2);
x2 = pt2(:,1);
y2 = pt2(:,2);

plot(x1,y1,'yo')
plot(x2,y2,'yo')
for i = 1: length(pt1)
stepX = (- x1(i) + x2(i))/40;
stepY = (- y1(i) + y2(i))/40;
x = x1(i) : stepX : x2(i);
y = y1(i) : stepY : y2(i);
if length(x) ~= length(y)
   if length(x) > length(y)
       difflen =length(x) - length(y);
       plot(x(1:end-difflen),y,linetype)
   else
       difflen =length(y) - length(x);
       plot(x,y(1:end-difflen),linetype)
   end
   
else
    plot(x,y,linetype)

end


end

end