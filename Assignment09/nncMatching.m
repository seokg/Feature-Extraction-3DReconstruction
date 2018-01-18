function index =nncMatching(img1,img2,pts1,pts2,windows)
iter = 0;
index=[];


% check if the points is matched
checkPts1 = zeros(1,length(pts1));
checkPts2 = zeros(1,length(pts2));



value =zeros(length(pts1),length(pts2));
template1 = cell(length(pts1),1);
for i = 1:length(pts1)
    template1{i} = segment(img1,pts1(:,i),windows);
end
template2 = cell(length(pts2),1);
for i = 1:length(pts2)
    template2{i} = segment(img2,pts2(:,i),windows);

end


for i = 1:length(pts1)
    for j = 1:length(pts2)

        tempImg1 = template1{i};
        tempImg2 = template2{j};
%         imshow(uint8([tempImg1,tempImg2]))

        meantempImg1 = sum(sum(tempImg1))/windows.^2;
        meantempImg2 = sum(sum(tempImg2))/windows.^2;
        
        numerator = sum(sum((tempImg1 - meantempImg1) .* (tempImg2-meantempImg2)));
        denominator = sum(sum((tempImg1 - meantempImg1).^2)) .* sum(sum((tempImg2-meantempImg2).^2));
        denominator = sqrt(denominator);
        
        if numerator/denominator > 0
            value(i,j) = numerator/denominator;
        end
        
        % count iter
        if mod(iter,100)==0
            disp(iter);
        end
        iter = iter +1;

    end
    [maxVal,maxIdx]= max(value(i,:));

    if  checkPts2(maxIdx) == 0 
        
        checkPts1(i) = 1;
        checkPts2(maxIdx) = 1;
        index =[index [i;maxIdx]];       
    else
        duplicateIdx = find(index(2,:) == maxIdx);
        if value(index(1,duplicateIdx),index(2,duplicateIdx)) > maxVal
            
            
        else
            % delete
            checkPts2(index(2,duplicateIdx)) = 0;            
            index(:,duplicateIdx) = [];
            
            % add
            checkPts1(i) = 1;
            checkPts2(maxIdx) = 1;
            index =[index [i;maxIdx]]; 

        end
    end
end




end

function output = segment(img,pts,windows)
    [r,c,~]=size(img);
    temp = zeros(r+windows-1,c+windows-1);
    temp(ceil(windows/2):ceil(windows/2)+r-1,ceil(windows/2):ceil(windows/2)+c-1) = img;
    output = zeros(windows);
    ceil(windows/2);
    for i = 1:windows
        for j = 1:windows
            output(i,j) = temp(pts(2)+i, ...
                pts(1)+j);
 
        end
    end
end
