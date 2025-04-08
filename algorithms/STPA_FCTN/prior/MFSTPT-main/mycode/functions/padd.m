function I1= padd(img,patchSize)
%UNTITLED3 此处显示有关此函数的摘要
%  扩充图像，边界进行复制patch进行填充
I1=padarray(img,[patchSize patchSize],'symmetric','both');
[n1,n2]=size(I1);
I2=I1(:,1:patchSize);
I2=flip(I2,2);%原图像的垂直镜像
I3=I1(:,n2-patchSize+1:n2);
I3=flip(I3,2);
I4=I1(1:patchSize,:);
I4=flip(I4,1);
I5=I1(n1-patchSize+1:n1,:);
I5=flip(I5,1);
I1(:,1:patchSize)=I2;
I1(:,n2-patchSize+1:n2)=I3;
I1(1:patchSize,:)=I4;
I1(n1-patchSize+1:n1,:)=I5;
end

