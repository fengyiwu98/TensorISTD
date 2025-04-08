function [ out_img ] = TLC(img_path,ext_name,cur_idx)
img_cur=imread([img_path,num2str(cur_idx,'%04d'),ext_name]);
if ndims(img_cur)==3
    img_cur=(rgb2gray(img_cur));
end

[m,n]=size(img_cur);
register_img=zeros(m,n,4);

for idx=-2:2
    if idx==0
        continue;
    end
    img2 = imread([img_path,num2str(cur_idx+idx,'%04d'),ext_name]);
    if ndims(img2)==3
        img2=rgb2gray(img2);
    end
    [matchLoc1,matchLoc2] = siftMatch(img_cur, img2);
    %[H,corrPtIdx] = findHomography(matchLoc2',matchLoc1');
    tform=fitgeotrans(matchLoc2,matchLoc1,'affine');
    img21= imwarp(img2,tform,'outputView',imref2d(size(img2)));
    mask=img21;
    mask(mask>0)=1;
    if idx<0
        register_img(:,:,idx+3)=mask.*max(img_cur-img21,0);
    else
        register_img(:,:,idx+2)= mask.*max(img_cur-img21,0);
    end
end
out_img=mean(register_img,3);
out_img=out_img./(max(out_img(:)));
imshow(out_img);
end


