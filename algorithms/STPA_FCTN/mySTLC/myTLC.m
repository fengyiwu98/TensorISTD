function TLC = myTLC(file_path, ext_name, img_idx)

path_cur = fullfile([file_path, num2str(img_idx,'%04d'), ext_name]);
img_cur = imread(path_cur); 
if ndims( img_cur ) == 3
    img_cur = double( rgb2gray(img_cur) );
end
[imgR, imgC] = size(img_cur);

frame_data = zeros(imgR, imgC, 4);
motionVect = zeros(2, 4);%每个相对帧的运动矢量为（2，1）（代价最小的块相对于宏块左上角的坐标的行列坐标）
%第一行为列，第二行为行坐标

for k = -2: 2
    
    temp = imread(fullfile([file_path, num2str(img_idx+k,'%04d'), ext_name]));
    if ndims(temp) == 3  
        temp = rgb2gray(temp);
    end
    
    if k < 0
        frame_data(:,:,k+3) = double(temp);
        motionVect(:,k+3) = get_motionVect_rand(img_cur, temp, 64, 10, -8*k);
    end
    if k > 0
        frame_data(:,:,k+2) = double(temp);
        motionVect(:,k+2) = get_motionVect_rand(img_cur, temp, 64, 10, 8*k);
    end
    
end

S = zeros(imgR, imgC);
 %%%%% 根据运动向量获取坐标？   
for y = 1:imgR
	for x = 1:imgC
        y1 = min( max(y+motionVect(1,1), 1), imgR);
        y2 = min( max(y+motionVect(1,2), 1), imgR);
        y3 = min( max(y+motionVect(1,3), 1), imgR);
        y4 = min( max(y+motionVect(1,4), 1), imgR);
        x1 = min( max(x+motionVect(2,1), 1), imgC);
        x2 = min( max(x+motionVect(2,2), 1), imgC);
        x3 = min( max(x+motionVect(2,3), 1), imgC);
        x4 = min( max(x+motionVect(2,4), 1), imgC);
        mt = 0.25 * ( frame_data(y1,x1,1) + frame_data(y4,x4,4) ) / 2  + ...
             0.75 * ( frame_data(y2,x2,2) + frame_data(y3,x3,3) ) / 2;
       % S(y,x) = abs( img_cur(y,x) - mt );
       S(y,x) = max( img_cur(y,x) - mt,0 );
	end
end
    
TLC = S  / max(S(:));

%imshow(TLC);
    
end