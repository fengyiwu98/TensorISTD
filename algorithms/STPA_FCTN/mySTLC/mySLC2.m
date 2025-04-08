function SLC = mySLC2(file_path, ext_name, img_idx)
%单核
path = fullfile([file_path, num2str(img_idx,'%04d'), ext_name]);
img = imread(path); 
if ndims( img ) == 3
    img = rgb2gray(img);
end
img = double(img);

x = [-1, 0, 1];

sigma_o = 1;
y_o = (1./(sqrt(2.*pi).*sigma_o)) .* exp(-x.^2./(2.*sigma_o.^2));
w_o = y_o ./ sum(y_o);

sigma_i = 0.3;
y_i = (1./(sqrt(2.*pi).*sigma_i)) .* exp(-x.^2./(2.*sigma_i.^2));
w_i = y_i ./ sum(y_i);

m1 = zeros(5, 5);
%m1(1,1)=-w_o(1); m1(2,2)=-w_o(2); m1(3,3)=-w_o(3);
m1(3,1)=w_i(1); m1(3,2)=w_i(2); m1(3,3)=w_i(3);

m2 = zeros(5, 5);
%m2(1,5)=-w_o(1); m2(2,5)=-w_o(2); m2(3,5)=-w_o(3);
m2(1,1)=w_i(1); m2(2, 2)=w_i(2); m2( 3,3)=w_i(3);

m3 = zeros(5, 5);
%m3(1,9)=-w_o(1); m3(2,8)=-w_o(2); m3(3,7)=-w_o(3);
m3(1,3)=w_i(1); m3(2,3)=w_i(2); m3( 3,3)=w_i(3);

m4 = zeros(5, 5);
%m4(5,9)=-w_o(1); m4(5,8)=-w_o(2); m4(5,7)=-w_o(3);
m4(1,5)=w_i(1); m4(2,4)=w_i(2); m4(3,3)=w_i(3);

m5 = zeros(5, 5);
%m5(9,9)=-w_o(1); m5(8,8)=-w_o(2); m5(7,7)=-w_o(3);
m5(3,5)=w_i(1); m5(3,4)=w_i(2); m5(3,3)=w_i(3);

m6 = zeros(5, 5);
%m6(9,5)=-w_o(1); m6(8,5)=-w_o(2); m6(7,5)=-w_o(3);
m6(5,5)=w_i(1); m6(4,4)=w_i(2); m6(3,3)=w_i(3);

m7 = zeros(5, 5);
%m7(9,1)=-w_o(1); m7(8,2)=-w_o(2); m7(7,3)=-w_o(3);
m7(5,3)=w_i(1); m7(4,3)=w_i(2); m7(3,3)=w_i(3);

m8 = zeros(5, 5);
%m8(5,1)=-w_o(1); m8(5,2)=-w_o(2); m8(5,3)=-w_o(3);
m8(5,1)=w_i(1); m8(4,2)=w_i(2); m8(3,3)=w_i(3);

S1 = abs(imfilter(img, m1, 'replicate'));
S2 = abs(imfilter(img, m2, 'replicate'));
S3 = abs(imfilter(img, m3, 'replicate'));
S4 = abs(imfilter(img, m4, 'replicate'));
S5 = abs(imfilter(img, m5, 'replicate'));
S6 = abs(imfilter(img, m6, 'replicate'));
S7 = abs(imfilter(img, m7, 'replicate'));
S8 = abs(imfilter(img, m8, 'replicate'));

S = min( min( min(S1,S2), min(S3,S4) ), min( min(S5,S6), min(S7,S8) ) );
SLC = S / max(S(:));

%imshow(SLC);

end