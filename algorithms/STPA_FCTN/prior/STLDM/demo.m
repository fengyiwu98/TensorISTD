clc;
clear ;
close all;
ext='bmp';
blk_rnk=1;
tstep=2;
segK=15;
for seq=1:3
  data=[num2str(seq) '\'];
  image_path=['D:\code\tensor_code\shiyan\benchmarkºÏ¼¯\benchmark_release(wyl)\dataset\data\sequence',data];
  STLDM(image_path, ext, blk_rnk, tstep, segK, seq);
end
%img=STLDM(path,ext,blk_rnk,tstep,segK);