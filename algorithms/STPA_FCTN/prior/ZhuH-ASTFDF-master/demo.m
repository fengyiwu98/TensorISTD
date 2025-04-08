clc;
clear;
close all;
ext='bmp';
ck=1;
tspan=3;
for seq=1:1
  data=[num2str(seq) '\'];
  path=['D:\code\tensor_code\shiyan\shiyan\benchmark合集\benchmark_release(wyl)\dataset\data\sequence',data];
  ASTFDF(path, ext, ck, tspan,seq);
end
