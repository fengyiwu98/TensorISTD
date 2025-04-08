clc;
clear;
close all;
ext='bmp';
tspan_rng=2;
swind_rng=3;
stlc_k=7.8;%7-10
for seq=1:12
  data=[num2str(seq) '\'];
  path=['D:\code\tensor_code\shiyan\benchmark合集\benchmark_release(wyl)\dataset\data\sequence',data];
  STLCF(path, ext, tspan_rng, swind_rng, stlc_k,seq);
end