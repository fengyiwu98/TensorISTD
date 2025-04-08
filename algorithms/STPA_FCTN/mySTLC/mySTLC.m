function Smap = mySTLC(file_path, ext_name, img_idx)
SLC = mySLC(file_path, ext_name, img_idx);
TLC = myTLC2(file_path, ext_name, img_idx);
Smap = SLC .* TLC;
Smap = ( Smap - min(Smap(:)) ) / (  max(Smap(:)) -  min(Smap(:)) );
imshow(Smap);
end

