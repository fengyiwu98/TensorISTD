function [gout, pix_seq, gs_mapping, gs_weight, gs_entropy] = pixSeq_LWE(cent_pix, pix_seq, gs_qty)
    pix_qty = size(pix_seq, 2);
    
    % calculate quantity weight of pixel grayscales in local neighborhood
    gs_weight = gs_qty/pix_qty;
    gs_entropy = zeros(1, 256);
    
    % filter valid grayscales with at least 1 pixel in local patch
    gs_mapping = zeros(1, pix_qty);
    gs_typ = 0;
    for i = 1:256
        pix_val = i - 1;
        if gs_qty(i)>0
            % enum valid grayscale
            gs_typ = gs_typ + 1;
            gs_mapping(gs_typ) = pix_val;
            
            % update enropy of current grayscale
            wgt = gs_weight(i);
            gs_entropy(i) = wgt*log2(wgt);
        end
    end
    gs_mapping = gs_mapping(1, 1:gs_typ);
    
    % calculate Weighted pixel value of all sampled pixels with formula given in 3
    H_pix = zeros(1, pix_qty);
    for i = 1:pix_qty
        pval = pix_seq(3, i);
        pval_idx = pval + 1;
        %gs_idx = gs_mapping(pval);
        for j = 1:gs_typ
            map_gs_idx = gs_mapping(j) + 1;
            H_pix(i) = H_pix(i) - (map_gs_idx - pval_idx)^2 * gs_entropy(map_gs_idx);
        end
    end
    pix_seq = [pix_seq; H_pix];
    
    % calculate enhanced central pixel
    gout = 0;
    for i = 1:pix_qty
        gout = gout + H_pix(i)^2;
    end
end