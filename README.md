# TensorISTD

[![MATLAB](https://img.shields.io/badge/MATLAB-2021b%2B-orange.svg)]()
[![License](https://img.shields.io/badge/license-MIT-green.svg)]()

**TensorISTD** is an **open-source MATLAB toolbox** designed for **optimization-based Infrared Small Target Detection (ISTD)** research and evaluation.

It provides a **unified, streamlined** evaluation pipeline covering:
- Single-frame / Multi-frame tasks
- Matrix / Tensor representations
- Mainstream evaluation metrics (SCRG, BSF, CG, and 3D-ROC series)

**Objectives:**
1. Enable rapid prototyping and verification of optimization-based detection methods.
2. Establish a standardized benchmark for fair comparison and reproducibility.
3. Encourage community contributions to enrich the algorithm library.

---
**Table of Contents**
- [TensorISTD](#tensoristd)
  - [Requirement](#requirement)
  - [Usage Instruction](#usage-instruction)
    - [Get Results](#-get-results)
    - [Evaluation](#-evaluation)
    - [Draw Visualization Images](#-draw-visualization-images)
  - [Evaluation Table](#evaluation-table)
  - [References](#references)
  - [Acknowledgements](#acknowledgements)
 
Note: 
1. This benchmark also includes the code of our profound and powerful baseline, STPA-FCTN, which can be found in [TensorISTD/algorithms/STPA_FCTN](https://github.com/fengyiwu98/TensorISTD/tree/main/algorithms/STPA_FCTN).
2. This repository will be updated regularly, so stay tuned for improvements and new features!
3. If you would like to contribute, please contact us!
---

## Requirement
Matlab 2021b or higher.


## Usage Instruction

<details>
<summary> 📂 Document Structure </summary>
  
```
TensorISTD/
├── algorithms/               # Core detection algorithms
│   ├── 4D_ISTD/              # [1]
│   ├── LogTFNN/              # [2]
│   ├── NPSTT/                # [3]
│   ├── PSTNN/                # [4]
│   └── STPA_FCTN/            # Ours
├── all_results/              # Raw detection outputs
│   └── sequence1/            # Sample sequence
│       └── STPA_FCTN/        # Algorithm-specific results
├── fig_results/              # Visualization outputs
│   └── sequence1/            # Sequence-specific figures
├── mat_results/              # Performance metrics
│   ├── curve_results/        # .mat documents
│   └── index_results/        # Metrics records(3D-ROC,SCR, etc.)
├── time_results/             # Runtime logs
│    └── STPA-FCTN.txt        # Time consumption record
├── utils/                    # Support utilities
│   ├── analyse_pts.m         # Point analysis script
│   ├── curves_draw.m         # Visualization generator
│   ├── get_algo_result.m     # Algorithm output collector
│   ├── get_curves.m          # Metrics calculator
│   ├── measure_calculator.m  # SCR/CG/BSF/BSR calculator
│   └── pt_nms.m              # Non-maximum suppression
├── 3D_Visualization          
├── LICENSE                   # open-source license
├── README.md
├── color1.mat           
└── evaluation.m              
```
  
</details>

### 🚀 Get Results

**Main command**
```matlab
evaluation.m
```
Use all the evaluation algorithms to get the result plots for all the evaluation datasets present in . /result in the mat file.

**1. Inititalization** 

Select the desired algorithm name and datasets：
```matlab
eval_algo_names = ...
    {
     'STPA_FCTN' %'TT','TR','LogTFNN','NPSTT','PSTNN','RIPT','STT'
    };

eval_data_names = ...
    {
     'sequence1' 
    };
```

Fill in the following configuration：
```matlab
img_types = {'*.jpg', '*.bmp', '*.png'}; 
algo_base_path = '.\algorithms\'; 
data_base_path = '.\dataset\data\';
res_base_path =  '.\all_result\'; 
time_path = '.\time_results\';
```

**2. Choose single-frame or multi-frame algorithm：**

Single-frame (PSTNN, etc.)
```matlab
get_algo_result(eval_algo_names, eval_data_names, ...
     img_types, algo_base_path, data_base_path, res_base_path, time_path );
```

Multi-frame (STT, NPSTT, 4-D-TT/TR, STPA-FCTN, etc.)
```matlab
get_algo_result_multiframe(eval_algo_names, eval_data_names, ...
     img_types, algo_base_path, data_base_path, res_base_path, time_path );
```

💡 `Note`

If the metrics are calculated directly from the existing test image, then comment out this section and go directly to **Evaluation**.

### 📈 Evaluation

In this repo, we include the following metrics:

 ✅ SCRG,  ✅ CG, ✅ BSF, 
 
 ✅ Diverse AUC Analysis ([6]):
 $AUC_{(FPR,TPR)}$</th>, $AUC_{(\tau,TPR)}$</th>, $AUC_{(\tau,FPR)}$</th>, $AUC_{ODP}$</th>, $AUC_{SNPR}$</th>, $AUC_{TD}$</th>, $AUC_{BS}$</th>, $AUC_{TDBS}$</th>

 
 **1. Calculate the corresponding .mat file** based on **Get Results** or the existing result plots with the target coordinates of the dataset and store it in the curve index folder.
   ```matlab
   get_curves(eval_algo_names, eval_data_names, thres_num, radius, res_base_path, ...
        mat_base_path, txt_base_path, mask_base_path, preimg_type);
   ```
This step will combine the result plots from step 1 and . GT under /dataset/ann/ to get the roc sequence results

   - Related Configurations：
   ```matlab
   %% evaluation.m
   mat_base_path = '.\mat_results\'; 
   txt_base_path =  '.\dataset\anno\'; 
   preimg_type = '*.png';
   ```
**2. Calculate the Multi-perspective AUC Analysis:** 
   
   ```matlab
   curves_drawer(1, eval_algo_names, eval_data_names, figure_base_path, mat_base_path, x_axis_ratio, FPR_thres);
   ```
   
**3. Calculate the SCRG gain, CG, BSF metrics:** 

   ```matlab
   measure_calculator(eval_algo_names, eval_data_names, data_base_path, res_base_path, ...
       mat_base_path, txt_base_path, img_types, preimg_type);
   ```

**4. Draw the 3D-ROC figures**

   - 📂The evaluation result will be saved in `index_results`:
   ```
   ├──./mat_result/
   │    ├── curve_results
   │    │    ├── sequence1_STPA-FCTN.mat
   │    │    ├── ...
   │    ├── index_results
   │    │    ├── sequence1.txt
   │    │    ├── ...
   ```

 - The 3D-ROC image can be obtained during the execution of step 3.

   - Related Configurations：
   ```matlab
   %% evaluation.m
   figure_base_path = '.\fig_results\';
   x_axis_ratio = 1e-4;
   FPR_thres = 1;
   ```
   ```matlab
   %% utils\curves_drawer.m
   % Line Color
       color_map = [ ...
   %          55/255   126/255  184/255;  % Blue
   %          77/255   175/255  74/255;   % Green
   %          228/255  26/255   28/255;   % Red
   %          255/255  217/255  47/255;   % Yellow
   ];
   % Line Type
   LineType = {':' }; %'-.'
   ```
   - 📂 The evaluation 3D-ROC result has the following structure:
   ```
   ├──./fig_results/
   │    ├── sequence1
   │    │    ├── SOTA_1
   │    │    ├── SOTA_2
   │    │    ├── SOTA_3
   │    │    ├── SOTA_4
   │    ├── sequence2
   │    │    ├── SOTA_1
   │    │    ├── ...
   │    ├── ...
   ```
   - The following figures are 3D-ROC results of STPA-FCTN in seq 1.
   
   <img src="https://github.com/fengyiwu98/TensorISTD/blob/main/fig_results/sequence1/SOTA_1.png" width="180px"><img src="https://github.com/fengyiwu98/TensorISTD/blob/main/fig_results/sequence1/SOTA_2.png" width="180px">
   <img src="https://github.com/fengyiwu98/TensorISTD/blob/main/fig_results/sequence1/SOTA_3.png" width="180px"><img src="https://github.com/fengyiwu98/TensorISTD/blob/main/fig_results/sequence1/SOTA_4.png" width="180px">

   - Comparison of multiple algorithms.
   
   <img src="https://github.com/fengyiwu98/TensorISTD/blob/main/fig_results/fig/algo_1.png" width="180px"><img src="https://github.com/fengyiwu98/TensorISTD/blob/main/fig_results/fig/algo_2.png" width="180px">
   <img src="https://github.com/fengyiwu98/TensorISTD/blob/main/fig_results/fig/algo_3.png" width="180px"><img src="https://github.com/fengyiwu98/TensorISTD/blob/main/fig_results/fig/algo_4.png" width="180px">

### 🎨 Draw Visualization Images
 - This following script provides a standardized pipeline for generating publication-quality 3D visualizations from 2D images. 

   - Key Configuration Parameters
   1. Figure Window
   ```matlab
   figure('Units', 'pixels', 'Position', [100 100 703 641]);
   ```
   2. Axis Configuration
   ```matlab
   xticks(0:50:250);          
   xtickformat('%d');          
   xlim([0 250]);             
   xtickangle(0);         
   ...
   ```
   3. View Perspective
   ```matlab
   view(-37.5,30);              % Set 3D view perspective (azimuth -37.5°, elevation 30°)
   ```
   4. Color Configuration
   
   The colormap is provided as color1.mat.
   - 3D Visualization
     
   <img src="https://github.com/fengyiwu98/TensorISTD/blob/main/fig_results/fig/31.png" width="180px"><img src="https://github.com/fengyiwu98/TensorISTD/blob/main/fig_results/fig/seqvis.png" width="180px">
   <img src="https://github.com/fengyiwu98/TensorISTD/blob/main/fig_results/fig/tar31.png" width="180px"><img src="https://github.com/fengyiwu98/TensorISTD/blob/main/fig_results/fig/tarvis.png" width="180px">


## Evaluation Table


<table class="tg"><thead>
  <tr>
    <th class="tg-9wq8" rowspan="2">Method</th>
    <th class="tg-c3ow" colspan="12">Averaged Values</th>
  </tr>
  <tr>
    <th class="tg-baqh">SCRG</th>
    <th class="tg-baqh">BSF</th>
    <th class="tg-baqh">CG</th>
    <th class="tg-baqh">$AUC_{FPR,TPR}$</th>
    <th class="tg-baqh">$AUC_{\tau,TPR}$</th>
    <th class="tg-baqh">$AUC_{\tau,FPR}$</th>
    <th class="tg-baqh">$AUC_{ODP}$</th>
    <th class="tg-baqh">$AUC_{SNPR}$</th>
    <th class="tg-baqh">$AUC_{TD}$</th>
    <th class="tg-baqh">$AUC_{BS}$</th>
    <th class="tg-baqh">$AUC_{TDBS}$</th>
    <th class="tg-baqh">Time(s)</th>
  </tr></thead>
<tbody>
<!--   <tr>
    <td class="tg-c3ow">TV-PCP</td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
  </tr>
  <tr>
    <td class="tg-c3ow">SMSL</td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
  </tr>
  <tr>
    <td class="tg-c3ow">NRAM</td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
  </tr>
  <tr>
    <td class="tg-c3ow">SRWS</td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
  </tr> -->
  <tr>
    <td class="tg-c3ow">PSTNN</td>
    <td class="tg-c3ow">Inf</td>
    <td class="tg-c3ow">Inf</td>
    <td class="tg-c3ow">1.4788</td>
    <td class="tg-c3ow">0.9429</td>
    <td class="tg-c3ow">0.7011</td>
    <td class="tg-c3ow">6.1038e-3</td>
    <td class="tg-c3ow">1.6249</td>
    <td class="tg-c3ow">1.2251e2</td>
    <td class="tg-c3ow">1.6310</td>
    <td class="tg-c3ow">0.9238</td>
    <td class="tg-c3ow">0.6950</td>
    <td class="tg-c3ow">0.1606</td>
  </tr>
  <tr>
    <td class="tg-c3ow">LogTFNN</td>
    <td class="tg-c3ow">4.5912</td>
    <td class="tg-c3ow">5.1813</td>
    <td class="tg-c3ow">1.0208</td>
    <td class="tg-c3ow">0.9756</td>
    <td class="tg-c3ow">0.5694</td>
    <td class="tg-c3ow">6.5000e-2</td>
    <td class="tg-c3ow">1.5385</td>
    <td class="tg-c3ow">0.9629e2</td>
    <td class="tg-c3ow">1.5450</td>
    <td class="tg-c3ow">0.9629</td>
    <td class="tg-c3ow">0.5629</td>
    <td class="tg-c3ow">1.3208</td>
  </tr>
<!--   <tr>
    <td class="tg-c3ow">Modek1k2</td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
    <td class="tg-c3ow"></td>
  </tr> -->
  <tr>
    <td class="tg-c3ow">NPSTT</td>
    <td class="tg-c3ow">Inf</td>
    <td class="tg-c3ow">Inf</td>
    <td class="tg-c3ow">1.8846</td>
    <td class="tg-c3ow">0.9699</td>
    <td class="tg-c3ow">0.7824</td>
    <td class="tg-c3ow">5.1423e-3</td>
    <td class="tg-c3ow">1.7517</td>
    <td class="tg-c3ow">1.5926e2</td>
    <td class="tg-c3ow">1.7523</td>
    <td class="tg-c3ow">0.9649</td>
    <td class="tg-c3ow">0.7818</td>
    <td class="tg-c3ow">5.1299</td>
  </tr>
  <tr>
    <td class="tg-c3ow">4-D TT</td>
    <td class="tg-c3ow">89.171</td>
    <td class="tg-c3ow">41.837</td>
    <td class="tg-c3ow">2.3124</td>
    <td class="tg-c3ow">0.9921</td>
    <td class="tg-c3ow">0.9277</td>
    <td class="tg-c3ow">5.1489e-3</td>
    <td class="tg-c3ow">1.9147</td>
    <td class="tg-c3ow">1.8021e2</td>
    <td class="tg-c3ow">1.9199</td>
    <td class="tg-c3ow">0.9869</td>
    <td class="tg-c3ow">0.9226</td>
    <td class="tg-c3ow">0.8694</td>
  </tr>
  <tr>
    <td class="tg-c3ow">4-D TR</td>
    <td class="tg-c3ow">100.62</td>
    <td class="tg-c3ow">49.084</td>
    <td class="tg-c3ow">2.2463</td>
    <td class="tg-c3ow">0.9976</td>
    <td class="tg-c3ow">0.9374</td>
    <td class="tg-c3ow">5.0947e-3</td>
    <td class="tg-c3ow">1.9300</td>
    <td class="tg-c3ow">1.8409e2</td>
    <td class="tg-c3ow">1.9351</td>
    <td class="tg-c3ow">0.9926</td>
    <td class="tg-c3ow">0.9324</td>
    <td class="tg-c3ow">2.0373</td>
  </tr>
  <tr>
    <td class="tg-c3ow">STPA-FCTN</td>
    <td class="tg-c3ow">149.26</td>
    <td class="tg-c3ow">68.905</td>
    <td class="tg-c3ow">2.3292</td>
    <td class="tg-c3ow">0.9996</td>
    <td class="tg-c3ow">0.9581</td>
    <td class="tg-c3ow">5.0778e-3</td>
    <td class="tg-c3ow">1.9527</td>
    <td class="tg-c3ow">1.8893e2</td>
    <td class="tg-c3ow">1.9579</td>
    <td class="tg-c3ow">0.9945</td>
    <td class="tg-c3ow">0.9530</td>
    <td class="tg-c3ow">0.9926</td>
  </tr>
</tbody></table>

💡 `Note`

1. In multi-frame evaluation, an "Inf" value does not necessarily indicate that all frames produced infinite results. This can occur if even a single frame yields an "Inf" value. Therefore, for multi-frame scenarios, we recommend using CG or 3-D ROC metrics for more robust and reliable assessment.

2. We used sequences from [[1]](http://www.csdata.org/en/p/387/).


## References
[1] F. Wu, H. Yu, A. Liu, J. Luo, and Z. Peng, “Infrared small target detection using spatiotemporal 4-d tensor train and ring unfolding,” IEEE Trans. Geosci. Remote Sens., vol. 61, pp. 1–22, 2023.

[![](https://img.shields.io/badge/Link-Paper-blue)](https://ieeexplore.ieee.org/document/10156866)


[2] X. Kong, C. Yang, S. Cao, C. Li and Z. Peng, "Infrared Small Target Detection via Nonconvex Tensor Fibered Rank Approximation," in IEEE Transactions on Geoscience and Remote Sensing, vol. 60, pp. 1-21, 2022.

[![](https://img.shields.io/badge/Link-Paper-blue)](https://ieeexplore.ieee.org/abstract/document/9394596)

[3] G. Wang, B. Tao, X. Kong, and Z. Peng, “Infrared small target detection using nonoverlapping patch spatial–temporal tensor factorization with capped nuclear norm regularization,” IEEE Trans. Geosci. Remote Sens., vol. 60,pp. 1–17, 2021.

[![](https://img.shields.io/badge/Link-Paper-blue)](https://ieeexplore.ieee.org/document/9606877)

[4] L. Zhang and Z. Peng, “Infrared small target detection based on partial sum of the tensor nuclear norm,” Remote Sens., vol. 11, no. 4, p. 382, 2019.

[![](https://img.shields.io/badge/Link-Paper-blue)](https://www.mdpi.com/2072-4292/11/4/382)

[5] C. -I. Chang, "An Effective Evaluation Tool for Hyperspectral Target Detection: 3D Receiver Operating Characteristic Curve Analysis," in IEEE Transactions on Geoscience and Remote Sensing, vol. 59, no. 6, pp. 5131-5153, June 2021,

[![](https://img.shields.io/badge/Link-Paper-blue)](https://ieeexplore.ieee.org/abstract/document/9205919)[![](https://img.shields.io/badge/Code-Matlab-orange)](https://umbc.atlassian.net/wiki/spaces/rssipl/pages/27885869/10.+Download)

**Note: If you used the above codes, please cite the relevant paper.**


## Acknowledgements

Despite the organizers of this repo, we would like to thank the former contributors--**Tianfang Zhang, Jian Li, Yuelu Wei, Guanghui Wang, Xuan Kong, Haiyang Yi, Ruochen Qie, Hang Yu, Anran Liu, Simin Liu, and Zhenming Peng--for this Toolbox.**


