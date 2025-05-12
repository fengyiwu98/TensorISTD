# TensorISTD
TensorISTD is a Matlab-based, open-source, and user-friendly toolbox designed for optimization-based infrared small target detection (ISTD).

This toolbox introduces a streamlined pipeline to test detection methods. It establishes a benchmark for comprehensively evaluating the performance of existing optimization-based approaches—whether single- or multi-frame, matrix- or tensor-based. It supports mainstream metrics such as SCRG, BSF, CG, and the 3-D ROC series.

TensorISTD empowers researchers with quick access to optimization-oriented infrared small target detection tools while encouraging the development of novel methods. We warmly invite contributors to enrich the benchmark by sharing their own techniques.

Note: This repository will be updated regularly, so stay tuned for improvements and new features!

## Requirement
Matlab 2021b or higher.


## 🛠️How to Use

<details>
<summary>📂Document Structure </summary>
  
```
TensorISTD/
├── algorithms/               # Core detection algorithms
│   ├── 4D_ISTD/              # 4D infrared small target detection
│   ├── ASTTV-NTLA/           # 
│   ├── IPI/                  # 
│   ├── LogTFNN/              # 
│   ├── MPCM/                 # 
│   ├── NPTST/                # 
│   ├── PSTNN/                # 
│   ├── RIPT/                 # 
│   ├── STPA_FCTN/            # 
│   ├── STT/                  # 
│   └── Tophat/               # 
├── all_results/              # Raw detection outputs
│   └── sequence1/            # Sample sequence
│       └── STPA_FCTN/        # Algorithm-specific results
├── fig_results/              # Visualization outputs
│   └── sequence1/            # Sequence-specific figures
├── mat_results/              # Performance metrics
│   ├── curve_results/        # .mat documents
│   └── index_results/        # Metrics records(3DROC,SCR, etc.)
├── time_results/             # Runtime logs
│    └── STPA-FCTN.txt        # Time consumption record
├── utils/                    # Support utilities
│   ├── analyse_pts.m         # Point analysis script
│   ├── curves_draw.m         # Visualization generator
│   ├── get_algo_result.m     # Algorithm output collector
│   ├── get_curves.m          # Metrics calculator
│   ├── measure_calculator.m  # SCR/CG/BSF/BSR calculator
│   └── pt_nms.m              # Non-maximum suppression
├── LICENSE                   # open-source license
├── README.md                  
└── evaluation.m              
```
  
</details>

### 🚀Get Results

- Select the desired algorithm name and datasets：
```matlab
%% evaluation.m
%% user configs
eval_algo_names = ...
    {
     'IPI' %'TT','TR','MPCM','IPI','LogTFNN','NPSTT','PSTNN','RIPT','STT','Tophat'
    };

eval_data_names = ...
    {
     'sequence1' %'sequence1','sequence4','sequence7','sequence9','sequence11','sequence12'
    };
```
- Execute the code in step 1：
```matlab
%% evaluation.m
%% step 1: Use all the evaluation algorithms to get the result plots for all the evaluation datasets present in . /result in the mat file
get_algo_result(eval_algo_names, eval_data_names, ...
     img_types, algo_base_path, data_base_path, res_base_path, time_path );
```
- Fill in the the following configuration：
```matlab
%% evaluation.m
img_types = {'*.jpg', '*.bmp', '*.png'}; % Image Type
algo_base_path = '.\algorithms\'; % Algorithm Path
data_base_path = '.\dataset/data\'; % Datasets Path
res_base_path =  '.\all_result\'; % Fig Results Path
time_path = '.\time_results\'; % Time Path
```
💡`Note`
If the metrics are calculated directly from the existing test image, then comment out this section and go directly to step 2.

### 📈Evaluation
 - The following metrics can be obtained：
   
 ✅SCRG, ✅CG, ✅BSF, ✅BSR

 ✅ Multi-perspective AUC Analysis($AUC_{FPR,TPR}$</th>, $AUC_{\tau,TPR}$</th>, $AUC_{\tau,FPR}$</th>, $AUC_{ODP}$</th>, $AUC_{SNPR}$</th>, $AUC_{TD}$</th>, $AUC_{BS}$</th>, $AUC_{TDBS}$</th>)

   - Calculate the corresponding .mat file based on step 1 or the existing result plots with the target coordinates of the dataset and store it in the curve index folder.
   ```matlab
   %% evaluation.m
   %% step 2: Combining the result plots from step 1 and . GT under /dataset/ann/ to get the roc sequence results
   get_curves(eval_algo_names, eval_data_names, thres_num, radius, res_base_path, ...
        mat_base_path, txt_base_path, mask_base_path, preimg_type);
   ```
   - Related Configurations：
   ```matlab
   %% evaluation.m
   mat_base_path = '.\mat_results\'; % Storage Path for .mat Files
   txt_base_path =  './dataset/anno/'; % Target coordinates Path
   preimg_type = '*.png'; % Result Image Format:.jpg&.png&.bmp...
   ```
   - Steps 3 and 4 are then performed to obtain the metrics and the 3DROC schematics：
   ```matlab
   %% evaluation.m
   %% step 3: Calculating metrics and plotting 3DROC
   curves_drawer(1, eval_algo_names, eval_data_names, figure_base_path, mat_base_path, x_axis_ratio, FPR_thres);

   %% step 4: Calculate SCRG gain,CG, BSF, BSR
   measure_calculator(eval_algo_names, eval_data_names, data_base_path, res_base_path, ...
       mat_base_path, txt_base_path, img_types, preimg_type);
   ```
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
 - The 3DROC image can be obtained during the execution of step 3.

   - Related Configurations：
   ```matlab
   %% evaluation.m
   figure_base_path = '.\fig_results\'; % 3DROC Figure Path
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
   - 📂The evaluation 3DROC result have the following structure:
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
### 🎨Draw Visualization Images
 - This following script provides a standardized pipeline for generating publication-quality 3D visualizations from 2D images. 

   - Key Configuration Parameters
   1. Figure Window
   ```matlab
   %% Create figure window
   figure('Units', 'pixels', 'Position', [100 100 703 641]); % Set window properties (units/position)
   ```
   2. Typography Settings
   ```matlab
   %% Axes configuration
   set(gca,...
       'linewidth', 1, ...      % Axis line width 1pt
       'Fontname', 'Times New Roman', ... % Western font typeface
       'FontSize', 27);         % Font size 27pt
   ```
   3. Axis Configuration
   ```matlab
   xticks(0:50:250);            % X-axis tick interval 50
   xtickformat('%d');           % X-axis integer format
   xlim([0 250]);               % X-axis display range [0-250]
   xtickangle(0);               % X-tick labels horizontal alignment
   ...
   ```
   4. View Perspective
   ```matlab
   view(-37.5,30);              % Set 3D view perspective (azimuth -37.5°, elevation 30°)
   ```
### Others


## Dataset
We used sequences from [[1]](http://www.csdata.org/en/p/387/).

<table class="tg"><thead>
  <tr>
    <th class="tg-9wq8" rowspan="2">Method</th>
    <th class="tg-c3ow" colspan="11">Averaged Values</th>
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
  </tr></thead>
<tbody>
  <tr>
    <td class="tg-c3ow">IPI</td>
    <td class="tg-c3ow">inf</td>
    <td class="tg-c3ow">inf</td>
    <td class="tg-c3ow">1.8051</td>
    <td class="tg-c3ow">0.9926</td>
    <td class="tg-c3ow">0.8562</td>
    <td class="tg-c3ow">5.6319e-2</td>
    <td class="tg-c3ow">1.8432</td>
    <td class="tg-c3ow">1.5384e2</td>
    <td class="tg-c3ow">1.8488</td>
    <td class="tg-c3ow">0.9870</td>
    <td class="tg-c3ow">0.8506</td>
  </tr>
  <tr>
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
  </tr>
  <tr>
    <td class="tg-c3ow">PSTNN</td>
    <td class="tg-c3ow">inf</td>
    <td class="tg-c3ow">inf</td>
    <td class="tg-c3ow">1.4857</td>
    <td class="tg-c3ow">0.9414</td>
    <td class="tg-c3ow">0.7511</td>
    <td class="tg-c3ow">5.9107e-3</td>
    <td class="tg-c3ow">1.6865</td>
    <td class="tg-c3ow">1.3064e2</td>
    <td class="tg-c3ow">1.6925</td>
    <td class="tg-c3ow">0.9354</td>
    <td class="tg-c3ow">0.7452</td>
  </tr>
  <tr>
    <td class="tg-c3ow">LogTFNN</td>
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
  </tr>
  <tr>
    <td class="tg-c3ow">NPSTT</td>
    <td class="tg-c3ow">inf</td>
    <td class="tg-c3ow">inf</td>
    <td class="tg-c3ow">1.8619</td>
    <td class="tg-c3ow">0.9749</td>
    <td class="tg-c3ow">0.7894</td>
    <td class="tg-c3ow">5.1547e-3</td>
    <td class="tg-c3ow">1.7591</td>
    <td class="tg-c3ow">1.5322e2</td>
    <td class="tg-c3ow">1.7643</td>
    <td class="tg-c3ow">0.9698</td>
    <td class="tg-c3ow">0.7842</td>
  </tr>
  <tr>
    <td class="tg-c3ow">4-D TT</td>
    <td class="tg-c3ow">92.588</td>
    <td class="tg-c3ow">44.404</td>
    <td class="tg-c3ow">2.1198</td>
    <td class="tg-c3ow">0.9994</td>
    <td class="tg-c3ow">0.9310</td>
    <td class="tg-c3ow">5.1442e-3</td>
    <td class="tg-c3ow">1.9253</td>
    <td class="tg-c3ow">1.8112e2</td>
    <td class="tg-c3ow">1.9304</td>
    <td class="tg-c3ow">0.9942</td>
    <td class="tg-c3ow">0.9258</td>
  </tr>
  <tr>
    <td class="tg-c3ow">4-D TR</td>
    <td class="tg-c3ow">105.63</td>
    <td class="tg-c3ow">50.668</td>
    <td class="tg-c3ow">2.1445</td>
    <td class="tg-c3ow">0.9969</td>
    <td class="tg-c3ow">0.9310</td>
    <td class="tg-c3ow">5.0919e-2</td>
    <td class="tg-c3ow">1.9228</td>
    <td class="tg-c3ow">1.8294e2</td>
    <td class="tg-c3ow">1.9279</td>
    <td class="tg-c3ow">0.9918</td>
    <td class="tg-c3ow">0.9259</td>
  </tr>
  <tr>
    <td class="tg-c3ow">STPA-FCTN</td>
    <td class="tg-c3ow">156.40</td>
    <td class="tg-c3ow">73.478</td>
    <td class="tg-c3ow">2.1593</td>
    <td class="tg-c3ow">0.9994</td>
    <td class="tg-c3ow">0.9587</td>
    <td class="tg-c3ow">5.0851e-3</td>
    <td class="tg-c3ow">1.9530</td>
    <td class="tg-c3ow">1.8866e2</td>
    <td class="tg-c3ow">1.9581</td>
    <td class="tg-c3ow">0.9944</td>
    <td class="tg-c3ow">0.9536</td>
  </tr>
</tbody></table>

## References
