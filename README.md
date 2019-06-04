# LungCancerDetection
FIT3081 Lung Cancer Detection using Matlab

LUNG CANCER DETECTION USING CT SCAN IMAGE PROCESSING (MATLAB 2018b)

Methodology

1. Image Acquisition
  - DICOM images are acquired from LIDC dataset (Cancer Imaging Archive)
  
2. Pre-processing:
  - Gabor Filter
  - Active Contour using masking to get Lung Volume
  - Extract nodules using Watershed Segmentation
  
3. Feature Extraction
  - Regionprops and GLCM to obtain nodules properties such as Areas, Perimeter, Mean Pixel Intensity.
  
4. Classification and Visualization of Nodules
  - Nodules are classified basd on TNM 8th edition.

  

