# Colocalization Image Analysis


## 1. Image Preparation

Images first have to be processed in ImageJ. In particular you have to:
- Choose one slice only from the Z-stack (do not do Z-projection because this creates false overlaps!)
- Create RGB images where each channel corresponds to a different protein (Merge Channels -> Make Composite)
- Create binary masks images of your cells (8-bit)


## 2. Image Analysis

Images can now be analysed in Python using the script "colocalization_analysis".


