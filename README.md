# CBCT-interactive-segmentation-tool-based-on-MITK
## 介绍

基于MITK的口腔CBCT牙齿交互式分割工具。

本工作调研了深度学习在医学图像分割领域的研究现状，抓住口腔医疗手术规划流程长、医师专业水平依赖性高的痛点，设计实现了口腔 CBCT 交互式分割工具，用于牙齿正畸等下游任务。

具体工作包括：基于 MITK 开源软件开发了交互式插件，对输入的口腔 CBCT 采用自动+手动调整的方式生成口腔全景图，并在全景图上通过手动框选抽取出单颗牙齿 ROI；基于 nnUNet 框架构建了牙齿分割数据集，并在 Ubuntu 环境下训练得到牙齿分割模型。

经过测试，训练过程 Dice 达到 94%，预测结果可以重建出精度较高的三维牙齿模型。



## 实现方法

主要参考了这篇论文的流程：

《A fully automated method for 3D individual tooth identification and segmentation in dental CBCT》

https://ieeexplore.ieee.org/abstract/document/9445658 



本工作流程图如下，主要包括三部分内容：（1）CBCT口腔全景图生成；（2）基于全景图的牙齿ROI提取；（3）基于nnU-Net的牙齿分割。本工作将前两点内容整合到了MITK插件中，第三点内容要使用Ubuntu操作系统，因此暂时没有整合到插件。

![](https://raw.githubusercontent.com/XiShuFan/picRepo/main/img/1680066967967.jpg)

### 口腔全景图生成

基于CBCT生成口腔全景图，大致分为两步：第一步从CBCT中提取出牙弓线，来表示牙齿的排列；第二步，基于牙弓线，考虑牙齿厚度插值生成全景图。下图是本方法处理流程，分成了5部分。

![](https://raw.githubusercontent.com/XiShuFan/picRepo/main/img/1680067300959.jpg)

### 牙齿ROI提取

基于全景图的牙齿ROI提取，大致分为两步：第一步，在全景图上框选单颗牙齿，获得宽松ROI；第二步，基于宽松ROI，投影分割+手动调整得到紧密ROI。下图是本方法处理流程，分为三部分。

![](https://raw.githubusercontent.com/XiShuFan/picRepo/main/img/1680067384593.jpg)

### 单颗牙齿分割

基于nnU-Net的深度学习牙齿分割，分为两步：第一步，处理训练数据；第二步，牙齿分割训练和预测。

![](https://raw.githubusercontent.com/XiShuFan/picRepo/main/img/1680067583523.jpg)



## 功能展示

最后实现的MITK插件界面如下：

<img src="https://raw.githubusercontent.com/XiShuFan/picRepo/main/img/1680067661314.jpg" style="zoom:80%;" />



全口牙分割完成之后的展示效果（itk-snap三维重建后展示）：

![](https://raw.githubusercontent.com/XiShuFan/picRepo/main/img/1680067747828.jpg)

本工作MITK部分的插件代码在文件夹`org.mitk.panorama`下，可以集成到MITK框架中使用。注意需要修改代码中的数据路径，并且需要手动添加ITK，VTK头文件和库文件。
