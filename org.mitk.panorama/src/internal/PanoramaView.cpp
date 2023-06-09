/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "PanoramaView.h"

const std::string PanoramaView::VIEW_ID = "org.mitk.panorama";
PanoramaView *PanoramaView::m_self = nullptr;

PanoramaView::PanoramaView() : m_Controls(nullptr), m_PointListModel(nullptr)
{
  std::string resultDir = "D:\\ctResult";
  std::string mipFile = resultDir + "\\mip.dcm";
  std::string otsuFile = resultDir + "\\otsu.dcm";
  std::string morphFile = resultDir + "\\morph.dcm";
  std::string thinFile = resultDir + "\\thin.png";
  std::string splineFile = resultDir + "\\spline.dcm";
  std::string panoramaFile = resultDir + "\\panorama.dcm";
  std::string looseROIFile = resultDir + "\\looseROI.nii.gz";
  std::string boxROIFile = resultDir + "\\boxROI.nii.gz";
  std::string toothMaskFile = resultDir + "\\toothMask.nii.gz";
  toothSurfaceFile = resultDir + "\\toothSurface.ply";
  m_self = this;
  this->panorama = new Panorama(resultDir,
                                mipFile,
                                otsuFile,
                                morphFile,
                                thinFile,
                                splineFile,
                                panoramaFile,
                                looseROIFile,
                                boxROIFile,
                                toothMaskFile);
}

PanoramaView::~PanoramaView()
{
  delete panorama;
}

void PanoramaView::SetFocus() {}

void PanoramaView::CreateQtPartControl(QWidget *parent)
{
  m_Controls = new Ui::PanoramaViewControls;
  m_Controls->setupUi(parent);

  m_Controls->dilateEdit->setText("15");
  m_Controls->erodeEdit->setText("20");
  m_Controls->thicknessEdit->setText("70");
  m_Controls->leftPaddingEdit->setText("0");
  m_Controls->rightPaddingEdit->setText("0");
  m_Controls->upPaddingEdit->setText("0");
  m_Controls->downPaddingEdit->setText("0");
  //添加一个节点集合
  m_PointSet = mitk::PointSet::New();
  m_PointSetNode = mitk::DataNode::New();
  m_PointSetNode->SetData(mitk::PointSet::Pointer(m_PointSet));
  m_PointSetNode->SetProperty("name", mitk::StringProperty::New("ControlPoints"));
  m_PointSetNode->SetProperty("opacity", mitk::FloatProperty::New(1));
  //显示点蓝色
  m_PointSetNode->SetColor(0.0, 0.0, 1.0);
  this->GetDataStorage()->Add(m_PointSetNode);
  m_PointListModel = new QmitkPointListModel(m_PointSetNode);
  //允许手动移动点
  if (m_PointSetNode.IsNotNull())
  {
    m_DataInteractor = m_PointSetNode->GetDataInteractor();
    if (m_DataInteractor.IsNull())
    {
      m_DataInteractor = mitk::PointSetDataInteractor::New();
      m_DataInteractor->LoadStateMachine("PointSet.xml");
      m_DataInteractor->SetEventConfig("PointSetConfig.xml");
      m_DataInteractor->SetDataNode(m_PointSetNode);
      qDebug() << "\n\nenable pick point\n\n";
    }
  }
  //信号槽函数
  connect(m_Controls->mipButton, &QPushButton::clicked, this, &PanoramaView::OnProcessMip);
  connect(m_Controls->otsuButton, &QPushButton::clicked, this, &PanoramaView::OnProcessOtsu);
  connect(m_Controls->morphButton, &QPushButton::clicked, this, &PanoramaView::OnProcessMorph);
  connect(m_Controls->splineButton, &QPushButton::clicked, this, &PanoramaView::OnProcessSpline);
  connect(m_Controls->archButton, &QPushButton::clicked, this, &PanoramaView::OnProcessArch);
  connect(m_Controls->panoramaButton, &QPushButton::clicked, this, &PanoramaView::OnProcessPanorama);
  connect(m_Controls->lowerSliceButton, &QPushButton::clicked, this, &PanoramaView::OnSetLowerSlice);
  connect(m_Controls->higherSliceButton, &QPushButton::clicked, this, &PanoramaView::OnSetHigherSlice);
  connect(m_Controls->looseROIButton, &QPushButton::clicked, this, &PanoramaView::OnSegLooseROI);
  connect(m_Controls->boxROIButton, &QPushButton::clicked, this, &PanoramaView::OnSegBoxROI);
  //重置按钮
  connect(m_Controls->resetMipButton, &QPushButton::clicked, this, &PanoramaView::OnResetMip);
  connect(m_Controls->resetMorphButton, &QPushButton::clicked, this, &PanoramaView::OnResetMorph);
  connect(m_Controls->resetArchButton, &QPushButton::clicked, this, &PanoramaView::OnResetArch);
  connect(m_Controls->resetPanoramaButton, &QPushButton::clicked, this, &PanoramaView::OnResetPanorama);
  connect(m_Controls->resetLooseROIButton, &QPushButton::clicked, this, &PanoramaView::OnResetLooseROI);
  connect(m_Controls->segmentButton, &QPushButton::clicked, this, &PanoramaView::OnToothSegment);
  connect(m_Controls->MCButton, &QPushButton::clicked, this, &PanoramaView::OnMarchingCube);
  //十字线设置
  auto renderWindowPart = this->GetRenderWindowPart();
  if (nullptr != renderWindowPart)
    this->RenderWindowPartActivated(renderWindowPart);
}

void PanoramaView::RenderWindowPartActivated(mitk::IRenderWindowPart *renderWindowPart)
{
  m_Sncs.insert(renderWindowPart->GetQmitkRenderWindow("axial")->GetSliceNavigationController());
  m_Sncs.insert(renderWindowPart->GetQmitkRenderWindow("sagittal")->GetSliceNavigationController());
  m_Sncs.insert(renderWindowPart->GetQmitkRenderWindow("coronal")->GetSliceNavigationController());
}

void PanoramaView::RenderWindowPartDeactivated(mitk::IRenderWindowPart *renderWindowPart)
{
  m_Sncs.erase(renderWindowPart->GetQmitkRenderWindow("axial")->GetSliceNavigationController());
  m_Sncs.erase(renderWindowPart->GetQmitkRenderWindow("sagittal")->GetSliceNavigationController());
  m_Sncs.erase(renderWindowPart->GetQmitkRenderWindow("coronal")->GetSliceNavigationController());
}

//获得牙齿数据，生成mip图像，加入ds
void PanoramaView::OnProcessMip()
{
  auto ds = this->GetDataStorage();
  mitk::DataNode::Pointer teethNode = ds->GetNamedNode("teeth");
  if (teethNode == nullptr)
  {
    qDebug() << "there is no teeth\n";
    return;
  }
  int lowerSlice = m_Controls->lowerSliceLabel->text().toInt();
  int higherSlice = m_Controls->higherSliceLabel->text().toInt();
  if (lowerSlice == higherSlice)
  {
    qDebug() << "there is no slice\n";
    return;
  }
  mitk::Image::Pointer teethImage = dynamic_cast<mitk::Image *>(teethNode->GetData());
  int thick = teethImage->GetLargestPossibleRegion().GetSize()[2];
  if (lowerSlice > thick || higherSlice > thick)
  {
    qDebug() << "slice out of range\n";
    return;
  }
  auto teethITKImage = mitk::ImageToItkImage<DCMPixelType, 3>(teethImage->Clone());
  auto mip = panorama->Mip(teethITKImage, lowerSlice, higherSlice);
  if (mip == nullptr)
  {
    qDebug() << "mip is nullptr\n";
    return;
  }
  auto mipMITKImage = mitk::ImportItkImage(mip);
  //加入到ds
  RemoveNode("mip");
  mitk::DataNode::Pointer node = mitk::DataNode::New();
  node->SetName("mip");
  node->SetData(mipMITKImage);
  ds->Add(node);
  std::vector<std::string> show = {"mip"};
  std::vector<std::string> hide = {"teeth",
                                   "otsu",
                                   "morph",
                                   "ControlPoints",
                                   "arch",
                                   "panorama",
                                   "frame",
                                   "looseROI",
                                   "boxROI",
                                   "toothMask",
                                   "toothSurface"};
  std::string locate = "mip";
  ResetViewBasicFunc(show, hide, locate);
}

//获得mip数据，阈值分割，加入ds
void PanoramaView::OnProcessOtsu()
{
  auto ds = this->GetDataStorage();
  mitk::DataNode::Pointer mipNode = ds->GetNamedNode("mip");
  if (mipNode == nullptr)
  {
    qDebug() << "there is no mip\n";
    return;
  }
  mitk::Image::Pointer mipImage = dynamic_cast<mitk::Image *>(mipNode->GetData());
  auto mipITKImage = mitk::ImageToItkImage<DCMPixelType, 2>(mipImage->Clone());
  auto otsu = panorama->Otsu(mipITKImage);
  if (otsu == nullptr)
  {
    qDebug() << "otsu is nullptr\n";
    return;
  }
  auto otsuMITKImage = mitk::ImportItkImage(otsu);
  //加入到ds
  RemoveNode("otsu");
  mitk::DataNode::Pointer node = mitk::DataNode::New();
  node->SetName("otsu");
  node->SetData(otsuMITKImage->Clone());
  node->SetColor(RED);
  ds->Add(node);
  std::vector<std::string> show = {"mip", "otsu"};
  std::vector<std::string> hide = {
    "teeth", "morph", "ControlPoints", "arch", "panorama", "frame", "looseROI", "boxROI", "toothMask", "toothSurface"};
  std::string locate = "otsu";
  ResetViewBasicFunc(show, hide, locate);
}

void PanoramaView::OnProcessMorph()
{
  auto ds = this->GetDataStorage();
  mitk::DataNode::Pointer otsuNode = ds->GetNamedNode("otsu");
  if (otsuNode == nullptr)
  {
    qDebug() << "there is no otsu\n";
    return;
  }
  mitk::Image::Pointer otsuImage = dynamic_cast<mitk::Image *>(otsuNode->GetData());
  auto otsuITKImage = mitk::ImageToItkImage<DCMPixelType, 2>(otsuImage->Clone());
  int dilateR = m_Controls->dilateEdit->text().toInt();
  int erodeR = m_Controls->erodeEdit->text().toInt();
  auto morph = panorama->Morph(otsuITKImage, dilateR, erodeR);
  if (morph == nullptr)
  {
    qDebug() << "morph is nullptr\n";
    return;
  }
  auto morphMITKImage = mitk::ImportItkImage(morph);
  //加入到ds
  RemoveNode("morph");
  mitk::DataNode::Pointer node = mitk::DataNode::New();
  node->SetName("morph");
  node->SetData(morphMITKImage->Clone());
  node->SetColor(GREEN);
  ds->Add(node);
  std::vector<std::string> show = {"mip", "otsu", "morph"};
  std::vector<std::string> hide = {
    "teeth", "ControlPoints", "arch", "panorama", "frame", "looseROI", "boxROI", "toothMask", "toothSurface"};
  std::string locate = "morph";
  ResetViewBasicFunc(show, hide, locate);
}

//获得样条曲线控制点
void PanoramaView::OnProcessSpline()
{
  auto ds = this->GetDataStorage();
  mitk::DataNode::Pointer morphNode = ds->GetNamedNode("morph");
  if (morphNode == nullptr)
  {
    qDebug() << "there is no morph\n";
    return;
  }
  mitk::Image::Pointer morphImage = dynamic_cast<mitk::Image *>(morphNode->GetData());
  auto morphITKImage = mitk::ImageToItkImage<DCMPixelType, 2>(morphImage->Clone());
  std::vector<Panorama::Point> controlPoints = panorama->SplineFitting(morphITKImage);
  if (controlPoints.size() != 11)
  {
    qDebug() << "controlPoints error\n";
    return;
  }
  //先清除一下
  m_PointSet->Clear();
  //要转换成物理坐标
  for (Panorama::Point controlPoint : controlPoints)
  {
    DCMImage2DType::IndexType index;
    index[0] = controlPoint.x;
    index[1] = controlPoint.y;
    itk::Point<double, 2> target;
    morphITKImage->TransformIndexToPhysicalPoint(index, target);
    mitk::PointSet::PointType point;
    point[0] = target[0];
    point[1] = target[1];
    point[2] = 0;

    m_PointSet->InsertPoint(point);
  }
  std::vector<std::string> show = {"mip", "ControlPoints"};
  std::vector<std::string> hide = {
    "teeth", "otsu", "morph", "arch", "panorama", "frame", "looseROI", "boxROI", "toothMask", "toothSurface"};
  std::string locate = "mip";
  ResetViewBasicFunc(show, hide, locate);
}

//生成牙弓线
void PanoramaView::OnProcessArch()
{
  auto ds = this->GetDataStorage();
  mitk::DataNode::Pointer morphNode = ds->GetNamedNode("morph");
  if (m_PointSet->GetSize() != 11 || morphNode == nullptr)
  {
    qDebug() << "there is no control point or morph\n";
    return;
  }
  mitk::Image::Pointer morphImage = dynamic_cast<mitk::Image *>(morphNode->GetData());
  auto morphITKImage = mitk::ImageToItkImage<DCMPixelType, 2>(morphImage->Clone());
  std::vector<Panorama::Point> controlPoints;
  mitk::PointSet::PointsContainer *pointsContainer = m_PointSet->GetPointSet()->GetPoints();
  for (auto it = pointsContainer->Begin(); it != pointsContainer->End(); it++)
  {
    mitk::PointSet::PointType p = it->Value();
    itk::Point<double, 2> point;
    point[0] = p[0];
    point[1] = p[1];
    DCMImage2DType::IndexType index;
    morphITKImage->TransformPhysicalPointToIndex(point, index);
    Panorama::Point target;
    target.x = index[0];
    target.y = index[1];
    controlPoints.push_back(target);
  }
  //先清除一下
  result.clear();
  auto arch = panorama->Arch(controlPoints, morphITKImage, result);
  if (arch == nullptr)
  {
    qDebug() << "arch is nullptr\n";
    return;
  }
  auto archMITKImage = mitk::ImportItkImage(arch);
  //加入到ds
  RemoveNode("arch");
  mitk::DataNode::Pointer node = mitk::DataNode::New();
  node->SetName("arch");
  node->SetData(archMITKImage->Clone());
  node->SetColor(BLUE);
  node->SetOpacity(255);
  ds->Add(node);

  std::vector<std::string> show = {"mip", "ControlPoints", "arch"};
  std::vector<std::string> hide = {
    "teeth", "otsu", "morph", "panorama", "frame", "looseROI", "boxROI", "toothMask", "toothSurface"};
  std::string locate = "mip";
  ResetViewBasicFunc(show, hide, locate);
}

//生成全景图
void PanoramaView::OnProcessPanorama()
{
  auto ds = this->GetDataStorage();
  auto teethNode = ds->GetNamedNode("teeth");
  auto morphNode = ds->GetNamedNode("morph");
  if (teethNode == nullptr || result.empty() || morphNode == nullptr)
  {
    qDebug() << "there is no result or teeth or morph\n";
    return;
  }
  mitk::Image::Pointer teethImage = dynamic_cast<mitk::Image *>(teethNode->GetData());
  auto teethITKImage = mitk::ImageToItkImage<DCMPixelType, 3>(teethImage->Clone());
  int thickness = m_Controls->thicknessEdit->text().toInt();
  interpoMap.clear();
  DCMImage2DType::Pointer panoramaImage = panorama->GeneratePanorama(teethITKImage, result, thickness, interpoMap);
  if (panoramaImage == nullptr)
  {
    qDebug() << "panorama is nullptr\n";
    return;
  }
  //这里拷贝信息很重要，不然坐标会出错
  mitk::Image::Pointer morphImage = dynamic_cast<mitk::Image *>(morphNode->GetData());
  auto morphITKImage = mitk::ImageToItkImage<DCMPixelType, 2>(morphImage->Clone());
  panoramaImage->CopyInformation(morphITKImage);
  auto panoramaMITKImage = mitk::ImportItkImage(panoramaImage);
  //加入到ds
  RemoveNode("panorama");
  mitk::DataNode::Pointer node = mitk::DataNode::New();
  node->SetName("panorama");
  node->SetData(panoramaMITKImage->Clone());
  ds->Add(node);
  std::vector<std::string> show = {"panorama"};
  std::vector<std::string> hide = {"teeth",
                                   "mip",
                                   "otsu",
                                   "morph",
                                   "ControlPoints",
                                   "arch",
                                   "frame",
                                   "looseROI",
                                   "boxROI",
                                   "toothMask",
                                   "toothSurface"};
  std::string locate = "panorama";
  ResetViewBasicFunc(show, hide, locate);
}

//移除旧节点数据
void PanoramaView::RemoveNode(std::string nodeName)
{
  auto ds = this->GetDataStorage();
  auto node = ds->GetNamedNode(nodeName);
  if (node != nullptr)
  {
    ds->Remove(node);
  }
}

void PanoramaView::OnSetLowerSlice()
{
  int slice =
    GetRenderWindowPart()->GetQmitkRenderWindow("axial")->GetSliceNavigationController()->GetSlice()->GetPos();
  m_Controls->lowerSliceLabel->setText(QString::number(slice));
}

void PanoramaView::OnSetHigherSlice()
{
  int slice =
    GetRenderWindowPart()->GetQmitkRenderWindow("axial")->GetSliceNavigationController()->GetSlice()->GetPos();
  m_Controls->higherSliceLabel->setText(QString::number(slice));
}

//分割出宽松ROI
void PanoramaView::OnSegLooseROI()
{
  auto ds = GetDataStorage();
  auto frame = ds->GetNamedNode("frame");
  auto panoramaNode = ds->GetNamedNode("panorama");
  auto teethNode = ds->GetNamedNode("teeth");
  if (frame == nullptr || panoramaNode == nullptr || teethNode == nullptr)
  {
    qDebug() << "there is no frame or panorama or teeth\n";
    return;
  }
  mitk::Image::Pointer panoramaImage = dynamic_cast<mitk::Image *>(panoramaNode->GetData());
  auto panoramaITKImage = mitk::ImageToItkImage<DCMPixelType, 2>(panoramaImage->Clone());
  mitk::Image::Pointer teethImage = dynamic_cast<mitk::Image *>(teethNode->GetData());
  auto teethITKImage = mitk::ImageToItkImage<DCMPixelType, 3>(teethImage->Clone());
  //框出的全景图牙齿范围
  mitk::PlanarFigure *frameRect = dynamic_cast<mitk::PlanarFigure *>(frame->GetData());
  int controlNum = frameRect->GetNumberOfControlPoints();
  std::vector<Panorama::Point> boxPoints;
  for (int idx = 0; idx < frameRect->GetNumberOfControlPoints(); idx++)
  {
    mitk::Point3D point = frameRect->GetWorldControlPoint(idx);
    itk::Point<double, 3> tmp;
    tmp[0] = point[0], tmp[1] = point[1], tmp[2] = point[2];
    DCMImage3DType::IndexType index;
    teethITKImage->TransformPhysicalPointToIndex(tmp, index);
    Panorama::Point result;
    result.x = index[0], result.y = index[1];
    boxPoints.push_back(result);
  }
  DCMImage3DType::Pointer looseROI = panorama->LooseROI(teethITKImage, boxPoints, interpoMap);
  if (looseROI == nullptr)
  {
    qDebug() << "looseROI is nullptr\n";
    return;
  }
  auto looseROIMITK = mitk::ImportItkImage(looseROI);
  //加入到ds
  RemoveNode("looseROI");
  mitk::DataNode::Pointer node = mitk::DataNode::New();
  node->SetName("looseROI");
  node->SetData(looseROIMITK->Clone());
  ds->Add(node);
  std::vector<std::string> show = {"teeth", "looseROI"};
  std::vector<std::string> hide = {
    "mip", "otsu", "morph", "ControlPoints", "arch", "panorama", "frame", "boxROI", "toothMask", "toothSurface"};
  std::string locate = "teeth";
  ResetViewBasicFunc(show, hide, locate);
  m_Controls->leftPaddingEdit->setText("0");
  m_Controls->rightPaddingEdit->setText("0");
  m_Controls->upPaddingEdit->setText("0");
  m_Controls->downPaddingEdit->setText("0");
}

//分割出全景图框选
void PanoramaView::OnSegBoxROI()
{
  auto ds = GetDataStorage();
  auto frame = ds->GetNamedNode("frame");
  auto panoramaNode = ds->GetNamedNode("panorama");
  auto teethNode = ds->GetNamedNode("teeth");
  if (frame == nullptr || panoramaNode == nullptr || teethNode == nullptr || interpoMap.empty())
  {
    qDebug() << "there is no frame or panorama or teeth or interpoMap\n";
    return;
  }
  mitk::Image::Pointer panoramaImage = dynamic_cast<mitk::Image *>(panoramaNode->GetData());
  auto panoramaITKImage = mitk::ImageToItkImage<DCMPixelType, 2>(panoramaImage->Clone());
  mitk::Image::Pointer teethImage = dynamic_cast<mitk::Image *>(teethNode->GetData());
  auto teethITKImage = mitk::ImageToItkImage<DCMPixelType, 3>(teethImage->Clone());
  mitk::PlanarFigure *rectangle = dynamic_cast<mitk::PlanarFigure *>(frame->GetData());
  int controlNum = rectangle->GetNumberOfControlPoints();
  //用户框出的牙齿四个角坐标
  std::vector<Panorama::Point> boxPoints;
  for (int idx = 0; idx < controlNum; idx++)
  {
    mitk::Point3D point = rectangle->GetWorldControlPoint(idx);

    itk::Point<double, 2> tmp;
    tmp[0] = point[0], tmp[1] = point[1];
    DCMImage2DType::IndexType index;
    panoramaITKImage->TransformPhysicalPointToIndex(tmp, index);
    Panorama::Point target;
    target.x = index[0];
    target.y = index[1];
    boxPoints.push_back(target);
  }
  int leftPadding = m_Controls->leftPaddingEdit->text().toInt();
  int rightPadding = m_Controls->rightPaddingEdit->text().toInt();
  int upPadding = m_Controls->upPaddingEdit->text().toInt();
  int downPadding = m_Controls->downPaddingEdit->text().toInt();
  DCMImage3DType::Pointer boxROI =
    panorama->BoxROI(teethITKImage, boxPoints, interpoMap, leftPadding, rightPadding, upPadding, downPadding);
  if (boxROI == nullptr)
  {
    qDebug() << "boxROI is nullptr\n";
    return;
  }
  auto boxROIMITK = mitk::ImportItkImage(boxROI);
  //加入到ds
  RemoveNode("boxROI");
  mitk::DataNode::Pointer node = mitk::DataNode::New();
  node->SetName("boxROI");
  node->SetData(boxROIMITK->Clone());
  ds->Add(node);
  std::vector<std::string> show = {"teeth", "boxROI"};
  std::vector<std::string> hide = {
    "mip", "otsu", "morph", "ControlPoints", "arch", "panorama", "frame", "looseROI", "toothMask", "toothSurface"};
  std::string locate = "teeth";
  ResetViewBasicFunc(show, hide, locate);
}

void PanoramaView::OnToothSegment()
{
  RemoveNode("toothMask");
  qDebug() << "tooth segment is not available on Windows system\n";
  return;
}

//三维表面重建
void PanoramaView::OnMarchingCube()
{
  auto ds = GetDataStorage();
  auto toothMaskNode = ds->GetNamedNode("toothMask");
  if (toothMaskNode == nullptr)
  {
    qDebug() << "there is no tooth mask\n";
    return;
  }
  mitk::Image::Pointer toothMaskImage = dynamic_cast<mitk::Image *>(toothMaskNode->GetData());
  //图割算法提取表面
  mitk::GraphcutSegmentationToSurfaceFilter::Pointer surfaceFilter = mitk::GraphcutSegmentationToSurfaceFilter::New();
  surfaceFilter->SetUseMedian(false);
  surfaceFilter->SetMedianKernelSize(1.0f, 1.0f, 1.0f);

  surfaceFilter->SetUseGaussianSmoothing(false);
  surfaceFilter->SetGaussianStandardDeviation(2.2f);
  surfaceFilter->SetGaussianRadius(0.2f);

  surfaceFilter->SetThreshold(1.0f);

  surfaceFilter->SetSmooth(true);
  surfaceFilter->SetSmoothIteration(50);
  surfaceFilter->SetSmoothRelaxation(0.1f);
  surfaceFilter->SetInput(toothMaskImage);
  try
  {
    surfaceFilter->UpdateLargestPossibleRegion();
    mitk::Surface::Pointer toothSurface = surfaceFilter->GetOutput();
    mitk::IOUtil::Save(toothSurface, toothSurfaceFile);
    //添加到节点中
    RemoveNode("toothSurface");
    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetName("toothSurface");
    node->SetData(toothSurface->Clone());
    ds->Add(node);
  }
  catch (mitk::Exception &e)
  {
    std::cerr << e << std::endl;
    return;
  }
  std::vector<std::string> show = {"toothSurface"};
  std::vector<std::string> hide = {
    "teeth", "mip", "otsu", "morph", "ControlPoints", "arch", "panorama", "frame", "looseROI", "boxROI", "toothMask"};
  std::string locate = "toothSurface";
  ResetViewBasicFunc(show, hide, locate);
}

void PanoramaView::OnResetMip()
{
  RemoveNode("frame");
  std::vector<std::string> show = {"teeth", "mip"};
  std::vector<std::string> hide = {
    "otsu", "morph", "ControlPoints", "arch", "panorama", "frame", "looseROI", "boxROI", "toothMask", "toothSurface"};
  std::string locate = "teeth";
  ResetViewBasicFunc(show, hide, locate);
}

void PanoramaView::OnResetMorph()
{
  RemoveNode("frame");
  std::vector<std::string> show = {"mip", "otsu", "morph"};
  std::vector<std::string> hide = {
    "teeth", "ControlPoints", "arch", "panorama", "frame", "looseROI", "boxROI", "toothMask", "toothSurface"};
  std::string locate = "morph";
  ResetViewBasicFunc(show, hide, locate);
}

void PanoramaView::OnResetArch()
{
  RemoveNode("frame");
  std::vector<std::string> show = {"mip", "arch", "ControlPoints"};
  std::vector<std::string> hide = {
    "teeth", "otsu", "morph", "panorama", "frame", "looseROI", "boxROI", "toothMask", "toothSurface"};
  std::string locate = "mip";
  ResetViewBasicFunc(show, hide, locate);
}

void PanoramaView::OnResetPanorama()
{
  RemoveNode("frame");
  std::vector<std::string> show = {"panorama"};
  std::vector<std::string> hide = {"teeth",
                                   "mip",
                                   "otsu",
                                   "morph",
                                   "ControlPoints",
                                   "arch",
                                   "frame",
                                   "looseROI",
                                   "boxROI",
                                   "toothMask",
                                   "toothSurface"};
  std::string locate = "panorama";
  ResetViewBasicFunc(show, hide, locate);
}

void PanoramaView::OnResetLooseROI()
{
  RemoveNode("frame");
  std::vector<std::string> show = {"panorama"};
  std::vector<std::string> hide = {"teeth",
                                   "mip",
                                   "otsu",
                                   "morph",
                                   "ControlPoints",
                                   "arch",
                                   "frame",
                                   "looseROI",
                                   "boxROI",
                                   "toothMask",
                                   "toothSurface"};
  std::string locate = "panorama";
  ResetViewBasicFunc(show, hide, locate);
}

//重新设置Image Navigator的显示
void PanoramaView::ResetViewBasicFunc(std::vector<std::string> show, std::vector<std::string> hide, std::string locate)
{
  auto ds = GetDataStorage();
  for (auto name : show)
  {
    auto node = ds->GetNamedNode(name);
    if (node != nullptr)
    {
      node->SetVisibility(true);
    }
  }
  for (auto name : hide)
  {
    auto node = ds->GetNamedNode(name);
    if (node != nullptr)
    {
      node->SetVisibility(false);
    }
  }
  auto locateNode = ds->GetNamedNode(locate);
  if (locateNode != nullptr)
  {
    if (locate == "toothSurface")
    {
      mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(ds);
    }
    else
    {
      auto image = dynamic_cast<mitk::Image *>(locateNode->GetData());
      mitk::RenderingManager::GetInstance()->InitializeViews(image->GetTimeGeometry());
    }
  }
  this->RequestRenderWindowUpdate();
}