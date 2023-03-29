/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#ifndef PanoramaView_h
#define PanoramaView_h

#include "CBCTPanorama.h"
#include "mitkGraphcutSegmentationToSurfaceFilter.h"
#include "qstackedwidget.h"
#include "ui_PanoramaViewControls.h"
#include <QString>
#include <QmitkAbstractView.h>
#include <QmitkIOUtil.h>
#include <QmitkPointListModel.h>
#include <QmitkRenderWindow.h>
#include <mitkDataNode.h>
#include <mitkIRenderWindowPartListener.h>
#include <mitkITKImageImport.h>
#include <mitkImageToItk.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkPlanarFigure.h>
#include <mitkPointSet.h>
#include <mitkPointSetDataInteractor.h>

class PanoramaView : public QmitkAbstractView, public mitk::IRenderWindowPartListener
{
  Q_OBJECT

public:
  static const std::string VIEW_ID;
  PanoramaView();
  ~PanoramaView();
  void SetFocus() override;
  void RenderWindowPartActivated(mitk::IRenderWindowPart *renderWindowPart) override;
  void RenderWindowPartDeactivated(mitk::IRenderWindowPart *renderWindowPart) override;

private:
  void CreateQtPartControl(QWidget *parent) override;
  void OnProcessMip();
  void OnProcessOtsu();
  void OnProcessMorph();
  void OnProcessSpline();
  void OnProcessArch();
  void OnProcessPanorama();
  void RemoveNode(std::string nodeName);
  void OnSetLowerSlice();
  void OnSetHigherSlice();
  void OnSegLooseROI();
  void OnSegBoxROI();
  void OnToothSegment();
  void OnMarchingCube();

  //重新设置参数
  void OnResetMip();
  void OnResetMorph();
  void OnResetArch();
  void OnResetPanorama();
  void OnResetLooseROI();

  void ResetViewBasicFunc(std::vector<std::string> show, std::vector<std::string> hide, std::string locate);
  static PanoramaView *m_self;
  Panorama *panorama;
  Ui::PanoramaViewControls *m_Controls;
  mitk::DataNode::Pointer m_PointSetNode;
  mitk::DataInteractor::Pointer m_DataInteractor;
  QmitkPointListModel *m_PointListModel;
  mitk::PointSet::Pointer m_PointSet;
  std::set<mitk::SliceNavigationController *> m_Sncs;
  //保存的样条曲线拟合点
  std::vector<Panorama::Point> result;
  //拟合点---每个拟合点的插值点
  std::vector<std::vector<Panorama::Point>> interpoMap;
  //三维重建牙齿模型文件
  std::string toothSurfaceFile;

  //节点颜色
  const float RED[3] = {255.0, 0, 0};
  const float GREEN[3] = {0, 255.0, 0};
  const float BLUE[3] = {0, 0, 255.0};
  const float GRAY[3] = {127.0, 127.0, 127.0};
};
#endif
