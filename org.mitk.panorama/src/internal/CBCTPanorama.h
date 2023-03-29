#pragma once
/*
生成口腔牙齿全景图
*/
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryThinningImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkMaximumProjectionImageFilter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkNumericTraits.h"
#include "itkOtsuMultipleThresholdsCalculator.h"
#include "itkPNGImageIO.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include <cmath>
#include <ep/include/ITK-4.13/itkExtractImageFilter.h>
#include <itkNiftiImageIO.h>
#include <mitkPoint.h>
#include <string>

// CBCT像素类型使用short
typedef int DCMPixelType;
// PNG像素类型使用unsigned char
typedef unsigned char PNGPixelType;
//标签类型数值使用unsigned short
typedef unsigned short LabelPixelType;
//图像类型
typedef itk::Image<DCMPixelType, 3> DCMImage3DType;
typedef itk::Image<DCMPixelType, 2> DCMImage2DType;
typedef itk::Image<PNGPixelType, 2> PNGImageType;
typedef itk::Image<LabelPixelType, 3> LabelImage3DType;
//大津阈值分割类型
typedef itk::Statistics::ScalarImageToHistogramGenerator<DCMImage2DType> ScalarImageToHistogramGeneratorType;
typedef ScalarImageToHistogramGeneratorType::HistogramType HistogramType;
typedef itk::OtsuMultipleThresholdsCalculator<HistogramType> CalculatorType;
typedef itk::BinaryThresholdImageFilter<DCMImage2DType, DCMImage2DType> BinaryFilterType;
typedef itk::SmoothingRecursiveGaussianImageFilter<DCMImage2DType, DCMImage2DType> Smoothing2DFilterType;
//读写文件类型
typedef itk::GDCMSeriesFileNames NamesGeneratorType;
typedef itk::GDCMImageIO DCMImageIO;
typedef itk::PNGImageIO PNGImageIO;
typedef itk::NiftiImageIO NiftiImageIO;
typedef itk::ImageSeriesReader<DCMImage3DType> seriesReaderType;
typedef itk::ImageSeriesWriter<DCMImage3DType, DCMImage2DType> seriesWriterType;
typedef itk::ImageFileReader<DCMImage2DType> DCMReaderType;
typedef itk::ImageFileWriter<DCMImage2DType> DCMWriterType;
typedef itk::ImageFileWriter<PNGImageType> PNGWriterType;
typedef itk::ImageFileReader<PNGImageType> PNGReaderType;
typedef itk::ImageRegionIterator<PNGImageType> PNGIteratorType;
typedef itk::ImageFileReader<DCMImage3DType> NiftiReaderType;
typedef itk::ImageFileWriter<DCMImage3DType> NiftiWriterType;
typedef itk::ImageFileWriter<LabelImage3DType> LabelWriterType;
//最大强度投影类型
typedef itk::MaximumProjectionImageFilter<DCMImage3DType, DCMImage2DType> MipFilterType;
typedef itk::ExtractImageFilter<DCMImage3DType, DCMImage3DType> ExtractFilterType;
//形态学细化类型
typedef itk::BinaryThinningImageFilter<DCMImage2DType, DCMImage2DType> BinaryThinningImageFilterType;
typedef itk::RescaleIntensityImageFilter<DCMImage2DType, PNGImageType> RescaleType;
//形态学开运算
typedef itk::BinaryBallStructuringElement<DCMPixelType, 2> StructuringElementType;
typedef itk::BinaryErodeImageFilter<DCMImage2DType, DCMImage2DType, StructuringElementType> ErodeFilterType;
typedef itk::BinaryDilateImageFilter<DCMImage2DType, DCMImage2DType, StructuringElementType> DilateFilterType;
typedef itk::ConnectedComponentImageFilter<DCMImage2DType, DCMImage2DType> ConnectedComponentFilter;
typedef itk::LabelShapeKeepNObjectsImageFilter<DCMImage2DType> LabelShapeKeepNObjectsFilter;
//全景图增强
typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<DCMImage2DType, DCMImage2DType> EnhancementFilterType;
//抽取感兴趣区域
typedef itk::RegionOfInterestImageFilter<DCMImage3DType, DCMImage3DType> ROIFilterType;
typedef itk::RegionOfInterestImageFilter<LabelImage3DType, LabelImage3DType> LabelFilterType;
typedef itk::SmoothingRecursiveGaussianImageFilter<DCMImage3DType, DCMImage3DType> Smoothing3DFilterType;
typedef itk::NumericSeriesFileNames SeriesFileNamesGeneratorType;

class Panorama
{
public:
  Panorama(std::string resultDir,
           std::string mipFile,
           std::string otsuFile,
           std::string morphFile,
           std::string thinFile,
           std::string splineFile,
           std::string panoramaFile,
           std::string looseROIFile,
           std::string boxROIFile,
           std::string toothMaskFile);
  ~Panorama();

  //像素点的坐标，小数，最后需要插值
  struct Point
  {
    double x = 0;
    double y = 0;
  };

  itk::SmartPointer<DCMImage2DType> Mip(itk::SmartPointer<DCMImage3DType> teeth, int lowerSlice, int higherSlice);
  itk::SmartPointer<DCMImage2DType> Otsu(itk::SmartPointer<DCMImage2DType> mip);
  itk::SmartPointer<DCMImage2DType> Morph(itk::SmartPointer<DCMImage2DType>, int dilateR, int erodeR);
  std::vector<Point> SplineFitting(itk::SmartPointer<DCMImage2DType> morph);
  DCMImage2DType::Pointer Arch(std::vector<Point> controlPoints,
                               DCMImage2DType::Pointer morph,
                               std::vector<Point> &result);
  DCMImage2DType::Pointer GeneratePanorama(DCMImage3DType::Pointer teeth,
                                           std::vector<Point> &result,
                                           int thickness,
                                           std::vector<std::vector<Panorama::Point>> &interpoMap);
  DCMImage3DType::Pointer LooseROI(DCMImage3DType::Pointer teeth,
                                   std::vector<Panorama::Point> boxPoints,
                                   std::vector<std::vector<Panorama::Point>> interpoMap);
  DCMImage3DType::Pointer BoxROI(DCMImage3DType::Pointer teeth,
                                 std::vector<Panorama::Point> boxPoints,
                                 std::vector<std::vector<Panorama::Point>> interpoMap,
                                 int leftPadding,
                                 int rightPadding,
                                 int upPadding,
                                 int downPadding);

private:
  //样条曲线基函数
  double BaseFunc(int i, int deg, double t, std::vector<double> knot)
  {
    if (deg == 0)
    {
      if (knot[i] <= t && t < knot[i + 1])
      {
        return 1;
      }
      return 0;
    }
    else
    {
      double result = 0;
      double base1 = BaseFunc(i, deg - 1, t, knot);
      double base2 = BaseFunc(i + 1, deg - 1, t, knot);
      result += base1 == 0 ? 0 : (t - knot[i]) == 0 ? 0 : (t - knot[i]) / (knot[i + deg] - knot[i]) * base1;
      result += base2 == 0                   ? 0 :
                (knot[i + deg + 1] - t) == 0 ? 0 :
                                               (knot[i + deg + 1] - t) / (knot[i + deg + 1] - knot[i + 1]) * base2;
      return result;
    }
  }

  // ImageIO不能重复使用
  void reset()
  {
    dcmIO = DCMImageIO::New();
    pngIO = PNGImageIO::New();
    niftiIO = NiftiImageIO::New();
  }

  //样条曲线拟合
  std::vector<Point> spline(std::vector<Point> controlPoints, const int degree, const double epsilon);

  //获得样条曲线的控制点
  std::vector<Point> getSplineControlPoints(std::vector<Point> pixels, int numberOfPoints);

  //右半牙弓线排序函数
  static bool rightCompare(Point &lhs, Point &rhs)
  {
    if (lhs.x == rhs.x)
      return lhs.y < rhs.y;
    else
      return lhs.x < rhs.x;
  }

  //左半牙弓线排序函数
  static bool leftCompare(Point &lhs, Point &rhs)
  {
    if (lhs.x == rhs.x)
      return lhs.y > rhs.y;
    else
      return lhs.x < rhs.x;
  }

  //双线性插值
  double interpo(DCMImage3DType::Pointer teeth, double x, double y, double z);

  //结果文件夹
  std::string resultDir;
  //最大强度投影文件
  std::string mipFile;
  //大津阈值分割文件
  std::string otsuFile;
  //形态学开运算文件
  std::string morphFile;
  //形态学细化文件
  std::string thinFile;
  //样条曲线
  std::string splineFile;
  //全景图文件
  std::string panoramaFile;
  // ROI文件
  std::string looseROIFile;
  // box文件
  std::string boxROIFile;
  //单颗牙齿mask文件
  std::string toothMaskFile;
  //细化牙弓线像素点个数
  int thinPixelNum = 500;

  ScalarImageToHistogramGeneratorType::Pointer scalarImageToHistogramGenerator;
  CalculatorType::Pointer calculator;
  BinaryFilterType::Pointer binaryFilter;
  seriesReaderType::Pointer seriesReader;
  seriesWriterType::Pointer seriesWriter;
  DCMImageIO::Pointer dcmIO;
  PNGImageIO::Pointer pngIO;
  NamesGeneratorType::Pointer namesGenerator;
  DilateFilterType::Pointer binaryDilate;
  ErodeFilterType::Pointer binaryErode;
  MipFilterType::Pointer mipFilter;
  ExtractFilterType::Pointer extractFilter;
  Smoothing2DFilterType::Pointer smoothing2DFilter;
  DCMReaderType::Pointer dcmReader;
  DCMWriterType::Pointer dcmWriter;
  BinaryThinningImageFilterType::Pointer binaryThinningImageFilter;
  RescaleType::Pointer rescaler;
  PNGWriterType::Pointer pngWriter;
  PNGReaderType::Pointer pngReader;
  ConnectedComponentFilter::Pointer connectedComponentFilter;
  LabelShapeKeepNObjectsFilter::Pointer labelShapeKeepNObjectsFilter;
  SeriesFileNamesGeneratorType::Pointer seriesFileNamesGenerator;
  NiftiImageIO::Pointer niftiIO;
  NiftiWriterType::Pointer niftiWriter;
  LabelWriterType::Pointer labelWriter;
};