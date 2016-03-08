#ifndef POINTSSOURCE_H
#define POINTSSOURCE_H

#include "vtkFiltersSourcesModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"
#include <vtkSmartPointer.h>

#define VTK_POINT_UNIFORM   1
#define VTK_POINT_SHELL     0

#ifndef _WIN32
class VTKFILTERSSOURCES_EXPORT PointSource : public vtkPolyDataAlgorithm
#else
class PointSource : public vtkPolyDataAlgorithm
#endif
{
public:
  static PointSource *New();
//  vtkTypeMacro(PointSource,vtkPolyDataAlgorithm);
//  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the number of points to generate.
//  vtkSetClampMacro(NumberOfPoints,vtkIdType,1,VTK_ID_MAX);
//  vtkGetMacro(NumberOfPoints,vtkIdType);

  // Description:
  // Set the center of the point cloud.
//  vtkSetVector3Macro(Center,double);
//  vtkGetVectorMacro(Center,double,3);

  // Description:
  // Set the radius of the point cloud.  If you are
  // generating a Gaussian distribution, then this is
  // the standard deviation for each of x, y, and z.
//  vtkSetClampMacro(Radius,double,0.0,VTK_DOUBLE_MAX);
//  vtkGetMacro(Radius,double);

  // Description:
  // Specify the distribution to use.  The default is a
  // uniform distribution.  The shell distribution produces
  // random points on the surface of the sphere, none in the interior.
//  vtkSetMacro(Distribution,int);
//  void SetDistributionToUniform() {
//    this->SetDistribution(VTK_POINT_UNIFORM);};
//  void SetDistributionToShell() {
//    this->SetDistribution(VTK_POINT_SHELL);};
//  vtkGetMacro(Distribution,int);

  void SetPoints(vtkSmartPointer<vtkPoints> *pts, int npts);
protected:
  PointSource(vtkIdType numPts=10);
  ~PointSource() {};

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  vtkIdType NumberOfPoints;
  vtkSmartPointer<vtkPoints> points;

//private:
//  PointSource(const PointSource&);  // Not implemented.
//  void operator=(const PointSource&);  // Not implemented.
};

#endif // POINTSSOURCE_H
