
#include "PointsSource.h"


#include "vtkPointSource.h"

#include "vtkCellArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"

#include <iostream>
#include <float.h>
#include <math.h>

vtkStandardNewMacro(PointSource);

//----------------------------------------------------------------------------
PointSource::PointSource(vtkIdType numPts)
{
  this->NumberOfPoints = (numPts > 0 ? numPts : 10);
  this->SetNumberOfInputPorts(0);
  this->points = vtkSmartPointer<vtkPoints>::New();
};

void PointSource::SetPoints(vtkSmartPointer<vtkPoints> *pts, int npts)
{
    vtkIdType ni = 0;
    double x[3];
    for(vtkIdType i = 0; i < npts; i++)
    {
        ni = pts[i]->GetNumberOfPoints();
        for(vtkIdType j = 0; j < ni; j++)
        {
            pts[i]->GetPoint(j, x);
            this->points->InsertNextPoint(x);
        }
    }
    this->NumberOfPoints = this->points->GetNumberOfPoints();
}

//----------------------------------------------------------------------------
int PointSource::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the ouptut
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkIdType i;
//  double theta, rho, cosphi, sinphi, radius;
  double x[3];
  vtkPoints *newPoints;
  vtkCellArray *newVerts;

  newPoints = vtkPoints::New();
  newPoints->Allocate(this->NumberOfPoints);
  newVerts = vtkCellArray::New();
  newVerts->Allocate(newVerts->EstimateSize(1,this->NumberOfPoints));

  newVerts->InsertNextCell(this->NumberOfPoints);
    for(i = 0; i < points->GetNumberOfPoints(); i++)
    {
        points->GetPoint(i, x);
        newVerts->InsertCellPoint(newPoints->InsertNextPoint(x));
    }

/*
  if (this->Distribution == VTK_POINT_SHELL)
    {  // only produce points on the surface of the sphere
    for (i=0; i<this->NumberOfPoints; i++)
      {
      cosphi = 1 - 2*vtkMath::Random();
      sinphi = sqrt(1 - cosphi*cosphi);
      radius = this->Radius * sinphi;
      theta = 2.0 * vtkMath::Pi() * vtkMath::Random();
      x[0] = this->Center[0] + radius*cos(theta);
      x[1] = this->Center[1] + radius*sin(theta);
      x[2] = this->Center[2] + this->Radius*cosphi;
      newVerts->InsertCellPoint(newPoints->InsertNextPoint(x));
      }
    }
  else
    { // uniform distribution throughout the sphere volume
    for (i=0; i<this->NumberOfPoints; i++)
      {
      cosphi = 1 - 2*vtkMath::Random();
      sinphi = sqrt(1 - cosphi*cosphi);
      rho = this->Radius*pow(vtkMath::Random(),0.33333333);
      radius = rho * sinphi;
      theta = 2.0 * vtkMath::Pi() * vtkMath::Random();
      x[0] = this->Center[0] + radius*cos(theta);
      x[1] = this->Center[1] + radius*sin(theta);
      x[2] = this->Center[2] + rho*cosphi;
      newVerts->InsertCellPoint(newPoints->InsertNextPoint(x));
      }
    }
  */
   //
   // Update ourselves and release memory
   //

  output->SetPoints(newPoints);
  newPoints->Delete();

  output->SetVerts(newVerts);
  newVerts->Delete();

  return 1;
}
