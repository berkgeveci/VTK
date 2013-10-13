/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkContourGrid.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSMPContourGrid.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkGenericCell.h"
#include "vtkNew.h"
#include "vtkNonMergingPointLocator.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPointLocator.h"

#include "vtkParallelUtilities.h"
#include "vtkFunctor.h"
#include "vtkThreadLocal.h"

#include <math.h>

vtkStandardNewMacro(vtkSMPContourGrid);

// Construct object with initial range (0,1) and single contour value
// of 0.0.
vtkSMPContourGrid::vtkSMPContourGrid()
{
}

vtkSMPContourGrid::~vtkSMPContourGrid()
{
}

namespace
{
vtkThreadLocal<bool> Initialized;

vtkThreadLocal<vtkGenericCell*> Cell;
vtkThreadLocal<vtkPoints*> NewPts;
vtkThreadLocal<vtkCellArray*> NewVerts;
vtkThreadLocal<vtkCellArray*> NewLines;
vtkThreadLocal<vtkCellArray*> NewPolys;
vtkThreadLocal<vtkDataArray*> CellScalars;
vtkThreadLocal<vtkPolyData*> Output;
vtkThreadLocal<vtkNonMergingPointLocator*> Locator;
}

class vtkContourGridFunctor : public vtkFunctor
{
  vtkSMPContourGrid* Filter;

  vtkUnstructuredGrid* Input;
  vtkDataArray* InScalars;

  double* Values;

  void Initialize() const
    {
      Cell.Local() = vtkGenericCell::New();

      vtkPoints*& newPts = NewPts.Local();
      newPts = vtkPoints::New();

      // set precision for the points in the output
      if(this->Filter->GetOutputPointsPrecision() == vtkAlgorithm::DEFAULT_PRECISION)
        {
        newPts->SetDataType(this->Input->GetPoints()->GetDataType());
        }
      else if(this->Filter->GetOutputPointsPrecision() == vtkAlgorithm::SINGLE_PRECISION)
        {
        newPts->SetDataType(VTK_FLOAT);
        }
      else if(this->Filter->GetOutputPointsPrecision() == vtkAlgorithm::DOUBLE_PRECISION)
        {
        newPts->SetDataType(VTK_DOUBLE);
        }

      vtkIdType numCells = this->Input->GetNumberOfCells();

      vtkIdType estimatedSize=static_cast<vtkIdType>(
        pow(static_cast<double>(numCells),.75));
      estimatedSize = estimatedSize / 1024 * 1024; //multiple of 1024
      if (estimatedSize < 1024)
        {
        estimatedSize = 1024;
        }

      newPts->Allocate(estimatedSize, estimatedSize);

      vtkNonMergingPointLocator*& locator = Locator.Local();
      locator = vtkNonMergingPointLocator::New();
      locator->SetPoints(newPts);

      vtkCellArray*& newVerts = NewVerts.Local();
      newVerts = vtkCellArray::New();
      newVerts->Allocate(estimatedSize,estimatedSize);

      vtkCellArray*& newLines = NewLines.Local();
      newLines = vtkCellArray::New();
      newLines->Allocate(estimatedSize,estimatedSize);

      vtkCellArray*& newPolys = NewPolys.Local();
      newPolys = vtkCellArray::New();
      newPolys->Allocate(estimatedSize,estimatedSize);

      vtkDataArray*& cellScalars = CellScalars.Local();
      cellScalars = this->InScalars->NewInstance();
      cellScalars->SetNumberOfComponents(this->InScalars->GetNumberOfComponents());
      cellScalars->Allocate(VTK_CELL_SIZE*this->InScalars->GetNumberOfComponents());

      vtkPolyData*& output = Output.Local();
      output = vtkPolyData::New();

      vtkPointData* outPd = output->GetPointData();
      vtkCellData* outCd = output->GetCellData();
      if (!this->Filter->GetComputeNormals())
        {
        outPd->CopyScalarsOff();
        }
      vtkPointData* inPd = this->Input->GetPointData();
      vtkCellData* inCd = this->Input->GetCellData();
      outPd->InterpolateAllocate(inPd, estimatedSize, estimatedSize);
      outCd->CopyAllocate(inCd, estimatedSize, estimatedSize);

      Initialized.Local() = true;
    }

public:

  vtkContourGridFunctor(vtkSMPContourGrid* filter,
                        vtkUnstructuredGrid* input,
                        vtkDataArray* inScalars,
                        double* values) : Filter(filter), Input(input), InScalars(inScalars), Values(values)
    {
    }

  ~vtkContourGridFunctor()
    {
      vtkThreadLocal<vtkGenericCell*>::iterator cellIter =
        Cell.begin();
      while(cellIter != Cell.end())
        {
        (*cellIter)->Delete();
        cellIter++;
        }
      vtkThreadLocal<vtkPoints*>::iterator newPtsIter =
        NewPts.begin();
      while(newPtsIter != NewPts.end())
        {
        (*newPtsIter)->Delete();
        newPtsIter++;
        }
      vtkThreadLocal<vtkDataArray*>::iterator cellScalarsIter =
        CellScalars.begin();
      while(cellScalarsIter != CellScalars.end())
        {
        (*cellScalarsIter)->Delete();
        cellScalarsIter++;
        }
      vtkThreadLocal<vtkNonMergingPointLocator*>::iterator locIter =
        Locator.begin();
      while(locIter != Locator.end())
        {
        (*locIter)->Delete();
        locIter++;
        }
    }

  void operator()(vtkIdType begin, vtkIdType end) const
    {
      if (!Initialized.Local())
        {
        this->Initialize();
        }

      vtkGenericCell* cell = Cell.Local();
      vtkDataArray* cs = CellScalars.Local();
      vtkPointData* inPd = this->Input->GetPointData();
      vtkCellData* inCd = this->Input->GetCellData();

      vtkPolyData* output = Output.Local();
      vtkPointData* outPd = output->GetPointData();
      vtkCellData* outCd = output->GetCellData();

      vtkCellArray* vrts = NewVerts.Local();
      vtkCellArray* lines = NewLines.Local();
      vtkCellArray* polys = NewPolys.Local();

      vtkNonMergingPointLocator* loc = Locator.Local();

      for (vtkIdType i=begin; i<end; i++)
        {
        this->Input->GetCell(i, cell);

        vtkIdList* cellPts = cell->GetPointIds();
        this->InScalars->GetTuples(cellPts, cs);

        cell->Contour(Values[0],
                      cs,
                      loc,
                      vrts,
                      lines,
                      polys,
                      inPd,
                      outPd,
                      inCd,
                      i,
                      outCd);
        }
    }
};

//
// Contouring filter for unstructured grids.
//
int vtkSMPContourGrid::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the input and output
  vtkUnstructuredGrid *input = vtkUnstructuredGrid::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  vtkPointData* inPd = input->GetPointData();
  vtkPointData* outPd = output->GetPointData();
  vtkCellData* inCd = input->GetCellData();
  vtkCellData* outCd = output->GetCellData();

  vtkDataArray* inScalars = this->GetInputArrayToProcess(0,inputVector);

  double *values=this->ContourValues->GetValues();

  vtkIdType numCells = input->GetNumberOfCells();

  vtkContourGridFunctor functor(this, input, inScalars, values);
  vtkParallelUtilities::ForEach(0, numCells, &functor);

  vtkThreadLocal<vtkPolyData*>::iterator outIter =
    Output.begin();
  vtkThreadLocal<vtkCellArray*>::iterator newVertsIter =
    NewVerts.begin();
  vtkThreadLocal<vtkCellArray*>::iterator newLinesIter =
    NewLines.begin();
  vtkThreadLocal<vtkCellArray*>::iterator newPolysIter =
    NewPolys.begin();
  while(outIter != Output.end())
    {
    vtkPolyData* output = *outIter;
    vtkCellArray* newVerts = *newVertsIter;
    vtkCellArray* newLines = *newLinesIter;
    vtkCellArray* newPolys = *newPolysIter;

    if (newVerts->GetNumberOfCells())
      {
      output->SetVerts(newVerts);
      }

    if (newLines->GetNumberOfCells())
      {
      output->SetLines(newLines);
      }

    if (newPolys->GetNumberOfCells())
      {
      output->SetPolys(newPolys);
      }

    output->Squeeze();

    newVerts->Delete();
    newLines->Delete();
    newPolys->Delete();
    newVertsIter++;
    newLinesIter++;
    newPolysIter++;
    outIter++;
    }

  vtkIdType totalNumCells = 0;
  outIter = Output.begin();
  while(outIter != Output.end())
    {
    vtkPolyData* anOutput = *outIter;
    totalNumCells += anOutput->GetNumberOfCells();
    anOutput->Delete();
    outIter++;
    }
  cout << "Total number of cells: " << totalNumCells << endl;

  return 1;

  //
  // Create objects to hold output of contour operation. First estimate
  // allocation size.
  //
  vtkIdType estimatedSize=static_cast<vtkIdType>(
    pow(static_cast<double>(numCells),.75));
  estimatedSize = estimatedSize / 1024 * 1024; //multiple of 1024
  if (estimatedSize < 1024)
    {
    estimatedSize = 1024;
    }

  if (!this->ComputeNormals)
    {
    outPd->CopyScalarsOff();
    }
  outPd->InterpolateAllocate(inPd, estimatedSize, estimatedSize);
  outCd->CopyAllocate(inCd, estimatedSize, estimatedSize);


  vtkNew<vtkPoints> newPts;

  // set precision for the points in the output
  if(this->GetOutputPointsPrecision() == vtkAlgorithm::DEFAULT_PRECISION)
    {
    newPts->SetDataType(input->GetPoints()->GetDataType());
    }
  else if(this->GetOutputPointsPrecision() == vtkAlgorithm::SINGLE_PRECISION)
    {
    newPts->SetDataType(VTK_FLOAT);
    }
  else if(this->GetOutputPointsPrecision() == vtkAlgorithm::DOUBLE_PRECISION)
    {
    newPts->SetDataType(VTK_DOUBLE);
    }

  newPts->Allocate(estimatedSize, estimatedSize);

  vtkNew<vtkCellArray> newVerts;
  newVerts->Allocate(estimatedSize,estimatedSize);

  vtkNew<vtkCellArray> newLines;
  newLines->Allocate(estimatedSize,estimatedSize);

  vtkNew<vtkCellArray> newPolys;
  newPolys->Allocate(estimatedSize,estimatedSize);

  vtkSmartPointer<vtkDataArray> cellScalars;
  cellScalars.TakeReference(inScalars->NewInstance());
  cellScalars->SetNumberOfComponents(inScalars->GetNumberOfComponents());
  cellScalars->Allocate(VTK_CELL_SIZE*inScalars->GetNumberOfComponents());

  vtkNew<vtkNonMergingPointLocator> locator;
  locator->SetPoints(newPts.GetPointer());

   /*
   vtkNew<vtkPointLocator> locator;
   locator->InitPointInsertion (newPts.GetPointer(),
                                input->GetBounds(),
                                input->GetNumberOfPoints());
   */


   vtkNew<vtkGenericCell> genericCell;

   vtkDataArray* cs = cellScalars.GetPointer();
   vtkPointLocator* loc = locator.GetPointer();
   vtkCellArray* vrts = newVerts.GetPointer();
   vtkCellArray* lines = newLines.GetPointer();
   vtkCellArray* polys = newPolys.GetPointer();
   vtkGenericCell* cell = genericCell.GetPointer();

   for (vtkIdType i=0; i<numCells; i++)
     {
     input->GetCell(i, cell);

     vtkIdList* cellPts = cell->GetPointIds();
     inScalars->GetTuples(cellPts, cs);

     cell->Contour(values[0],
                   cs,
                   loc,
                   vrts,
                   lines,
                   polys,
                   inPd,
                   outPd,
                   inCd,
                   i,
                   outCd);
     }

  if (newVerts->GetNumberOfCells())
    {
    output->SetVerts(newVerts.GetPointer());
    }

  if (newLines->GetNumberOfCells())
    {
    output->SetLines(newLines.GetPointer());
    }

  if (newPolys->GetNumberOfCells())
    {
    output->SetPolys(newPolys.GetPointer());
    }

  locator->Initialize();//releases leftover memory
  output->Squeeze();

  return 1;
}

void vtkSMPContourGrid::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
