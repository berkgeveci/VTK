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
#include "vtkInitializableFunctor.h"
#include "vtkThreadLocal.h"
#include "vtkThreadLocalObject.h"

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

class vtkContourGridFunctor : public vtkInitializableFunctor
{
  vtkSMPContourGrid* Filter;

  vtkUnstructuredGrid* Input;
  vtkDataArray* InScalars;

  double* Values;

  mutable vtkThreadLocal<vtkDataArray*> CellScalars;

  mutable vtkThreadLocalObject<vtkGenericCell> Cell;
  mutable vtkThreadLocalObject<vtkPoints> NewPts;
  mutable vtkThreadLocalObject<vtkCellArray> NewVerts;
  mutable vtkThreadLocalObject<vtkCellArray> NewLines;
  mutable vtkThreadLocalObject<vtkCellArray> NewPolys;
  mutable vtkThreadLocalObject<vtkPolyData> Output;
  mutable vtkThreadLocalObject<vtkNonMergingPointLocator> Locator;

public:

  vtkContourGridFunctor(vtkSMPContourGrid* filter,
                        vtkUnstructuredGrid* input,
                        vtkDataArray* inScalars,
                        double* values) : Filter(filter), Input(input), InScalars(inScalars), Values(values)
    {
    }

  ~vtkContourGridFunctor()
    {
      vtkThreadLocal<vtkDataArray*>::iterator cellScalarsIter =
        this->CellScalars.begin();
      while(cellScalarsIter != this->CellScalars.end())
        {
        (*cellScalarsIter)->Delete();
        ++cellScalarsIter;
        }
    }

  void Initialize() const
    {
      vtkPoints*& newPts = this->NewPts.Local();

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

      vtkNonMergingPointLocator*& locator = this->Locator.Local();
      locator->SetPoints(newPts);

      vtkCellArray*& newVerts = this->NewVerts.Local();
      newVerts->Allocate(estimatedSize,estimatedSize);

      vtkCellArray*& newLines = this->NewLines.Local();
      newLines->Allocate(estimatedSize,estimatedSize);

      vtkCellArray*& newPolys = this->NewPolys.Local();
      newPolys->Allocate(estimatedSize,estimatedSize);

      vtkDataArray*& cellScalars = this->CellScalars.Local();
      cellScalars = this->InScalars->NewInstance();
      cellScalars->SetNumberOfComponents(this->InScalars->GetNumberOfComponents());
      cellScalars->Allocate(VTK_CELL_SIZE*this->InScalars->GetNumberOfComponents());

      vtkPolyData*& output = this->Output.Local();

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
    }


  void operator()(vtkIdType begin, vtkIdType end) const
    {
      vtkGenericCell* cell = this->Cell.Local();
      vtkDataArray* cs = this->CellScalars.Local();
      vtkPointData* inPd = this->Input->GetPointData();
      vtkCellData* inCd = this->Input->GetCellData();

      vtkPolyData* output = this->Output.Local();
      vtkPointData* outPd = output->GetPointData();
      vtkCellData* outCd = output->GetCellData();

      vtkCellArray* vrts = this->NewVerts.Local();
      vtkCellArray* lines = this->NewLines.Local();
      vtkCellArray* polys = this->NewPolys.Local();

      vtkNonMergingPointLocator* loc = this->Locator.Local();

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

  void Finalize()
    {
      vtkThreadLocalObject<vtkPolyData>::iterator outIter =
        this->Output.begin();
      vtkThreadLocalObject<vtkCellArray>::iterator newVertsIter =
        this->NewVerts.begin();
      vtkThreadLocalObject<vtkCellArray>::iterator newLinesIter =
        this->NewLines.begin();
      vtkThreadLocalObject<vtkCellArray>::iterator newPolysIter =
        this->NewPolys.begin();
      while(outIter != this->Output.end())
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

        ++newVertsIter;
        ++newLinesIter;
        ++newPolysIter;
        ++outIter;
        }

      vtkIdType totalNumCells = 0;
      outIter = this->Output.begin();
      while(outIter != this->Output.end())
        {
        vtkPolyData* anOutput = *outIter;
        totalNumCells += anOutput->GetNumberOfCells();
        ++outIter;
        }
      cout << "Total number of cells: " << totalNumCells << endl;

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

  return 1;
}

void vtkSMPContourGrid::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
