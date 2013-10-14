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
#include "vtkMultiBlockDataSet.h"
#include "vtkMultiPieceDataSet.h"
#include "vtkParallelUtilities.h"
#include "vtkInitializableFunctor.h"
#include "vtkThreadLocal.h"
#include "vtkThreadLocalObject.h"
#include "vtkInformation.h"

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

  vtkMultiBlockDataSet* Output;

  double* Values;

  mutable vtkThreadLocal<vtkDataArray*> CellScalars;

  mutable vtkThreadLocalObject<vtkGenericCell> Cell;
  mutable vtkThreadLocalObject<vtkPoints> NewPts;
  mutable vtkThreadLocalObject<vtkCellArray> NewVerts;
  mutable vtkThreadLocalObject<vtkCellArray> NewLines;
  mutable vtkThreadLocalObject<vtkCellArray> NewPolys;
  mutable vtkThreadLocalObject<vtkPolyData> Outputs;
  mutable vtkThreadLocalObject<vtkNonMergingPointLocator> Locator;
  //mutable vtkThreadLocalObject<vtkPointLocator> Locator;

public:

  vtkContourGridFunctor(vtkSMPContourGrid* filter,
                        vtkUnstructuredGrid* input,
                        vtkDataArray* inScalars,
                        double* values,
                        vtkMultiBlockDataSet* output) : Filter(filter),
                                                        Input(input),
                                                        InScalars(inScalars),
                                                        Values(values),
                                                        Output(output)
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
      vtkPolyData*& output = this->Outputs.Local();

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

      output->SetPoints(newPts);

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

      // vtkPointLocator*& locator = this->Locator.Local();
      // locator->InitPointInsertion (newPts,
      //                              this->Input->GetBounds(),
      //                              this->Input->GetNumberOfPoints());

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

      vtkPolyData* output = this->Outputs.Local();
      vtkPointData* outPd = output->GetPointData();
      vtkCellData* outCd = output->GetCellData();

      vtkCellArray* vrts = this->NewVerts.Local();
      vtkCellArray* lines = this->NewLines.Local();
      vtkCellArray* polys = this->NewPolys.Local();

      vtkNonMergingPointLocator* loc = this->Locator.Local();
      //vtkPointLocator* loc = this->Locator.Local();

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
      vtkNew<vtkMultiPieceDataSet> mp;
      int count = 0;

      vtkThreadLocalObject<vtkPolyData>::iterator outIter =
        this->Outputs.begin();
      vtkThreadLocalObject<vtkCellArray>::iterator newVertsIter =
        this->NewVerts.begin();
      vtkThreadLocalObject<vtkCellArray>::iterator newLinesIter =
        this->NewLines.begin();
      vtkThreadLocalObject<vtkCellArray>::iterator newPolysIter =
        this->NewPolys.begin();
      while(outIter != this->Outputs.end())
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

        mp->SetPiece(count++, output);

        ++newVertsIter;
        ++newLinesIter;
        ++newPolysIter;
        ++outIter;
        }

      this->Output->SetBlock(0, mp.GetPointer());

      vtkIdType totalNumCells = 0;
      outIter = this->Outputs.begin();
      while(outIter != this->Outputs.end())
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
  vtkMultiBlockDataSet *output = vtkMultiBlockDataSet::GetData(outputVector);

  // Not thread safe so calculate first.
  input->GetBounds();

  vtkPointData* inPd = input->GetPointData();
  vtkCellData* inCd = input->GetCellData();

  vtkDataArray* inScalars = this->GetInputArrayToProcess(0,inputVector);

  double *values=this->ContourValues->GetValues();

  vtkIdType numCells = input->GetNumberOfCells();

  vtkContourGridFunctor functor(this, input, inScalars, values, output);
  vtkParallelUtilities::ForEach(0, numCells, &functor);

  return 1;
}

int vtkSMPContourGrid::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // now add our info
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  return 1;
}

void vtkSMPContourGrid::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
