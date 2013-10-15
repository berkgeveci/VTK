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
#include "vtkMergePoints.h"
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
  //mutable vtkThreadLocalObject<vtkNonMergingPointLocator> Locator;
  mutable vtkThreadLocalObject<vtkMergePoints> Locator;

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

      // vtkNonMergingPointLocator*& locator = this->Locator.Local();
      // locator->SetPoints(newPts);

      vtkMergePoints*& locator = this->Locator.Local();
      locator->InitPointInsertion (newPts,
                                   this->Input->GetBounds(),
                                   this->Input->GetNumberOfPoints());

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
      //cout << begin << " " << end << endl;

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

      //vtkNonMergingPointLocator* loc = this->Locator.Local();
      vtkMergePoints* loc = this->Locator.Local();

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

class vtkContourGridFunctor2 : public vtkInitializableFunctor
{
  vtkSMPContourGrid* Filter;

  vtkUnstructuredGrid* Input;
  vtkDataArray* InScalars;

  vtkMultiBlockDataSet* Output;

  double* Values;

  mutable vtkThreadLocal<std::vector<vtkPolyData*> > Outputs;

public:

  vtkContourGridFunctor2(vtkSMPContourGrid* filter,
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

  ~vtkContourGridFunctor2()
    {
    }

  void Initialize() const
    {
    }


  void operator()(vtkIdType begin, vtkIdType end) const
    {
      //cout << begin << " " << end << endl;
      vtkNew<vtkPolyData> output;

      vtkNew<vtkPoints> newPts;

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

      output->SetPoints(newPts.GetPointer());

      vtkIdType numCells = this->Input->GetNumberOfCells();

      vtkIdType estimatedSize=static_cast<vtkIdType>(
        pow(static_cast<double>(numCells),.75));
      estimatedSize = estimatedSize / 1024 * 1024; //multiple of 1024
      if (estimatedSize < 1024)
        {
        estimatedSize = 1024;
        }

      newPts->Allocate(estimatedSize, estimatedSize);

      // vtkNew<vtkNonMergingPointLocator> locator;
      // locator->SetPoints(newPts.GetPointer());

      vtkNew<vtkMergePoints> locator;
      locator->InitPointInsertion (newPts.GetPointer(),
                                   this->Input->GetBounds(),
                                   this->Input->GetNumberOfPoints());

      // vtkNew<vtkPointLocator> locator;
      // locator->InitPointInsertion (newPts.GetPointer(),
      //                              this->Input->GetBounds(),
      //                              this->Input->GetNumberOfPoints());

      vtkNew<vtkCellArray> newVerts;
      newVerts->Allocate(estimatedSize,estimatedSize);

      vtkNew<vtkCellArray> newLines;
      newLines->Allocate(estimatedSize,estimatedSize);

      vtkNew<vtkCellArray> newPolys;
      newPolys->Allocate(estimatedSize,estimatedSize);

      vtkSmartPointer<vtkDataArray> cellScalars;
      cellScalars.TakeReference(this->InScalars->NewInstance());
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

      vtkNew<vtkGenericCell> cell;

      for (vtkIdType i=begin; i<end; i++)
        {
        this->Input->GetCell(i, cell.GetPointer());

        vtkIdList* cellPts = cell->GetPointIds();
        this->InScalars->GetTuples(cellPts, cellScalars.GetPointer());

        cell->Contour(Values[0],
                      cellScalars.GetPointer(),
                      locator.GetPointer(),
                      newVerts.GetPointer(),
                      newLines.GetPointer(),
                      newPolys.GetPointer(),
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

        output->Squeeze();

      output->Register(0);
      this->Outputs.Local().push_back(output.GetPointer());
    }

  void Finalize()
    {
      vtkNew<vtkMultiPieceDataSet> mp;
      int count = 0, count2 = 0;
      vtkIdType totalNumCells = 0;

      vtkThreadLocal<std::vector<vtkPolyData*> >::iterator outIter =
        this->Outputs.begin();
      while(outIter != this->Outputs.end())
        {
        std::vector<vtkPolyData*>& outs = *outIter;
        std::vector<vtkPolyData*>::iterator iter = outs.begin();
        while (iter != outs.end())
          {
          //cout << count << " : " << (*iter)->GetNumberOfCells() << endl;
          mp->SetPiece(count++, *iter);
          totalNumCells += (*iter)->GetNumberOfCells();
          (*iter)->Delete();
          iter++;
          }
        ++outIter;
        }

      this->Output->SetBlock(0, mp.GetPointer());

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

  //vtkContourGridFunctor functor(this, input, inScalars, values, output);
  vtkContourGridFunctor2 functor(this, input, inScalars, values, output);
  // When using vtkContourGridFunctor2, it is crucial to set the grain
  // right. When the grain is too small, which tends to be the default,
  // the overhead of allocating data structures, building locators etc.
  // ends up being too big. When using vtkContourGridFunctor, it doesn't
  // matter as much.
  vtkParallelUtilities::ForEach(0, numCells, &functor, 50000);

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
