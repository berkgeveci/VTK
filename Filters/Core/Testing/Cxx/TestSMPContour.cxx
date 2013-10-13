/*=========================================================================

  Program:   Visualization Toolkit
  Module:    TestCutter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkNew.h"
#include "vtkRTAnalyticSource.h"
#include "vtkPolyData.h"
#include "vtkDataSetTriangleFilter.h"
#include "vtkSMPContourGrid.h"
#include "vtkContourGrid.h"
#include "vtkContourFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkTimerLog.h"
#include "vtkNonMergingPointLocator.h"
#include "vtkParallelUtilities.h"

const int EXTENT = 50;
int TestSMPContour(int, char *[])
{
  //vtkParallelUtilities::Initialize(1);

  vtkNew<vtkTimerLog> tl;

  vtkNew<vtkRTAnalyticSource> imageSource;
  imageSource->SetWholeExtent(-EXTENT, EXTENT, -EXTENT, EXTENT, -EXTENT, EXTENT);

  vtkNew<vtkDataSetTriangleFilter> tetraFilter;
  tetraFilter->SetInputConnection(imageSource->GetOutputPort());

  tl->StartTimer();
  tetraFilter->Update();
  tl->StopTimer();
  cout << "Tetrahedralize time: " << tl->GetElapsedTime() << endl;

  vtkNew<vtkContourGrid> cg;
  cg->SetInputData(tetraFilter->GetOutput());
  cg->SetInputArrayToProcess(0, 0, 0, 0, "RTData");
  cg->SetValue(0, 200);
  tl->StartTimer();
  cg->Update();
  tl->StopTimer();

  cout << "Contour grid: " << endl;
  cout << "Number of cells: " << cg->GetOutput()->GetNumberOfCells() << endl;
  cout << "Time: " << tl->GetElapsedTime() << endl;

  vtkNew<vtkContourFilter> cf;
  vtkNew<vtkNonMergingPointLocator> nmpl;
  cf->SetLocator(nmpl.GetPointer());
  cf->SetInputData(tetraFilter->GetOutput());
  cf->SetInputArrayToProcess(0, 0, 0, 0, "RTData");
  cf->SetValue(0, 200);
  tl->StartTimer();
  cf->Update();
  tl->StopTimer();

  cout << "Contour filter: " << endl;
  cout << "Number of cells: " << cf->GetOutput()->GetNumberOfCells() << endl;
  cout << "Time: " << tl->GetElapsedTime() << endl;

  vtkNew<vtkSMPContourGrid> cg2;
  cg2->SetInputData(tetraFilter->GetOutput());
  cg2->SetInputArrayToProcess(0, 0, 0, 0, "RTData");
  cg2->SetValue(0, 200);
  tl->StartTimer();
  cg2->Update();
  tl->StopTimer();

  cout << "SMP Contour grid: " << endl;
  cout << "Number of cells: " << cg2->GetOutput()->GetNumberOfCells() << endl;
  cout << "Time: " << tl->GetElapsedTime() << endl;

  return EXIT_SUCCESS;
}
