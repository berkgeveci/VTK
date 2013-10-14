/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSMPContourGrid.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSMPContourGrid -

#ifndef __vtkSMPContourGrid_h
#define __vtkSMPContourGrid_h

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkContourGrid.h"

class VTKFILTERSCORE_EXPORT vtkSMPContourGrid : public vtkContourGrid
{
public:
  vtkTypeMacro(vtkSMPContourGrid,vtkContourGrid);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct object with initial range (0,1) and single contour value
  // of 0.0.
  static vtkSMPContourGrid *New();

protected:
  vtkSMPContourGrid();
  ~vtkSMPContourGrid();

  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  virtual int FillOutputPortInformation(int port, vtkInformation* info);

private:
  vtkSMPContourGrid(const vtkSMPContourGrid&);  // Not implemented.
  void operator=(const vtkSMPContourGrid&);  // Not implemented.
};

#endif
