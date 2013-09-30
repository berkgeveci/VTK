/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkFunctor.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkParallelUtilities.h"

#include "vtkFunctor.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkParallelUtilities);

// Simple implementation that runs everything sequentially.

//--------------------------------------------------------------------------------
vtkParallelUtilities::vtkParallelUtilities()
{
}

//--------------------------------------------------------------------------------
vtkParallelUtilities::~vtkParallelUtilities()
{
}

//--------------------------------------------------------------------------------
void vtkParallelUtilities::Initialize(int)
{
}

//--------------------------------------------------------------------------------
void vtkParallelUtilities::ForEach(vtkIdType first,
                                   vtkIdType last,
                                   const vtkFunctor* op,
                                   int)
{
  vtkIdType n = last - first;
  if (!n)
    {
    return;
    }

  (*op)(first, last);
}

//--------------------------------------------------------------------------------
void vtkParallelUtilities::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf( os, indent );
}
