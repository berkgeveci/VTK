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
#include "vtkInitializableFunctor.h"
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

template <typename T>
void vtkParallelUtilitiesForEach(vtkIdType first,
                                 vtkIdType last,
                                 const T* op,
                                 int grain)
{
  vtkIdType n = last - first;
  if (!n)
    {
    return;
    }

  if (grain == 0 || grain >= n)
    {
    op->Execute(first, last);
    }
  else
    {
    vtkIdType b = first;
    while (b < last)
      {
      vtkIdType e = b + grain;
      if (e > last)
        {
        e = last;
        }
      op->Execute(b, e);
      b = e;
      }
    }
}

//--------------------------------------------------------------------------------
void vtkParallelUtilities::ForEach(vtkIdType first,
                                   vtkIdType last,
                                   const vtkFunctor* op,
                                   int grain)
{
  vtkParallelUtilitiesForEach(first, last, op, grain);
}

//--------------------------------------------------------------------------------
void vtkParallelUtilities::ForEach(vtkIdType first,
                                   vtkIdType last,
                                   vtkInitializableFunctor* op,
                                   int grain)
{
  vtkParallelUtilitiesForEach(first, last, op, grain);
  op->Finalize();
}

//--------------------------------------------------------------------------------
void vtkParallelUtilities::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf( os, indent );
}
