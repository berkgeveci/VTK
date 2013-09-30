/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkFunctor.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkFunctor - A simple base-class that defined operator (). Used by parallel execution code.
// .SECTION Description
// vtkFunctor is a simple base-class that can be subclassed to defined
// an execution operator that can process a range of elements. This
// class is used by vtkParallelUtilities to process elements in parallel.
// Note that the parallel execution functions do not create any copies
// of their functor arguments. So it is up to the developer to make sure
// that functors do not change their state when execution operator(). If
// you need to store execution specific data, consider using vtkThreadLocal.

#ifndef __vtkFunctor_h__
#define __vtkFunctor_h__

#include "vtkType.h"
#include "vtkCommonCoreModule.h" // For export macro

class VTKCOMMONCORE_EXPORT vtkFunctor
{
public:
  // Description:
  // This is the execution operator that needs to be defined
  // by subclasses in order to implement an operation that can
  // be run in parallel. All implementations of this method
  // need to be thread safe.
  virtual void operator() (vtkIdType begin, vtkIdType end) const = 0;

  // Description:
  // Default construtor and destructor.
  vtkFunctor();
  virtual ~vtkFunctor();
};

#endif
