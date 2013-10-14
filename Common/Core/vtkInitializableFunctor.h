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
// .NAME vtkFunctor -
// .SECTION Description

#ifndef __vtkInitializableFunctor_h__
#define __vtkInitializableFunctor_h__

#include "vtkType.h"
#include "vtkCommonCoreModule.h" // For export macro
#include "vtkThreadLocal.h" // For Initialized

class VTKCOMMONCORE_EXPORT vtkInitializableFunctor
{
public:
  // Description:
  virtual void Initialize() const = 0;

  // Description:
  virtual void Finalize() {}

  // Description:
  // This is the execution operator that needs to be defined
  // by subclasses in order to implement an operation that can
  // be run in parallel. All implementations of this method
  // need to be thread safe.
  virtual void operator() (vtkIdType begin, vtkIdType end) const = 0;

  // Description:
  // Default construtor and destructor.
  vtkInitializableFunctor();
  virtual ~vtkInitializableFunctor();

  // Description:
  void Execute(vtkIdType begin, vtkIdType end) const;

private:
  mutable vtkThreadLocal<unsigned char> Initialized;

};

#endif
