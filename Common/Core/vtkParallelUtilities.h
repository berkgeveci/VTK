/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkParallelUtilities.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkParallelUtilities - A set of parallel (multi-threaded) utility functions.
// .SECTION Description
// vtkParallelUtilities provides a set of utility functions that can
// be used to parallelize parts of VTK code using multiple threads.
// There are several back-end implementations of parallel functionality
// (currently Sequential, TBB and X-Kaapi) that actual execution is
// delegated to.

#ifndef __vtkParallelUtilities_h__
#define __vtkParallelUtilities_h__

#include "vtkCommonCoreModule.h" // For export macro
#include "vtkObject.h"

class vtkFunctor;
class vtkInitializableFunctor;

class VTKCOMMONCORE_EXPORT vtkParallelUtilities : public vtkObject
{
public:
  vtkTypeMacro(vtkParallelUtilities,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkParallelUtilities* New();

  // Description:
  // Execute a for each operation in parallel. First and last
  // define the range over which to operate (which is defined
  // by the operator). The operation executed is defined by
  // operator() of the op object. The grain gives the parallel
  // engine a hint about the coarseness over which to parallelize
  // the function (as defined by last-first of each execution of
  // operator() ).
  static void ForEach(vtkIdType first,
                      vtkIdType last,
                      const vtkFunctor* op,
                      int grain=0);

  // Description:
  static void ForEach(vtkIdType first,
                      vtkIdType last,
                      vtkInitializableFunctor* op,
                      int grain=0);

  // Description:
  // Initialize the underlying libraries for execution. This is
  // not required as it is automatically called before the first
  // execution of any parallel code. However, it can be used to
  // control the maximum number of threads used when the back-end
  // supports it (currently TBB only). Make sure to call it before
  // any other parallel operation.
  static void Initialize(int numThreads=-1);

protected:
  vtkParallelUtilities();
  ~vtkParallelUtilities();

private:
  vtkParallelUtilities(const vtkParallelUtilities&);  // Not implemented.
  void operator=(const vtkParallelUtilities&);  // Not implemented.
};

#endif
