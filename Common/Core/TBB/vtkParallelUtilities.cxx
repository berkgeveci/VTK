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

#include "vtkCriticalSection.h"
#include "vtkFunctor.h"
#include "vtkObjectFactory.h"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

vtkStandardNewMacro(vtkParallelUtilities);

struct vtkParallelUtilitiesInit
{
  tbb::task_scheduler_init Init;

  vtkParallelUtilitiesInit(int numThreads) : Init(numThreads)
    {
    }
};

static bool vtkParallelUtilitiesInitialized = 0;
static vtkSimpleCriticalSection vtkParallelUtilitiesCS;


namespace
{
class FuncCall
{
  const vtkFunctor* o;

public:
  void operator() (const tbb::blocked_range<vtkIdType>& r) const
    {
      (*o)(r.begin(), r.end());
    }

  FuncCall (const vtkFunctor* _o) : o(_o)
    {
    }
  ~FuncCall ()
    {
    }
};
}

//--------------------------------------------------------------------------------
vtkParallelUtilities::vtkParallelUtilities()
{
}

//--------------------------------------------------------------------------------
vtkParallelUtilities::~vtkParallelUtilities()
{
}

//--------------------------------------------------------------------------------
void vtkParallelUtilities::Initialize(int numThreads)
{
  vtkParallelUtilitiesCS.Lock();
  if (!vtkParallelUtilitiesInitialized)
    {
    static vtkParallelUtilitiesInit aInit(numThreads);
    vtkParallelUtilitiesInitialized = true;
    }
  vtkParallelUtilitiesCS.Unlock();
}

//--------------------------------------------------------------------------------
void vtkParallelUtilities::ForEach(vtkIdType first,
                                   vtkIdType last,
                                   const vtkFunctor* op,
                                   int grain)
{
  vtkIdType n = last - first;
  if (!n)
    {
    return;
    }
  if (grain > 0)
    {
    tbb::parallel_for(tbb::blocked_range<vtkIdType>(first, last, grain), FuncCall(op));
    }
  else
    {
    tbb::parallel_for(tbb::blocked_range<vtkIdType>(first, last), FuncCall(op));
    }
}

//--------------------------------------------------------------------------------
void vtkParallelUtilities::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf( os, indent );
}
