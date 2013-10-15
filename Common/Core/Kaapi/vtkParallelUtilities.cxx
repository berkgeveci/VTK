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
#include "vtkInitializableFunctor.h"
#include "vtkObjectFactory.h"

#include <kaapic.h>

vtkStandardNewMacro(vtkParallelUtilities);

struct vtkParallelUtilitiesInit
{
  vtkParallelUtilitiesInit()
    {
      kaapic_init(KAAPIC_START_ONLY_MAIN);
    }

  ~vtkParallelUtilitiesInit()
    {
      kaapic_finalize();
    }
};

static bool vtkParallelUtilitiesInitialized = 0;
static vtkSimpleCriticalSection vtkParallelUtilitiesCS;

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
  vtkParallelUtilitiesCS.Lock();
  if (!vtkParallelUtilitiesInitialized)
    {
    static vtkParallelUtilitiesInit aInit;
    vtkParallelUtilitiesInitialized = true;
    }
  vtkParallelUtilitiesCS.Unlock();
}

namespace
{
inline void vtkParallelUtilitiesDoFor1(int32_t b, int32_t e, int32_t, const vtkFunctor* o )
{
  o->Execute(b, e);
}

inline void vtkParallelUtilitiesDoFor2(int32_t b, int32_t e, int32_t, const vtkInitializableFunctor* o )
{
  o->Execute(b, e);
}

}

//--------------------------------------------------------------------------------
void vtkParallelUtilities::ForEach(vtkIdType first,
                                   vtkIdType last,
                                   const vtkFunctor* op,
                                   int grain)
{
  vtkParallelUtilities::Initialize();

  vtkIdType n = last - first;
  if (!n)
    {
    return;
    }

  vtkIdType g = grain ? grain : sqrt(n);

  kaapic_begin_parallel(KAAPIC_FLAG_DEFAULT);
  kaapic_foreach_attr_t attr;
  kaapic_foreach_attr_init(&attr);
  kaapic_foreach_attr_set_grains(&attr, g, g);
  kaapic_foreach( first, last, &attr, 1, vtkParallelUtilitiesDoFor1, op );
  kaapic_end_parallel(KAAPIC_FLAG_DEFAULT);
  kaapic_foreach_attr_destroy(&attr);
}

//--------------------------------------------------------------------------------
void vtkParallelUtilities::ForEach(vtkIdType first,
                                   vtkIdType last,
                                   vtkInitializableFunctor* op,
                                   int grain)
{
  vtkParallelUtilities::Initialize();

  vtkIdType n = last - first;
  if (!n)
    {
    return;
    }

  vtkIdType g = grain ? grain : sqrt(n);

  kaapic_begin_parallel(KAAPIC_FLAG_DEFAULT);
  kaapic_foreach_attr_t attr;
  kaapic_foreach_attr_init(&attr);
  kaapic_foreach_attr_set_grains(&attr, g, g);
  kaapic_foreach( first, last, &attr, 1, vtkParallelUtilitiesDoFor2, op );
  kaapic_end_parallel(KAAPIC_FLAG_DEFAULT);
  kaapic_foreach_attr_destroy(&attr);

  op->Finalize();
}

//--------------------------------------------------------------------------------
void vtkParallelUtilities::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf( os, indent );
}
