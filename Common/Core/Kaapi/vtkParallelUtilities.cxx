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
  if (!vtkParallelUtilitiesInitialized)
    {
    vtkParallelUtilitiesCS.Lock();
    if (!vtkParallelUtilitiesInitialized)
      {
      static vtkParallelUtilitiesInit aInit;
      vtkParallelUtilitiesInitialized = true;
      }
    vtkParallelUtilitiesCS.Unlock();
    }
}

namespace
{
inline void vtkParallelUtilitiesDoFor(int32_t b, int32_t e, int32_t, const vtkFunctor* o )
  {
    (*o)(b, e);
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

  kaapic_begin_parallel(KAAPIC_FLAG_DEFAULT);
  kaapic_foreach_attr_t attr;
  kaapic_foreach_attr_init(&attr);
  if (grain > 0)
    {
    kaapic_foreach_attr_set_grains(&attr, grain, grain);
    }
  kaapic_foreach( first, last, &attr, 1, vtkParallelUtilitiesDoFor, op );
  kaapic_end_parallel(KAAPIC_FLAG_DEFAULT);
  kaapic_foreach_attr_destroy(&attr);
}

//--------------------------------------------------------------------------------
void vtkParallelUtilities::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf( os, indent );
}
