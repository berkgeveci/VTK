/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTimeStamp.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkTimeStamp.h"

#include "vtkObjectFactory.h"
#include "vtkWindows.h"

#include "vtkAtomicInt.h"

static bool vtkTimeStampTimeFrozen = false;

//-------------------------------------------------------------------------
vtkTimeStamp* vtkTimeStamp::New()
{
  // If the factory was unable to create the object, then create it here.
  return new vtkTimeStamp;
}

//-------------------------------------------------------------------------
void vtkTimeStamp::FreezeTime()
{
  vtkTimeStampTimeFrozen = true;
}

//-------------------------------------------------------------------------
void vtkTimeStamp::UnfreezeTime()
{
  vtkTimeStampTimeFrozen = false;
}

//-------------------------------------------------------------------------
void vtkTimeStamp::Modified()
{
#if VTK_SIZEOF_VOID_P == 8
  static vtkAtomicInt<vtkTypeInt64> GlobalTimeStamp(0);
#else
  static vtkAtomicInt<vtkTypeInt32> GlobalTimeStamp(0);
#endif

  if (!vtkTimeStampTimeFrozen)
    {
    this->ModifiedTime = (unsigned long)++GlobalTimeStamp;
    }
}
