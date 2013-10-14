/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkInitializableFunctor.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkInitializableFunctor.h"

//--------------------------------------------------------------------------------
vtkInitializableFunctor::vtkInitializableFunctor() : Initialized(0)
{
}

//--------------------------------------------------------------------------------
vtkInitializableFunctor::~vtkInitializableFunctor()
{
}

//--------------------------------------------------------------------------------
void vtkInitializableFunctor::Execute(vtkIdType begin, vtkIdType end) const
{
  unsigned char& inited = this->Initialized.Local();
  if (!inited)
    {
    this->Initialize();
    inited = 1;
    }
  this->operator()(begin, end);
}
