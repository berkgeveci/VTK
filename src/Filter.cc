/*=========================================================================

  Program:   Visualization Library
  Module:    Filter.cc
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

This file is part of the Visualization Library. No part of this file or its 
contents may be copied, reproduced or altered in any way without the express
written consent of the authors.

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen 1993, 1994 

=========================================================================*/
#include "Filter.hh"

void vlFilter::SetStartMethod(void (*f)())
{
  if ( f != this->StartMethod )
    {
    this->StartMethod = f;
    this->Modified();
    }
}

void vlFilter::SetEndMethod(void (*f)())
{
  if ( f != this->EndMethod )
    {
    this->EndMethod = f;
    this->Modified();
    }
}

void vlFilter::PrintSelf(ostream& os, vlIndent indent)
{
  if (this->ShouldIPrint(vlFilter::GetClassName()))
    {
    vlObject::PrintSelf(os,indent);

    if ( this->StartMethod )
      {
      os << indent << "StartMethod: (" << this->StartMethod << ")\n";
      }
    else
      {
      os << indent << "StartMethod: (none)\n";
      }

    if ( this->EndMethod )
      {
      os << indent << "EndMethod: (" << this->EndMethod << ")\n";
      }
    else
      {
      os << indent << "EndMethod: (none)\n";
      }

    os << indent << "Execute time: " <<this->ExecuteTime.GetMtime() << "\n";
   }
}
