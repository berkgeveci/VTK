 /*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkThreadLocalObject.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkThreadLocalObject -
// .SECTION Description

#ifndef __vtkThreadLocalObject_h
#define __vtkThreadLocalObject_h

#include <vtkThreadLocal.h>

template <typename T>
class vtkThreadLocalObject
{
  typedef vtkThreadLocal<T*> TLS;
  typedef typename vtkThreadLocal<T*>::iterator TLSIter;
public:
  // Description:
  // Default constructor.
  vtkThreadLocalObject()
    {
    }

  virtual ~vtkThreadLocalObject()
    {
      iterator iter = this->begin();
      while (iter != this->end())
        {
        if (*iter)
          {
          (*iter)->Delete();
          }
        ++iter;
        }
    }

  // Description:
  T*& Local()
    {
      T*& vtkobject = this->Internal.Local();
      if (!vtkobject)
        {
        vtkobject = T::New();
        }
      return vtkobject;
    }

  // Description:
  // Subset of the standard iterator API.
  // The most common design patter is to use iterators in a sequential
  // code block and to use only the thread local objects in parallel
  // code blocks.
  class iterator
  {
  public:
    iterator& operator++()
      {
        ++this->Iter;
        return *this;
      }

    bool operator!=(const iterator& other)
      {
        return this->Iter != other.Iter;
      }

    T*& operator*()
      {
        return *this->Iter;
      }

  private:
    TLSIter Iter;

    friend class vtkThreadLocalObject<T>;
  };

  iterator begin()
    {
      iterator iter;
      iter.Iter = this->Internal.begin();
      return iter;
    };

  iterator end()
    {
      iterator iter;
      iter.Iter = this->Internal.end();
      return iter;
    }

private:
  TLS Internal;
};

#endif
// VTK-HeaderTest-Exclude: vtkThreadLocalObject.h
