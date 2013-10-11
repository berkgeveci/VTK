/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkThreadedCompositeDataPipeline.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkThreadedCompositeDataPipeline.h"

#include "vtkAlgorithm.h"
#include "vtkAlgorithmOutput.h"
#include "vtkExecutive.h"
#include "vtkExecutiveCollection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkInformationExecutivePortKey.h"
#include "vtkInformationExecutivePortVectorKey.h"
#include "vtkInformationIntegerKey.h"
#include "vtkInformationObjectBaseKey.h"
#include "vtkInformationRequestKey.h"
#include "vtkInformationVector.h"

#include "vtkMultiThreader.h"
#include "vtkCompositeDataIterator.h"
#include "vtkCompositeDataSet.h"
#include "vtkTimerLog.h"
#include "vtkNew.h"
#include "vtkSmartPointer.h"
#include "vtkObjectFactory.h"
#include "vtkDebugLeaks.h"
#include "vtkImageData.h"

#include "vtkFunctor.h"
#include "vtkParallelUtilities.h"

#include <vector>
#include <assert.h>

int vtkThreadedCompositeDataPipeline::NumChunks  = 4;

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkThreadedCompositeDataPipeline);
vtkInformationKeyMacro(vtkThreadedCompositeDataPipeline, REQUEST_DIVIDE, Request);
vtkInformationKeyMacro(vtkThreadedCompositeDataPipeline, REQUEST_MERGE, Request);

//----------------------------------------------------------------------------
namespace
{
  inline bool Find(vtkInformation* info, vtkInformationVector** infoVec, int n)
  {
    for(int i=0; i<n; i++)
      {
      vtkInformationVector* infoVec_i = infoVec[i];
      for(int j=0; j<infoVec_i->GetNumberOfInformationObjects(); j++)
        {
        if(infoVec_i->GetInformationObject(j)==info)
          {
          return true;
          }
        }
      }
    return false;
  }

  static vtkInformationVector** Clone(vtkInformationVector** src, int n)
  {
    vtkInformationVector** dst = new vtkInformationVector*[n];
    for(int i=0; i<n; ++i)
      {
      dst[i] = vtkInformationVector::New();
      dst[i]->Copy(src[i],1);
      }
    return dst;
  }
  static void DeleteAll(vtkInformationVector** dst, int n)
  {
    for(int i=0; i<n; ++i)
      {
      dst[i]->Delete();
      }
    delete []dst;
  }
};

//----------------------------------------------------------------------------
class ProcessBlockData: public vtkObjectBase
{
public:
  vtkTypeMacro(ProcessBlockData, vtkObjectBase);
  vtkInformationVector** In;
  vtkInformationVector* Out;
  int InSize;

  static ProcessBlockData* New()
  {
    // This is required everytime we're implementing our own New() to avoid
    // "Deleting unknown object" warning from vtkDebugLeaks.
#ifdef VTK_DEBUG_LEAKS
    vtkDebugLeaks::ConstructClass("ProcessBlockData");
#endif
    return new ProcessBlockData();
  }

  void Construct(vtkInformationVector** inInfoVec,
                 int inInfoVecSize,
                 vtkInformationVector* outInfoVec)
  {
    this->InSize  = inInfoVecSize;
    this->In = Clone(inInfoVec, inInfoVecSize);
    this->Out = vtkInformationVector::New();
    this->Out->Copy(outInfoVec,1);
  }

  ~ProcessBlockData()
  {
    DeleteAll(this->In, this->InSize);
    this->Out->Delete();
  }

protected:
  ProcessBlockData():
    In(NULL),
    Out(NULL)
  {

  }
};
//----------------------------------------------------------------------------
class ProcessBlock : public vtkFunctor
{
public:
  ProcessBlock(vtkThreadedCompositeDataPipeline* exec,
               vtkInformationVector** inInfoVec,
               vtkInformationVector* outInfoVec,
               int compositePort,
               int connection,
               vtkInformation* request,
               const std::vector<vtkDataObject*>& inObjs,
               std::vector<vtkDataObject*>& outObjs
  )
    :Exec(exec),
     InInfoVec(inInfoVec),
     OutInfoVec(outInfoVec),
     CompositePort(compositePort),
     Connection(connection),
     Request(request),
     InObjs(inObjs)
  {
    int numInputPorts = this->Exec->GetNumberOfInputPorts();
    this->OutObjs = &outObjs[0];
    this->InfoPrototype = vtkSmartPointer<ProcessBlockData>::New();
    this->InfoPrototype->Construct(this->InInfoVec, numInputPorts, this->OutInfoVec);
  }

  ~ProcessBlock()
  {
  }

  void operator() (vtkIdType begin, vtkIdType end) const
  {
    bool IamFirstChunk=  begin==0;
    vtkInformationVector** inInfoVec;
    vtkInformationVector* outInfoVec;

    if(IamFirstChunk)
      {
      inInfoVec = this->InInfoVec;
      outInfoVec = this->OutInfoVec;
      }
    else
      {
      inInfoVec = Clone(this->InfoPrototype->In, this->InfoPrototype->InSize);
      outInfoVec = vtkInformationVector::New();
      outInfoVec->Copy(this->InfoPrototype->Out,1);
      }

    vtkSmartPointer<vtkInformation> request = vtkSmartPointer<vtkInformation>::New();
    request->Copy(this->Request,1);

    vtkInformation* inInfo = inInfoVec[this->CompositePort]->GetInformationObject(this->Connection);
    vtkInformation* outInfo = outInfoVec->GetInformationObject(0);

    for(vtkIdType i= begin; i<end; ++i)
      {
      vtkDataObject* outObj =
        this->Exec->ExecuteSimpleAlgorithmForBlock(&inInfoVec[0],
                                                   outInfoVec,
                                                   inInfo,
                                                   outInfo,
                                                   request,
                                                   this->InObjs[i]);
      this->OutObjs[i] = outObj;
      }

    if(!IamFirstChunk)
      {
      DeleteAll(inInfoVec, this->InfoPrototype->InSize);
      outInfoVec->Delete();
      }
  }

protected:
  vtkThreadedCompositeDataPipeline* Exec;
  vtkInformationVector** InInfoVec;
  vtkInformationVector* OutInfoVec;
  vtkSmartPointer<ProcessBlockData> InfoPrototype;
  int CompositePort;
  int Connection;
  vtkInformation* Request;
  const std::vector<vtkDataObject*>& InObjs;
  vtkDataObject** OutObjs;
};


//----------------------------------------------------------------------------
vtkThreadedCompositeDataPipeline::vtkThreadedCompositeDataPipeline()
{
}

//----------------------------------------------------------------------------
vtkThreadedCompositeDataPipeline::~vtkThreadedCompositeDataPipeline()
{
}

//-------------------------------------------------------------------------
void vtkThreadedCompositeDataPipeline::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//-------------------------------------------------------------------------
void vtkThreadedCompositeDataPipeline::ExecuteEach(vtkCompositeDataIterator* iter,
                                                   vtkInformationVector** inInfoVec,
                                                   vtkInformationVector* outInfoVec,
                                                   int compositePort,
                                                   int connection,
                                                   vtkInformation* request,
                                                   vtkCompositeDataSet* compositeOutput)
{
  vtkInformation* inInfo  =inInfoVec[compositePort]->GetInformationObject(connection);

  assert(Find(inInfo,inInfoVec, this->GetNumberOfInputPorts()));

  // from input data objects  itr -> (inObjs, indices)
  // inObjs are the non-null objects that we will loop over.
  // indices map the input objects to inObjs
  std::vector<vtkDataObject*> inObjs;
  std::vector<int> indices;
  for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem())
    {
    vtkDataObject* dobj = iter->GetCurrentDataObject();
    if (dobj)
      {
      inObjs.push_back(dobj);
      indices.push_back(inObjs.size()-1);
      }
    else
      {
      indices.push_back(-1);
      }
    }

  // instantiate outObjs, the output objects that will be created from inObjs
  std::vector<vtkDataObject*> outObjs;
  outObjs.resize(indices.size(),NULL);

  // create the parallel task processBlock
  ProcessBlock processBlock(this,
                            inInfoVec,
                            outInfoVec,
                            compositePort,
                            connection,
                            request,
                            inObjs,outObjs);

  vtkParallelUtilities::ForEach(0, inObjs.size(), &processBlock);

  size_t i =0;
  for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem(), i++)
    {
    int j = indices[i];
    if(j>=0)
      {
      vtkDataObject* outObj = outObjs[j];
      compositeOutput->SetDataSet(iter, outObj);
      if(outObj)
        {
        outObj->FastDelete();
        }
      }
    }
}
