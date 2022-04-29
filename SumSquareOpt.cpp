//---------------------------------------------------------------------------

#include "ImageTreatment.h"
#include "PrincipalAxis.h"
#include "SumSquareOpt.h"

//---------------------------------------------------------------------------

bool getSurfel(Image<float> *im,Image<float> *imMove,float *A,Image<bool> *mask,int type)
{
  int posX,posY;
  float x,y,deltaX,deltaY;

  int index=0;
  for(int j=0;j<im->height;j++)
    for(int i=0;i<im->width;i++)
    {
      (*imMove)[i]=0;
      if (mask==NULL || (*mask)[index]==true)
      {
        x=A[0]*(i-im->width/2)+A[1]*(j-im->height/2)+A[2]+im->width/2;
        y=A[3]*(i-im->width/2)+A[4]*(j-im->height/2)+A[5]+im->height/2;
        (*imMove)[index]=getValue(im,x,y,type);
      }
      index++;
    }

  return true;
}

//---------------------------------------------------------------------------

float *ssoRealign2D(int nbParams,int *U,Image<float> *imageRef,
                    Image<float> *imageToRegister)
{
  int nbIter=25;
  float a=1;         // translation factor (pixels)
  float b=M_PI/180;  // rotation factor (radians)
  float c=1.0/100;   // scale factor
  float d=1.0/100;   // cisaillement factor
  float paramsThreshold = 0.01;

  float *dQ = new float[6];
  Image<float> *d1 = new Image<float>(imageRef->width,imageRef->height,1,1);
  Image<float> *d2 = new Image<float>(imageRef->width,imageRef->height,1,1);
  Image<float> *X = new Image<float>(imageRef->width,imageRef->height,1,1);
  Image<float> *Y = new Image<float>(imageRef->width,imageRef->height,1,1);
  float *Q = new float[nbParams];
  float *q0 = new float[nbParams];
  float *matA = new float[nbParams*nbParams];
  float *matB = new float[nbParams*nbParams];
  float *dXdQ = new float[imageRef->width*imageRef->height*nbParams];
  float *P = new float[6];
  bool *change = new bool[nbParams];

  dQ[0]=a;  dQ[1]=a;  dQ[2]=b;  dQ[3]=c;  dQ[4]=c;  dQ[5]=d;

  for(int i=0;i<nbParams;i++)
  {
    Q[i]=0;
    q0[i]=0.0;
  }

  float *imRefData=imageRef->getData();
  float *imageToRegisterData=imageToRegister->getData();

  // Build a region excluding noise area on both reference image and the
  // image to register in order to optimize the computation time
  Image<bool> *mask = new Image<bool>(imageRef->width,imageRef->height,1,1);
  for(int i=0;i<imageRef->width*imageRef->height;i++)
    if (imRefData[i]>=Magnitude_Threshold_Noise ||
        imageToRegisterData[i]>Magnitude_Threshold_Noise)
      (*mask)[i]=true;
    else
      (*mask)[i]=false;

  // Analyse principal axis in a first estimation
  float posX1,posY1,angle1,dimensionX1,dimensionY1;
  analysePrincipalAxis(imageRef,&posX1,&posY1,&angle1,&dimensionX1,&dimensionY1,mask->getData());
  float posX2,posY2,angle2,dimensionX2,dimensionY2;
  analysePrincipalAxis(imageToRegister,&posX2,&posY2,&angle2,&dimensionX2,&dimensionY2,mask->getData());
  Q[0]=-(posX2-posX1);
  Q[1]=-(posY2-posY1);

  for(int i=0;i<6;i++) P[i]=0;

  for(int iBis=0;iBis<6;iBis++) P[iBis]=0;
  getSurfel(imageToRegister,X,transformationMatrix2D(P),mask,BILINEAR);
  for(int j=0;j<nbParams;j++)
  {
    for(int iBis=0;iBis<6;iBis++) P[iBis]=0;
    P[U[j]]=-dQ[U[j]];
    getSurfel(imageToRegister,d1,transformationMatrix2D(P),mask,BILINEAR);
    for(int iBis=0;iBis<6;iBis++) P[iBis]=0;
    P[U[j]]=dQ[U[j]];
    getSurfel(imageToRegister,d2,transformationMatrix2D(P),mask,BILINEAR);
    for(int k=0;k<imageRef->width*imageRef->height;k++)
      dXdQ[k*nbParams+j]=((*d2)[k]-(*d1)[k])/2;
  }

  for(int k=0;k<nbParams*nbParams;k++) matA[k]=0;
  int index=0,indexI=0,indexJ=0,indexK=0;
  for(int k=0;k<imageRef->width*imageRef->height;k++)
  {
    if ((*mask)[k])
    {
      indexJ=indexK;
      index=0;
      for(int jBis=0;jBis<nbParams;jBis++)
      {
        indexI=indexK;
        for(int iBis=0;iBis<nbParams;iBis++)
        {
          matA[index]+=dXdQ[indexI]*dXdQ[indexJ];
          index++;
          indexI++;
        }
        indexJ++;
      }
    }
    indexK+=nbParams;
  }

  for(int i=0;i<nbParams;i++) change[i]=true;
  bool find=false;
  int i=0;

  do
  {
    for(int iBis=0;iBis<6;iBis++) P[iBis]=0;
    for(int iBis=0;iBis<nbParams;iBis++) P[U[iBis]]=Q[iBis];

    getSurfel(imageRef,Y,transformationMatrix2D(P),mask,BILINEAR);

    for(int iBis=0;iBis<nbParams;iBis++) matB[iBis]=0;
    int index=0;
    for(int k=0;k<imageRef->width*imageRef->height;k++)
      for(int iBis=0;iBis<nbParams;iBis++)
      {
        matB[iBis]+=dXdQ[index]*((*Y)[k]-(*X)[k]);
        index++;
      }


    q0=solveSyst(matA,matB,nbParams);

    find=true;
    for(int iBis=0;iBis<nbParams;iBis++)
      if (fabs(q0[iBis])>=paramsThreshold)
      {
        change[iBis]=true;
        find=false;
      }
      else
      {
        q0[iBis]=0;
        change[iBis]=false;
      }

    for(int iBis=0;iBis<nbParams;iBis++)
      Q[iBis]-=q0[iBis]*dQ[U[iBis]];

    i++;
  } while(find==false && i<nbIter);

  float params[6];
  for(int i=0;i<6;i++) params[i]=0;
  for(int i=0;i<nbParams;i++) params[U[i]]=-Q[i];

  delete[] dQ;
  delete d1;
  delete d2;
  delete X;
  delete Y;
  delete[] Q;
  delete[] q0;
  delete[] matA;
  delete[] matB;
  delete[] dXdQ;
  delete[] P;
  delete[] change;
  delete mask;

  return transformationMatrix2D(params);
}

//---------------------------------------------------------------------------

float *ssoRealign2D(int mode,Image<float> *imageRef,
                    Image<float> *imageToRegister)
{
  int nbParams;

  switch(mode)
  {
    case 0:
      nbParams=3; // TranslationX+TranslationY+Rotation
    break;
    case 1:
      nbParams=5; // TranslationX+TranslationY+Rotation+ZoomX+ZoomY
    break;
    case 2:
      nbParams=6; // TranslationX+TranslationY+Rotation+ZoomX+ZoomY+Cisaillement
    break;
  }
  int *U = new int[nbParams];
  for(int i=0;i<nbParams;i++) U[i]=i;
  float *transformation = ssoRealign2D(nbParams,U,imageRef,imageToRegister);
  delete[] U;

  return transformation;
}

//---------------------------------------------------------------------------

 
