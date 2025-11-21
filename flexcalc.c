#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"

#define MAXBUFF 256
#define MAXFNM  1024
typedef struct _frame
{
   REAL x, y, z;
   struct _frame *next;
}  FRAME;


BOOL ParseCmdLine(int argc, char **argv, char *inFile);
FRAME *CalculateMeanCoords(FILE *in, ULONG frameCount);
FRAME *FindClosestToMean(FILE *in, FRAME *meanFrame, char *header);
REAL CalculateMeanRMSD(FILE *in, FRAME *closestFrame, ULONG frameCount);
ULONG CountFrames(FILE *fp);
FRAME *ReadFrame(FILE *in);
REAL RMSFrame(FRAME *frame1, FRAME *frame2);
void Usage(void);
FRAME *CopyFrame(FRAME *frame);
void AddFrame(FRAME *meanFrame, FRAME *frame, ULONG frameCount);

int main(int argc, char **argv)
{
   FILE *in;
   char inFile[MAXFNM];
   inFile[0] = '\0';
   
   if(ParseCmdLine(argc, argv, inFile))
   {
      if((in=fopen(inFile, "r"))!=NULL)
      {
         char  header[MAXBUFF];
         FRAME *meanFrame = NULL,
               *closestFrame = NULL;
         REAL  meanRMSD;
         ULONG frameCount = 0;
         frameCount = CountFrames(in);
         meanFrame    = CalculateMeanCoords(in, frameCount);
         closestFrame = FindClosestToMean(in, meanFrame, header);
         meanRMSD     = CalculateMeanRMSD(in, closestFrame, frameCount);
         fclose(in);

         printf("%.4f\n", meanRMSD);
      }
   }
   else
   {
      Usage();
   }

   return(0);
}

REAL CalculateMeanRMSD(FILE *in, FRAME *closestFrame, ULONG frameCount)
{
   REAL meanRMSD = 0.0;
   FRAME *frame = NULL;

   ReadFrame(NULL);
   
   while((frame = ReadFrame(in))!=NULL)
   {
      meanRMSD += RMSFrame(closestFrame, frame);
      FREELIST(frame, FRAME);
   }
   meanRMSD /= frameCount;
   rewind(in);
   return(meanRMSD);
}

FRAME *FindClosestToMean(FILE *in, FRAME *meanFrame, char *header)
{
   FRAME *frame = NULL,
         *closestFrame = NULL;
   BOOL  firstFrame = TRUE;
   REAL  lowestRMSD;

   ReadFrame(NULL);

   while((frame = ReadFrame(in))!=NULL)
   {
      REAL rmsd;
      rmsd = RMSFrame(meanFrame, frame);
      if(firstFrame)
      {
         closestFrame = CopyFrame(frame);
         lowestRMSD = rmsd;
         firstFrame = FALSE;
      }
      else
      {
         if(rmsd < lowestRMSD)
         {
            lowestRMSD = rmsd;
            FREELIST(closestFrame, FRAME);
            closestFrame = CopyFrame(frame);
         }
      }
      FREELIST(frame, FRAME);
   }
   rewind(in);
   return(closestFrame);
}

FRAME *CalculateMeanCoords(FILE *in, ULONG frameCount)
{
   FRAME *frame = NULL,
         *meanFrame = NULL,
         *f;

   ReadFrame(NULL);

   /* Initialize the meanFram */
   meanFrame = ReadFrame(in);
   for(f=meanFrame; f!=NULL; NEXT(f))
   {
      f->x = 0.0;
      f->y = 0.0;
      f->z = 0.0;
   }
   rewind(in);
   ReadFrame(NULL);
   
   while((frame = ReadFrame(in))!=NULL)
   {
      AddFrame(meanFrame, frame, frameCount);
      FREELIST(frame, FRAME);
   }
   rewind(in);
   return(meanFrame);
}

BOOL ParseCmdLine(int argc, char **argv, char *inFile)
{
   argc--; argv++;
   if(argc > 0)
   {
      strcpy(inFile, argv[0]);
      argc--; argv++;
   }
   return(TRUE);
}

FRAME *ReadFrame(FILE *in)
{
   static char buffer[MAXBUFF];
   static BOOL firstEntry = TRUE;
   char header[MAXBUFF];
   FRAME *frame = NULL,
         *f = NULL;

   if(in==NULL)
   {
      firstEntry = TRUE;
      return(NULL);
   }
   
   /* Copy the existing buffer - which should be the next header        */
   if(firstEntry)
   {
      buffer[0] = '\0';
   }
   else
   {
      strncpy(header, buffer, MAXBUFF-1);
   }

   while(fgets(buffer, MAXBUFF-1, in))
   {
      TERMINATE(buffer);
      
      if(buffer[0] == '>')  /* A header                                 */
      {
         if(!firstEntry)
         {
            strncpy(header, buffer, MAXBUFF-1);
            break;
         }
         else
         {
            firstEntry = FALSE;
         }
      }
      else
      {
         if(frame == NULL)
         {
            INIT(frame, FRAME);
            f = frame;
         }
         else
         {
            ALLOCNEXT(f, FRAME);
         }
         if(f==NULL)
         {
            FREELIST(frame, FRAME);
            return(NULL);
         }

         sscanf(buffer, "%lf %lf %lf", &(f->x), &(f->y), &(f->z));
      }
   }
   
   return(frame);
}


/*****************************************************************/


ULONG CountFrames(FILE *fp)
{
   ULONG frameCount = 0;
   char buffer[MAXBUFF];

   while(fgets(buffer, MAXBUFF-1, fp))
   {
      if(buffer[0] == '>')
         frameCount++;
   }
   rewind(fp);
   return(frameCount);
}


REAL RMSFrame(FRAME *frame1, FRAME *frame2)
{
   REAL  rmsd;
   int nCoor;
   FRAME *p, *q;
   p = frame1;
   q = frame2;

   rmsd  = 0.0;
   nCoor = 0;

   while((p!=NULL) && (q!=NULL))
   {
      rmsd += (p->x - q->x) * (p->x - q->x) +
              (p->y - q->y) * (p->y - q->y) +
              (p->z - q->z) * (p->z - q->z);
      nCoor++;
      NEXT(p);
      NEXT(q);
   }

   if((p != NULL) || (q != NULL))
   {
      return(-1.0);
   }
   
   return(sqrt(rmsd/nCoor));
}

FRAME *CopyFrame(FRAME *frame)
{
   FRAME *copy = NULL,
         *c, *f;

   for(f=frame; f!=NULL; NEXT(f))
   {
      if(copy==NULL)
      {
         INIT(copy, FRAME);
         c = copy;
      }
      else
      {
         ALLOCNEXT(c, FRAME);
      }
      if(c==NULL)
      {
         FREELIST(copy, FRAME);
         return(NULL);
      }

      c->x = f->x;
      c->y = f->y;
      c->z = f->z;
   }

   return(copy);
}

void AddFrame(FRAME *meanFrame, FRAME *frame, ULONG frameCount)
{
   meanFrame->x += (frame->x / frameCount);
   meanFrame->y += (frame->y / frameCount);
   meanFrame->z += (frame->z / frameCount);
}

void Usage(void)
{
}

