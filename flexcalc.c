/*************************************************************************

   Program:    flexcalc
   File:       flexcalc.c
   
   Version:    V1.0
   Date:       24.11.25
   Function:   Calculate a flexibility score from an MD trajectory
   
   Copyright:  (c) Prof. Andrew C. R. Martin 2025
   Author:     Prof. Andrew C. R. Martin
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   Licensed under the GPL V3.0. See the LICENCE file.

**************************************************************************

   Description:
   ============
   Calculates a flexibility score for a trajectory supplied in a simple
   format:
      >header
      x y z
      x y z
      ...
      >header
      x y z
      x y z
      ...
      (etc)

   The program minimizes memory usage by making multiple passes through
   the file.

   The code procedes as follows:
   1. Read the file through to obtain the number of frames
   2. Read a second time to calculate the average position for each atom
   3. Read a third time to find the frame closest to the average
   4. Read a fouth time to calculate the RMSD of each frame from the
      frame closest to the average
   5. Calculate and display the average of the RMSDs

**************************************************************************

   Usage:
   ======
   flexcalc in.pdb

**************************************************************************

   Revision History:
   =================
   V1.0   24.11.25 Original

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"

/***********************************************************************/
/* Defines and macros
 */
#define MAXBUFF 256
#define MAXFNM  1024
typedef struct _frame
{
   REAL x, y, z;
   struct _frame *next;
}  FRAME;


/***********************************************************************/
/* Prototypes
 */
BOOL  ParseCmdLine(int argc, char **argv, char *inFile);
FRAME *CalculateMeanCoords(FILE *in, ULONG frameCount);
FRAME *FindClosestToMean(FILE *in, FRAME *meanFrame, char *header);
REAL  CalculateMeanRMSD(FILE *in, FRAME *closestFrame, ULONG frameCount);
ULONG CountFrames(FILE *fp);
FRAME *ReadFrame(FILE *in);
REAL  RMSFrame(FRAME *frame1, FRAME *frame2);
void  Usage(void);
FRAME *CopyFrame(FRAME *frame);
BOOL  AddFrame(FRAME *meanFrame, FRAME *frame, ULONG frameCount);
void Die(char *msg, char *submsg);

/***********************************************************************/
/*>main(int argc, char **argv)
   ---------------------------
*//**
   Main program

-  24.11.25 Original   By: ACRM
*/
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
         if((frameCount   = CountFrames(in)) < 1)
            Die("flexcalc Error: no frames in trajectory","");

         if((meanFrame    = CalculateMeanCoords(in, frameCount))==NULL)
            Die("flexcalc Error: unable to calculate mean coordinates",
                "");

         if((closestFrame =
             FindClosestToMean(in, meanFrame, header))==NULL)
            Die("flexcalc Error: couldn't find closest frame","");

         if((meanRMSD     = CalculateMeanRMSD(in, closestFrame,
                                              frameCount)) < 0.0)
            Die("flexcalc Error: unable to calculate mean RMSD",
                "\n   Number of coordinates in frame doesn't match \
first frame");
         
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


/***********************************************************************/
/*>void Die(char *msg, char *submsg)
   ---------------------------------
*//**
   \param[in]  *msg            Main messsage
   \param[in]  *submsg         Submessage

   Prints an error message and exits

-  24.11.25 Original   By: ACRM
*/
void Die(char *msg, char *submsg)
{
   fprintf(stderr, "%s%s\n", msg, submsg);
   exit(1);
}


/***********************************************************************/
/*>REAL CalculateMeanRMSD(FILE *in, FRAME *closestFrame,
                          ULONG frameCount)
   -----------------------------------------------------
*//**
   \param[in]  *in             file pointer to trajectory
   \param[in]  *closestFrame   The frame closest to the mean coordinates
   \param[in]  frameCount      The number of frames
   \return                     The mean RMSD

   Reads through the frames and calculate an RMSD for each to the frame
   closest to the mean coordinates. Returns the average of these RMSD
   values.

-  24.11.25 Original   By: ACRM
*/
REAL CalculateMeanRMSD(FILE *in, FRAME *closestFrame, ULONG frameCount)
{
   REAL meanRMSD = 0.0;
   FRAME *frame = NULL;

   /* Reset the frame reading                                           */
   ReadFrame(NULL);

   /* Read frames, one at a time                                        */
   while((frame = ReadFrame(in))!=NULL)
   {
      REAL rmsd;
      /* Add the RMSD of this frame to the closest-to-mean frame and
         the free this frame
      */
      if((rmsd = RMSFrame(closestFrame, frame)) < 0.0)
         return(-1.0);
      
      meanRMSD += rmsd;
      FREELIST(frame, FRAME);
   }

   /* Divide by number of frames                                        */
   meanRMSD /= frameCount;
   rewind(in);
   return(meanRMSD);
}


/***********************************************************************/
/*>FRAME *FindClosestToMean(FILE *in, FRAME *meanFrame, char *header)
   ------------------------------------------------------------------
   
*//**
   \param[in]  *in             file pointer to trajectory
   \param[in]  *meanFrame      A pretend trajectory frame containing
                               the averaged coordinates
   \param[out] *header         The header for the frame closest to the
                               mean coordinates
   \return                     The frame closest to the mean coordinates

   Finds the frame closest to the averaged (pretend) frame.
   
-  24.11.25 Original   By: ACRM
*/
/* TODO - use header */

FRAME *FindClosestToMean(FILE *in, FRAME *meanFrame, char *header)
{
   FRAME *frame = NULL,
         *closestFrame = NULL;
   BOOL  firstFrame = TRUE;
   REAL  lowestRMSD;

   /* Reset the frame reading                                           */
   ReadFrame(NULL);

   /* Read frames, one at a time                                        */
   while((frame = ReadFrame(in))!=NULL)
   {
      REAL rmsd;

      /* Calculate RMSD to the mean frame                               */
      rmsd = RMSFrame(meanFrame, frame);

      /* If it was the first frame, make a copy of it and assume it's
         the best
      */
      if(firstFrame)
      {
         closestFrame = CopyFrame(frame);
         lowestRMSD = rmsd;
         firstFrame = FALSE;
      }
      else
      {
         /* Not the first, so if better, free up the copy of the current
            best and copy this to be the best
         */
         if(rmsd < lowestRMSD)
         {
            lowestRMSD = rmsd;
            FREELIST(closestFrame, FRAME);
            closestFrame = CopyFrame(frame);
         }
      }
      /* Free the current frame                                         */
      FREELIST(frame, FRAME);
   }
   rewind(in);
   return(closestFrame);
}


/***********************************************************************/
/*>FRAME *CalculateMeanCoords(FILE *in, ULONG frameCount)
   ------------------------------------------------------
*//**
   \param[in]  *in             file pointer to trajectory
   \param[in]  frameCount      the nnumber of frames
   \return                     a pretend frame containing coordinates
                               averaged across the real frames

   Generates a new frame containing coordinates averaged across the
   other frames.
                               
-  24.11.25 Original   By: ACRM
*/
FRAME *CalculateMeanCoords(FILE *in, ULONG frameCount)
{
   FRAME *frame = NULL,
         *meanFrame = NULL,
         *f;

   /* Reset the frame reading                                           */
   ReadFrame(NULL);

   /* Initialize the meanFrame - just read in the first frame then
      reset to zero coordinates, then rewind the file so we will
      re-read it.
   */
   meanFrame = ReadFrame(in);
   for(f=meanFrame; f!=NULL; NEXT(f))
   {
      f->x = 0.0;
      f->y = 0.0;
      f->z = 0.0;
   }
   rewind(in);
   
   /* Reset the frame reading                                           */
   ReadFrame(NULL);
   
   /* Read frames, one at a time                                        */
   while((frame = ReadFrame(in))!=NULL)
   {
      if(!AddFrame(meanFrame, frame, frameCount))
         Die("flexcalc Error: number of coordinates in frame doesn't \
match the first frame"," TODO: Frame header here!");
      FREELIST(frame, FRAME);
   }
   rewind(in);
   return(meanFrame);
}


/***********************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *inFile)
   ------------------------------------------------------
*//**
   \param[in]  argc            Argument count
   \param[in]  argv            Argument array
   \param[out] *inFile         Input filename from command line

   Parses the command line

-  24.11.25 Original   By: ACRM
*/
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

/***********************************************************************/
/*>FRAME *ReadFrame(FILE *in)
   --------------------------
*//**
   \param[in]  *in             file pointer to trajectory
   \return                     the next frame from the file

   Reads a frame allocating a FRAME linked list

-  24.11.25 Original   By: ACRM
*/

/* TODO output the header */
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


/***********************************************************************/
/*>ULONG CountFrames(FILE *fp)
   ---------------------------
*//**
   \param[in]  *in             file pointer to trajectory
   \return                     the number of frames in the file

   Read the number of frames from the file

-  24.11.25 Original   By: ACRM
*/
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


/***********************************************************************/
/*>REAL RMSFrame(FRAME *frame1, FRAME *frame2)
   -------------------------------------------
*//**
   \param[in]  *frame1         a frame linked list
   \param[in]  *frame2         a frame linked list

   Calculates the RMSD between two frames

-  24.11.25 Original   By: ACRM
*/
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

/* HERE */

/***********************************************************************/
/*>FRAME *CopyFrame(FRAME *frame)
   ------------------------------
*//**
   \param[in]  *frame          a frame linked list
   \return                     a newly allocated copy of the input
                               frame

   Copies a frame linked list                               

-  24.11.25 Original   By: ACRM
*/
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

/***********************************************************************/
/*>BOOL AddFrame(FRAME *meanFrame, FRAME *frame, ULONG frameCount)
   ---------------------------------------------------------------
*//**
   \param[out] *meanFrame      a pre-allocated pretend frame containing
                               averaged coordinates
   \param[in]  *frame          another frame linked list
   \param[in]  frameCount      the total number of frames
   \return                     Do the number of coordinates in the
                               frame match the meanFrame?

   Adds the coordinates for `frame` to the `meanFrame`.

   To increase the accuracy, the coordinates for each frame are
   divided by the number of frames before adding rather than adding
   everything first and dividing by the number of frames.

-  24.11.25 Original   By: ACRM
*/
BOOL AddFrame(FRAME *meanFrame, FRAME *frame, ULONG frameCount)
{
   FRAME *f, *g;
   f = meanFrame;
   g = frame;
   while((f!=NULL) && (g!=NULL))
   {
      f->x += (g->x / frameCount);
      f->y += (g->y / frameCount);
      f->z += (g->z / frameCount);

      NEXT(f);
      NEXT(g);
   }

   if((p != NULL) || (q != NULL))
   {
      return(FALSE);
   }
   
   return(TRUE);
}

/***********************************************************************/
/*>void Usage(void)
   ----------------
*//**
   Print usage message

-  24.11.25 Original   By: ACRM
*/
void Usage(void)
{
}

