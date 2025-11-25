/*************************************************************************

   Program:    flexcalc
   File:       flexcalc.c
   
   Version:    V1.0
   Date:       24.11.25
   Function:   Calculate a flexibility score from an MD trajectory
   
   Copyright:  (c) Prof. Andrew C. R. Martin, abYinformatics, 2025
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

   The code proceeds as follows:
   1. Read the file through to obtain the number of frames
   2. Read a second time to calculate the average position for each atom
   3. Read a third time to find the frame closest to the average
   4. Read a fourth time to calculate the RMSD of each frame from the
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
#define MAXBUFF 512
#define MAXFNM  1024
#define PROGNAME "flexcalc"
#define MSG_ATOMMISMATCH "Number of coordinates in frame doesn't match \
first frame.\n  Frame Header: "
#define MSG_NOMEM "No memory"
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
FRAME *ReadFrame(FILE *in, char *header);
REAL  RMSFrame(FRAME *frame1, FRAME *frame2);
void  Usage(void);
FRAME *CopyFrame(FRAME *frame);
BOOL  AddFrame(FRAME *meanFrame, FRAME *frame, ULONG frameCount);
void  Die(char *msg, char *submsg);
void  Msg(char *msg, char *submsg);
void  PrintFrame(char *header, FRAME *frame);

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
            Die("No frames in trajectory", "");

         if((meanFrame    = CalculateMeanCoords(in, frameCount))==NULL)
            Die("Unable to calculate mean coordinates", "");

         if((closestFrame = FindClosestToMean(in, meanFrame,
                                              header))==NULL)
            Die("Couldn't find closest frame", header);

         if((meanRMSD     = CalculateMeanRMSD(in, closestFrame,
                                              frameCount)) < 0.0)
            Die("Unable to calculate mean RMSD", "");
         
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
   Msg(msg, submsg);
   exit(1);
}

/***********************************************************************/
/*>void Msg(char *msg, char *submsg)
   ---------------------------------
*//**
   \param[in]  *msg            Main messsage
   \param[in]  *submsg         Submessage

   Prints an error message

-  25.11.25 Original   By: ACRM
*/
void Msg(char *msg, char *submsg)
{
   fprintf(stderr, "%s error: %s%s\n", PROGNAME, msg, submsg);
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
   char  header[MAXBUFF];

   /* Reset the frame reading                                           */
   ReadFrame(NULL, NULL);

   /* Read frames, one at a time                                        */
   while((frame = ReadFrame(in, header))!=NULL)
   {
      REAL rmsd;
      /* Add the RMSD of this frame to the closest-to-mean frame and
         the free this frame
      */
      if((rmsd = RMSFrame(closestFrame, frame)) < 0.0)
      {
         Msg(MSG_ATOMMISMATCH, header);
         return(-1.0);
      }

#ifdef DEBUG
      PrintFrame(header, frame);
#endif
      
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
FRAME *FindClosestToMean(FILE *in, FRAME *meanFrame, char *header)
{
   FRAME *frame = NULL,
         *closestFrame = NULL;
   BOOL  firstFrame = TRUE;
   REAL  lowestRMSD;
   char  thisHeader[MAXBUFF];

   /* Reset the frame reading                                           */
   ReadFrame(NULL, NULL);

   /* Read frames, one at a time                                        */
   while((frame = ReadFrame(in, thisHeader))!=NULL)
   {
      REAL rmsd;

      /* Calculate RMSD to the mean frame                               */
      if((rmsd = RMSFrame(meanFrame, frame)) < 0.0)
      {
         Msg(MSG_ATOMMISMATCH, header);
         return(NULL);
      }

      /* If it was the first frame, make a copy of it and assume it's
         the best
      */
      if(firstFrame)
      {
         if((closestFrame = CopyFrame(frame))==NULL)
         {
            Msg(MSG_NOMEM, header);
            return(NULL);
         }
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
            if((closestFrame = CopyFrame(frame))==NULL)
               return(NULL);
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
   \param[in]  frameCount      the number of frames
   \return                     a pretend frame containing coordinates
                               averaged across the real frames

   Generates a new frame containing coordinates averaged across the
   other frames.
                               
-  24.11.25 Original   By: ACRM
*/
FRAME *CalculateMeanCoords(FILE *in, ULONG frameCount)
{
   FRAME *frame     = NULL,
         *meanFrame = NULL,
         *f;
   char  header[MAXBUFF];

   /* Reset the frame reading                                           */
   ReadFrame(NULL, NULL);

   /* Initialize the meanFrame - just read in the first frame then
      reset to zero coordinates, then rewind the file so we will
      re-read it.
   */
   meanFrame = ReadFrame(in, header);
   for(f=meanFrame; f!=NULL; NEXT(f))
   {
      f->x = 0.0;
      f->y = 0.0;
      f->z = 0.0;
   }
   rewind(in);
   
   /* Reset the frame reading                                           */
   ReadFrame(NULL, NULL);
   
   /* Read frames, one at a time                                        */
   while((frame = ReadFrame(in, header))!=NULL)
   {
      BOOL ok;
      ok = AddFrame(meanFrame, frame, frameCount);
      FREELIST(frame, FRAME);
      if(!ok)
      {
         Msg(MSG_ATOMMISMATCH, header);
         FREELIST(meanFrame, FRAME);
         meanFrame = NULL;
         break;
      }
   }
   rewind(in);

#ifdef DEBUG
   PrintFrame("average", meanFrame);
#endif
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
/*>FRAME *ReadFrame(FILE *in, char *header)
   ----------------------------------------
*//**
   \param[in]  *in             file pointer to trajectory
   \param[out] *header         the frame header
   \return                     the next frame from the file

   Reads a frame allocating a FRAME linked list

-  24.11.25 Original   By: ACRM
*/

FRAME *ReadFrame(FILE *in, char *header)
{
   static char buffer[MAXBUFF];
   static BOOL firstEntry = TRUE;
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
            break;
         else
            firstEntry = FALSE;
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

   if((f != NULL) || (g != NULL))
   {
      return(FALSE);
   }
   
   return(TRUE);
}


/***********************************************************************/
void PrintFrame(char *header, FRAME *frame)
{
   FRAME *f;
   
   printf("%s\n", header);
   for(f=frame; f!=NULL; NEXT(f))
   {
      printf("%.3f %.3f %.3f\n",f->x, f->y, f->z);
   }
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
   printf("\nflexcalc V1.0 (c) Andrew C.R. Martin, abYinformatics\n");

   printf("\nUsage: flexcalc trajectoryfile\n");

   printf("\nTakes a simple trajectory file in the format:\n");
   printf("      >frame header\n");
   printf("      x y z\n");
   printf("      x y z\n");
   printf("      ...\n");
   printf("      >frame header\n");
   printf("      x y z\n");
   printf("      x y z\n");
   printf("      ...\n");
   printf("      (etc)\n");

   printf("\nand calculates a flexibility score. This is done by \
calculating the\n");
   printf("mean coordinate postions across the frames, finding the \
frame closest\n");
   printf("to the mean positions and then calculating the RMSD of each \
frame to\n");
   printf("that closest-to-mean frame. These RMSD values are then \
averaged.\n");

   printf("\nThe code makes multiple passes of the file in order to \
minimize\n");
   printf("memory usage for large trajectories.\n\n");
}

