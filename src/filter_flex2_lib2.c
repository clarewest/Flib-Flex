#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>   
#include<assert.h>

int main(int argc, char *argv[])
{
  int length,start1,match_score,seq_score,start2,ss_score,total,old=-1;
  int mode = -1;
  char Header[82],Chain,type,ALN_Seq[30],c;
  char Fasta_Seq[1000];
  char AUX[82];
  float torsion,rmsd,score,rmsd_threshold;
  double resolution;
  FILE *library_a, *library_b, *library_c;
  FILE *input_fasta;
  FILE *out;

  if(argc != 4)
  {
    printf("USAGE: %s PDB_ID FLIB_FILE.lib3000_rmsd rmsd_threshold\n",argv[0]);
    return 0;
  }

  /*** SET RMSD THRESHOLD ***/
  rmsd_threshold = atof(argv[3]);

  /*** FILE HANDLING: ***/
  sprintf(AUX, "%s_a", argv[2]);
  library_a = fopen(AUX,"r");
  if(library_a==NULL) { fprintf(stderr,"Library file: %s not found!\n",AUX); 	return 0; 	}
  sprintf(AUX, "%s_b", argv[2]);
  library_b = fopen(AUX,"r");
  if(library_b==NULL) { fprintf(stderr,"Library file: %s not found!\n",AUX); 	return 0; 	}
  sprintf(AUX, "%s_c", argv[2]);
  library_c = fopen(AUX,"r");
  if(library_c==NULL) { fprintf(stderr,"Library file: %s not found!\n",AUX); 	return 0; 	}

  strcpy(AUX,argv[1]);
  input_fasta = fopen(strcat(AUX,".fasta.txt"),"r");
  if(input_fasta==NULL) { fprintf(stderr,"Fasta file: %s not found!\n",AUX); 	return 0;	}

  /* READ INPUT FASTA SEQUENCE */
  /* Remove Header from FASTA file */
  for(c = fgetc(input_fasta); c!='\n' ; c=fgetc(input_fasta));
  fscanf(input_fasta,"%s",Fasta_Seq);

  while(fscanf(input_fasta,"%s",AUX)!=EOF) 	
    strcat(Fasta_Seq,AUX);



  strcpy(AUX,argv[1]);
  sprintf(AUX, "%s.lib3000_flex%.1f", AUX, rmsd_threshold);
  out = fopen(AUX,"w");

  for(total=0;fscanf(library_a,"%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%lf\t%d\t%f\t%f\n",Header,&Chain,&start1,&length,ALN_Seq,&type,&match_score,&seq_score,&length,&start2,&resolution,&ss_score,&torsion,&rmsd)!=EOF;total++)
  {

    assert(old <= start2);

    if(old != start2) /* Fragments for a new position... */
    {
      total=0;
      old = start2;
      
      //decide which library to output depending on rmsd at the position
      if(rmsd > -1)
      {
        assert(fscanf(library_c,"%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%lf\t%d\t%f\t%f\n",Header,&Chain,&start1,&length,ALN_Seq,&type,&match_score,&seq_score,&length,&start2,&resolution,&ss_score,&torsion,&rmsd));
        if(rmsd > -1)
        {
          mode = 1;
        }else
        {
          mode = 2;
        }
      }else
      {
        assert(fscanf(library_c,"%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%lf\t%d\t%f\t%f\n",Header,&Chain,&start1,&length,ALN_Seq,&type,&match_score,&seq_score,&length,&start2,&resolution,&ss_score,&torsion,&rmsd));
        if(rmsd == -1)
        {
          mode = 3;
        }else
        {
          fprintf(stderr, "rmsd not sorted correctly in files.\n");
          return 1;
        }
      }

    }else
    {
      assert(fscanf(library_c,"%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%lf\t%d\t%f\t%f\n",Header,&Chain,&start1,&length,ALN_Seq,&type,&match_score,&seq_score,&length,&start2,&resolution,&ss_score,&torsion,&rmsd));
    }
    //at this point I have read library_c
 
    if(mode == 1)
    {
      // if rmsd are all > -1, then output those with rmsd <= rmsd_threshold and fill up if necessary
      // i.e., output c

      if(rmsd <= rmsd_threshold)
      {
        fprintf(out,"%s\t%c\t%3d\t%3d\t%s\t%c\t%3d\t%3d\t%3d\t%3d\t%.2lf\t%d\t%.2f\t%.2f\n",Header,Chain,start1,length,ALN_Seq,type,match_score,seq_score,length,start2,resolution,ss_score,torsion,rmsd);
      }else if(total < 20)
      {
        fprintf(out,"%s\t%c\t%3d\t%3d\t%s\t%c\t%3d\t%3d\t%3d\t%3d\t%.2lf\t%d\t%.2f\t%.2f\n",Header,Chain,start1,length,ALN_Seq,type,match_score,seq_score,length,start2,resolution,ss_score,torsion,rmsd);
      }

      assert(fscanf(library_b,"%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%lf\t%d\t%f\t%f\n",Header,&Chain,&start1,&length,ALN_Seq,&type,&match_score,&seq_score,&length,&start2,&resolution,&ss_score,&torsion,&rmsd));
    }else if(mode == 2 || mode == 3)
    {

      //if rmsd are not all > -1, then output all with rmsc <= rmsd_threshold, sorted by torsion angle.
      //i.e., output b
      assert(fscanf(library_b,"%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%lf\t%d\t%f\t%f\n",Header,&Chain,&start1,&length,ALN_Seq,&type,&match_score,&seq_score,&length,&start2,&resolution,&ss_score,&torsion,&rmsd));


      if(rmsd <= rmsd_threshold)
      {
        fprintf(out,"%s\t%c\t%3d\t%3d\t%s\t%c\t%3d\t%3d\t%3d\t%3d\t%.2lf\t%d\t%.2f\t%.2f\n",Header,Chain,start1,length,ALN_Seq,type,match_score,seq_score,length,start2,resolution,ss_score,torsion,rmsd);
      }else if(total < 20)
      {
        fprintf(out,"%s\t%c\t%3d\t%3d\t%s\t%c\t%3d\t%3d\t%3d\t%3d\t%.2lf\t%d\t%.2f\t%.2f\n",Header,Chain,start1,length,ALN_Seq,type,match_score,seq_score,length,start2,resolution,ss_score,torsion,rmsd);
      }


    }else
    {
      fprintf(stderr, "Error.\n");
    }

  }
  fclose(library_a);
  fclose(library_b);
  fclose(library_c);
  fclose(input_fasta);
  fclose(out);

  return 0;
}
