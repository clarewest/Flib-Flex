#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>     

int main(int argc, char *argv[])
{
	int length,start1,match_score,seq_score,start2,ss_score,total,old=-1;
	char Header[82],Chain,type,ALN_Seq[30],c;
	char Fasta_Seq[1000];
	char AUX[82];
	float torsion,rmsd,score,rmsd_threshold;
	double resolution;
	FILE *library;
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
	library = fopen(argv[2],"r");
	if(library==NULL) { fprintf(stderr,"Library file: %s not found!\n",argv[2]); 	return 0; 	}

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

	for(total=0;fscanf(library,"%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%lf\t%d\t%f\t%f\n",Header,&Chain,&start1,&length,ALN_Seq,&type,&match_score,&seq_score,&length,&start2,&resolution,&ss_score,&torsion,&rmsd)!=EOF;total++)
	{
		if(old != start2) /* Fragments for a new position... */
		{
			total=0;
			old = start2;
		}
                
                if(rmsd <= rmsd_threshold)
                {
                        fprintf(out,"%s\t%c\t%3d\t%3d\t%s\t%c\t%3d\t%3d\t%3d\t%3d\t%.2lf\t%d\t%.2f\t%.2f\n",Header,Chain,start1,length,ALN_Seq,type,match_score,seq_score,length,start2,resolution,ss_score,torsion,rmsd);
                }else if(total < 50)
                {
		        fprintf(out,"%s\t%c\t%3d\t%3d\t%s\t%c\t%3d\t%3d\t%3d\t%3d\t%.2lf\t%d\t%.2f\t%.2f\n",Header,Chain,start1,length,ALN_Seq,type,match_score,seq_score,length,start2,resolution,ss_score,torsion,rmsd);
                }
	}
	fclose(library);
	fclose(input_fasta);
	fclose(out);
    
	return 0;
}
