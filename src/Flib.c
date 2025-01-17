#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<ctype.h>
#include"slink.h"
#include<getopt.h>

#define RANDOM 1
#define EXHAUSTIVE 1
#define MaxFrag 10000
#define TRUE 1
#define FALSE 0
#define CUTOFF 1.5
#define MATCH 2
#define MISMATCH -2
#define GAP -100
#define VERBOSE 1
#define MAXLEN 10000
#define MAX_FIT_COORDS 200

#define max(A,B) (((A)>(B))?(A):(B))
#define min(A,B) (((A)<(B))?(A):(B))
#define score(A,B) (((A)==(B))?(MATCH):(MISMATCH))

/* Collect the coordinates from the pdb file starting from the correct residue. */
int coordcol( char bum[FILEL],int u, int length, int start_res,char ChainQ,char* Seq);
int coordcol2( FRAGMENTS * Piece,char PDB[82],char Chain, char* Seq,int *start );
int MAX_FRAG=3000;
int COEVO=0;
int COEVO_ONLY=0;
int VALIDATE=0;
int PRINT_DIST=0;
int TRUE_SS=0;

/* Performs the fragment superposition, returns the RMSD */
float super4(FRAGMENTS *One, FRAGMENTS *Two, int nres);
float Calc_dih(RESIDUE * OneP,RESIDUE * One,RESIDUE * OneA, float * temp);
BITS Atmrec;
FRAGMENTS *all;

typedef struct confidence
{
	float coil;
	float helix;
	float beta;
}CONF;

typedef struct Fragment 
{
	char line[250];
	struct Fragment *next;
}fragment;

int   A[600][MAXLEN];			  /* The Alignment Matrix 								   */
int   A_len[600][MAXLEN];		  /* The Alignment Length Matrix						   */
int   A_SS[600][MAXLEN];		  /* The Alignment SS score Matrix						   */
int   A_ramach[600][MAXLEN];	  /* The Alignment Ramach score Matrix	 				   */
float A_SSconf[600][MAXLEN];	  /* The Alignment SS Conf score Matrix					   */
fragment *Frag_lib[1000][100]; 	  /* A linked list for each query position (max=1000). Each sequence score is a key for the list (max=100). */
fragment *Frag_lib_rnd[1000][100]; 	  /* A linked list for each query position (max=1000). Each sequence score is a key for the list (max=100). */

/* Given a one-letter code for an aa, returns the index */
/* corresponding to that aa on the dictionary.          */
int find(char A,char *Dict)
{
  int i;
  for(i=0;i<strlen(Dict);i++)
    if(Dict[i]==A) 
      return i; 
  return strlen(Dict)+1;  
}

/* SS Alignment scoring function */
int score_ss(char A, char B) 
{
  if( A=='L' )
    A=' ';
//  printf("A: %c B: %c\n", A,B);
//  printf("Match: %d, Mismatch: %d\n", MATCH, MISMATCH);
//  printf("Score: %d\n", score(A,B));
	return score(A,B);
}

float score_ssconf(CONF confidence, char B) 
{

	if( B!='H' && B !='E' )
		return confidence.coil;
	if( B=='H' )
		return confidence.helix;
	if( B=='E' )
		return confidence.beta;
}

int main(int argc,char* argv[])
{
	short int LOGODS[7][20][20];  /* FREAD Environment Tables                              */
	short int Blossum[24][24];    /* BLOSSUM62 Env. Table                                  */
	int length,seq_score,ramach_score,ss_score,top=0, printed;
	int i,j,k,k2,total,start1,start2,index1,index2,start_res; /* Counters                  */
	int loop,helix,beta;		  /* Counters											   */
	int *Total;			 		  /* Total # of frags for each query position.	  		   */
	int *Total_rnd;		 		  /* Total # of frags for each query position.	  		   */
	int *Worst;					  /* Score of worst frag for each query position.		   */	
	int *Worst_rnd;				  /* Score of worst frag for each query position.		   */	
	int m;			 			  /* The length of the Query's fasta sequence              */
	int n;			  		      /* The length of the DB protein fasta sequence           */
  int num_header_lines;
	float Phi[MAXLEN], Psi[MAXLEN];
	char AUX[600],SS[2],Res[2];	  /* Multi-purpose Auxiliary Strings 		               */
	char c;						  /* Auxiliary Char 						               */
	char DB_Seq[MAXLEN];		  /* The DB protein fasta sequence                         */
  	char DB_SS[MAXLEN];		      /* The DB protein's secondary structure        		   */
	int DB_Ramach[MAXLEN];		  /* The DB protein's Ramachandran grouping.		       */	
	char Fasta_Seq[MAXLEN];		  /* The Query's fasta sequence                            */
	char Fasta_SS[MAXLEN];		  /* The Query's predicted secondary structure             */
    char Fasta_True_SS[MAXLEN];   /* The Query's true secondary structure                  */
	char Chain;                   /* The Protein chain.                            		   */
	char Header[82];              /* Information extracted from the header on the Database */
	char Query[82];
	char Dict2[21]={'G','A','V','L','M','I','F','Y','W','S','T','C','P','N','Q','K','R','H','D','E','\0'};              /* A one-letter code residue dictionary. */	
	char PATH[1000],pdb_file[1000],Query_Chain=0;
	int path_len;
	float rmsd,aux;
	float probability;
	float Angles[MAXLEN][4];
	float True_Angles[MAXLEN][4];
  float tPhi,tPsi;
	double resolution;			  /* The resolution of the crystal structure on the DB.	   */	
	double new_score,dist1,dist2;			  /* The predicted torsion angle score					   */	
    int min_length = 6;
    int max_length = 14;
	CONF Fasta_Conf[MAXLEN];      /* The confidence in the predicted secondary structure   */
	fragment *new_frag;
	FILE *input_fasta,*input_ss,*input_true_ss,*input_pdb,*input_phipsi,*input_true_phipsi,*blossum_file,*logods,*input_contact_file;
    int con1, con2, num_con, ct, c1, c2,flag,found_cb;
    double XORT[MAX_FIT_COORDS],YORT[MAX_FIT_COORDS],ZORT[MAX_FIT_COORDS];
    int Contacts[MAXLEN][2];
    int Flag[MAXLEN];
    float dist;

	srand (time(NULL));

    while (1)
    {
      static struct option long_options[] =
      {
          {"coevo_only",   no_argument, &COEVO_ONLY, 1},
          {"true_ss",      no_argument, &TRUE_SS, 1},
          {"min_length",   required_argument,       0, 'l'},
          {"max_length",   required_argument,       0, 'L'},
          {"max_frags",    required_argument,       0, 'M'},
          {"help",         no_argument,             0, 'h'},
          {"contact_file", required_argument,       0, 'C'},
          {"validate",     required_argument,       0, 'v'},
          {"chain",        required_argument,       0, 'N'},
          {"print_dist",   no_argument, &PRINT_DIST, 1},
          {"input_file",   required_argument,       0, 'i'},
          {0, 0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;

      c = getopt_long (argc, argv, "b:e:O:l:L:M:h:C:v:N:i:",
                       long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1)
        break;

      switch (c)
      {
         case 0:
              if (long_options[option_index].flag != 0)
                break;

        case 'l':
          min_length = atoi(optarg);
          break;

        case 'L':
          max_length = atoi(optarg);
          break;

        case 'M':
          MAX_FRAG = atoi(optarg);
          break;

        case 'C':
          input_contact_file = fopen(optarg,"r");
          if (input_contact_file == NULL) {printf("Predicted contact file not found: %s\n", optarg); return 0;} 
          COEVO=1;
          break;

        case 'v':
          strcpy(pdb_file,optarg);
          VALIDATE=1;
          break;

        case 'N':
            Query_Chain=optarg[0];
            break;

        case 'i':
          strcpy(Query,optarg);
          break;

        case 'h':
            printf("\n\nUSAGE: %s [-options] -i PDB_ID\n\n",argv[0]);
            printf("\t-M, --max_frags int_value [3000] --- The number of candidate fragments considered.\n");
            printf("\t-l, --min_length int_value [6]  --- The minimum fragment length\n");
            printf("\t-L, --max_length int_value [14] --- The maximum fragment length\n");
            printf("\t-C, --contact_file filename --- A predicted contact file in three-column format (res1,res2,score). \n");
            printf("\t-v, --validation filename --- A PDB file for your target if benchmarking/validating fragments. \n");
            printf("\t-N, --chain char [A] --- If validating fragments, the protein chain of your target. \n");
            printf("\t--print-dist --- If using co-evolution fragments, output the distance between the residues predicted to be in contact.\n");
            printf("\t--true_ss    --- Uses the true secondary structure, as output by DSSP.\n");
            printf("\t--coevo_only --- Only output fragments that satisfy predicted contacts \n"); 
            printf("\n\nDESCRIPTION:\n\n\tFlib Co-evo v1.0\n\n");
            break;

        default:
            abort ();
      }
    }

	if(strlen(Query)<1)
	{
        fprintf(stderr,"\n\nUSAGE: %s [-options] -i PDB_ID\n\nUse '-help' to print a detailed description of command line options\n\n",argv[0]);
		return 0;
	}

    if(VALIDATE && !Query_Chain)
    {
        fprintf(stderr,"[WARNING] No chain provided for fragment validation. Using chain A as default. \n\t  Please, make sure that this is correct.\n");
        Query_Chain='A';    
    }

    if(TRUE_SS && !Query_Chain)
    {
        fprintf(stderr,"[WARNING] No chain provided for true secondary structure. Using chain A as default. \n\t  Please, make sure that this is correct.\n");
        Query_Chain='A';
    }

    if( !VALIDATE && Query_Chain && !TRUE_SS )
    {
        fprintf(stderr,"[WARNING] Chain argument will be ignored.\n");
    }

    if(getenv("FLIB")==NULL) 
    {
        printf("[ERROR] Please, set the environment variable FLIB to point to the directory where Flib was installed. e.g.\n\nsh: export FLIB=/path/to/Flib/\n\n");
        return 0;
    }

    if(getenv("PDB")==NULL)
    {
        printf("[ERROR] Please, set the environment variable PDB to point to your local copy of the PDB. e.g.\n\nsh: export PDB=/path/to/PDB/\n\n");
        return 0;
    }

  	/*****      FILE HANDLING   *****/
	strcpy(AUX,Query);
	input_fasta = fopen(strcat(AUX,".fasta.txt"),"r");
	if (input_fasta == NULL) {printf("Fasta file not found: %s\n", AUX); return 0;}	
	
	strcpy(AUX,Query);
	input_ss = fopen(strcat(AUX,".fasta.ss"),"r");
	if (input_ss == NULL) {printf("Predicted Secondary Structure file not found: %s\n", AUX); return 0;}
  //skip header
  num_header_lines=0;
  while(fgetc(input_ss) == '#')
  {
    for(c=getc(input_ss);c!='\n' && c!= EOF;c=getc(input_ss));
    num_header_lines++;
  }
  fclose(input_ss);
	strcpy(AUX,Query);
  input_ss = fopen(strcat(AUX,".fasta.ss"),"r");

  while(num_header_lines > 0) 
  {
    fgets(AUX, 600, input_ss);
    num_header_lines--;
  }
    
    if(TRUE_SS)
    {
        sprintf(AUX, "validator_%s.dssp", Query);
  //      printf("Using true SS file: %s\n", AUX);
        input_true_ss = fopen(AUX,"r");
        if (input_true_ss == NULL) {printf("True Secondary Structure file not found: %s\n", AUX); return 0;}
    }

	sprintf(PATH,"%s/data/parsedPDB_new.txt",getenv("FLIB"));
	input_pdb = fopen(PATH,"r");
	if (input_pdb == NULL) {printf("Protein database file is missing: parsedPDB_new.txt\n"); return 0;}

	strcpy(AUX,Query);
    	input_phipsi = fopen(strcat(AUX,".spd3"),"r");
	if (input_phipsi == NULL) {printf("SPIDER2 input file (predicted torsion angles) not found!\n"); return 0;}

  if(TRUE_SS)
  {
    sprintf(AUX, "validator_%s.angles", Query);
    input_true_phipsi = fopen(AUX,"r");
    if (input_true_phipsi == NULL) {printf("Validator angles file not found!\n"); return 0;}
  }

	sprintf(PATH,"%s/data/blossum62.txt",getenv("FLIB"));
	blossum_file = fopen(PATH,"r");
	if (blossum_file == NULL) {printf("BLOSSUM file not found!\n"); return 0;}
	/***** END OF FILE HANDLING *****/

	/***** READ QUERY FASTA SEQUENCE *****/
	/* Remove Header from FASTA file */
	for(c = fgetc(input_fasta); c!='\n' ; c=fgetc(input_fasta)); 
	/* Read the sequence and store it in the Fasta_Seq string */
	for(fscanf(input_fasta,"%s",Fasta_Seq); fscanf(input_fasta,"%s",AUX)!=EOF ; strcat(Fasta_Seq,AUX) );
	/* Compute the length of the sequence and store it in m */
	m = strlen(Fasta_Seq);
	/***** END OF READ QUERY FASTA SEQUENCE *****/

	/***** READ QUERY'S PREDICTED SECONDARY STRUCTURE SEQUENCE *****/
	/* READ INPUT SECONDARY STRUCTURE SEQUENCE */

	while(fscanf(input_ss,"%d",&i)!=EOF && i<= m  )
		fscanf(input_ss,"%s %c %f %f %f",AUX,&Fasta_SS[i-1],&Fasta_Conf[i-1].coil,&Fasta_Conf[i-1].helix,&Fasta_Conf[i-1].beta);
	if(i!=m) { printf("ERROR: Fasta sequence and pred. secondary structure have different lengths!\n"); return 0; }
	Fasta_SS[m]='\0'; 
	/***** END OF READ QUERY'S PREDICTED SECONDARY STRUCTURE SEQUENCE *****/
  //      printf("%s\n",Fasta_SS);
    /***** IF SPECIFIED AT INPUT, READ QUERY'S TRUE SECONDARY STRUCTURE SEQUENCE *****/
    if(TRUE_SS)
    {   
  //      printf("%s\n",Fasta_SS);
        while(1)
        {
            /* Disregard DSSP header */
            for(c=getc(input_true_ss);c!='#' && c!= EOF;c=getc(input_true_ss));
            if(c==EOF)
                break;
            fscanf(input_true_ss,"%s",AUX);
            if(!strcmp(AUX,"RESIDUE"))
            {
                for(c=fgetc(input_true_ss);c!='\n' && c!=EOF;c=fgetc(input_true_ss));
                break;
            }
            
        }
        for(i=0;i<=m;i++) Fasta_True_SS[i]=0;

        /* Start reading true SS */
        for(c2=0; c2<m; c2++)
        {
            Fasta_True_SS[c2]='-';
        }
        Fasta_True_SS[c2]='\0';

        for(;fscanf(input_true_ss,"%d %s %c %s",&c1,&AUX,&Chain,&Res) != EOF;)
        {   
             if(AUX[0]!='!' && Chain == Query_Chain)
             {
               /* Assign real SS only to resolved residues */
                sscanf(AUX, "%d", &c2);
                c2--;
                fgetc(input_true_ss);
                fgetc(input_true_ss);
                Fasta_True_SS[c2]=fgetc(input_true_ss);    
                if(Fasta_True_SS[c2]==' ')
                    Fasta_True_SS[c2]='C';
             }
             for(c=fgetc(input_true_ss);c!='\n' && c!=EOF;c=fgetc(input_true_ss)); 
        }
        
    //    printf("%s\n",Fasta_True_SS);

        /* Copy the  relevant section of the true SS across */
        for(i=0; i<m; i++)
            if(Fasta_True_SS[i] != '-')
                Fasta_SS[i] = Fasta_True_SS[i];
        
  //      printf("%s\n",Fasta_SS);
    }

    /***** END OF READ QUERY'S TRUE SECONDARY STRUCTURE SEQUENCE *****/


	/***** READ QUERY'S PREDICTED TORSION ANGLES *****/
	/* Remove Header from SPIDER2 output file */
	for(c=fgetc(input_phipsi);c!='\n';c=fgetc(input_phipsi));

	for(i=0;fscanf(input_phipsi,"%d %s %s %f %f %f",&k,Res,SS,&aux,&Angles[i][0],&Angles[i][1])!=EOF;i++){
		for(c=fgetc(input_phipsi);c!='\n';c=fgetc(input_phipsi));
  }


        if(i!=m) { printf("ERROR: Fasta sequence and pred. torsion angles have different lengths!\n"); return 0; }
	
	
    if(TRUE_SS)
    {   
        /* Start reading true torsion angles */
        // Remove header from angles file
	      for(c=fgetc(input_true_phipsi);c!='\n';c=fgetc(input_true_phipsi));
        for(c2=0; c2<m; c2++)
        {
          for(j=0; j<2; j++)
            True_Angles[c2][j]=-999;
        }

        for(;fscanf(input_true_phipsi,"%d %c %f %f",&c2,Res,&tPhi,&tPsi) != EOF;)
        {   
               /* Assign real TA only to resolved residues */
                c2--;
                True_Angles[c2][0]=tPhi;
                True_Angles[c2][1]=tPsi;
         //       printf("%d %f %f\n", c2, tPhi, tPsi);
             for(c=fgetc(input_true_phipsi);c!='\n' && c!=EOF;c=fgetc(input_true_phipsi)); 
        }
        
        /* Copy relevant section of the true torsion angles across */
        for(i=0; i<m; i++){
          for(j=0; j<2; j++)
            if(True_Angles[i][j] != -999)
                Angles[i][j] = True_Angles[i][j];
      //    printf("%d %f %f\n", i, Angles[i][0], Angles[i][1]);
          }
    }
        /***** END OF READ QUERY'S PREDICTED TORSION ANGLES *****/

	/***** READ BLOSSUM MATRIX FOR THE ALIGNMENT *****/
	for(i=0;i<24;i++)
		for(j=0;j<24;j++)
			fscanf(blossum_file,"%hd",&Blossum[i][j]);  

	/***** END OF READ BLOSSUM MATRIX FOR THE ALIGNMENT *****/
	
	/***** READ FREAD ENV. TABLES *****/
	sprintf(AUX,"%s/data/logods0.txt",getenv("FLIB"));
	path_len = strlen(AUX);
	for(;AUX[path_len-5]<48+7;AUX[path_len-5]++)
	{	
		logods=fopen(AUX,"r");
		if (logods == NULL) {printf("%s file for FREAD environment matrix not found\n",AUX); return 0;}

		for(i=0;i<20;i++)
			for(j=0;j<20;j++)
				fscanf(logods,"%hd",&LOGODS[(int)AUX[path_len-5]-48][i][j]);
		fclose(logods);
	}
	/***** END OF READ FREAD ENV. TABLES *****/

    /***** IF PROVIDED, READ CONTACT FILE *****/
    if(COEVO)
    {
        for(i=0; fscanf(input_contact_file,"%d %d %lf",&Contacts[i][0],&Contacts[i][1],&aux)!=EOF; i++);
        num_con=i;

    }

	/***** INITIALIZE THE FRAGMENT LIBRARY LINKED LIST *****/
	Total=malloc(sizeof(int)*m);
	Total_rnd=malloc(sizeof(int)*m);
	Worst=malloc(sizeof(int)*m);
	Worst_rnd=malloc(sizeof(int)*m);
	for(i=0;i<m;i++)
	{
		Total[i]=0;
		Total_rnd[i]=0;
		Worst[i]=0;
		Worst_rnd[i]=-1;
		for(j=0;j<100;j++)
		{
			Frag_lib[i][j]=NULL;
			Frag_lib_rnd[i][j]=NULL;
		}
	}

	/***** END OF INITIALIZE THE FRAGMENT LIBRARY LINKED LIST *****/

	/***** BEGINNING OF FRAGMENT EXTRACTION *****/

	probability = (3.0*(float)(m)*10000)/(120916.0/4); /* Define a probability for random fragment extraction */
	/* Iterate over every sequence in the input database */
	for(total=0; fscanf(input_pdb,"%s %s %c %lf",AUX,Header,&Chain,&resolution) != EOF  ;total++)
	{
		/* Read DB protein fasta sequence */
		fscanf(input_pdb,"%s",DB_Seq);			/* Extract the protein sequence from the database. */	
		fgetc(input_pdb); 				/* Get the line break char (chomp) 				   */
		n=strlen(DB_Seq);				/* Get the length of the DB protein sequence 	   */

		//printf("%s %d\n",DB_Seq,n);
		
		/* Read DB protein secondary structure  */
		for(j=0 , DB_SS[j]= fgetc(input_pdb); DB_SS[j]!=(int)'\n' && DB_SS[j] != EOF ; j++, DB_SS[j]= fgetc(input_pdb));
		DB_SS[j]='\0';
		if(j != n)
		{
			for(c=fgetc(input_pdb);c!='>';c=fgetc(input_pdb));
			ungetc(c,input_pdb);
			fprintf(stderr,"Sequence and secondary structure uncompatible lengths for protein: %s\t Seq: %d\tSS: %d\n",Header,n,j);
			continue;
		}
		
		/* Read DB protein Ramachandran Groups */
		for(j=0 , DB_Ramach[j]= (int)fgetc(input_pdb); DB_Ramach[j]!=(int)'\n' && DB_Ramach[j] != EOF ; j++, DB_Ramach[j]= (int)fgetc(input_pdb))
				DB_Ramach[j]-=48;	
		if(j!= n)
		{
			for(c=fgetc(input_pdb);c!='>';c=fgetc(input_pdb));
			ungetc(c,input_pdb);
			fprintf(stderr,"Sequence and Ramachandran Grouping uncompatible lengths for protein: %s\t Seq: %d\tRamach: %d\n",Header,n,j);
			continue;
		}

		/***** RANDOM EXTRACTION *****/
		#if RANDOM
		for(j=probability;j >= 0 && n > max_length; j--)
		{	
			if(!j && (float)(rand()%10000)/10000.0 > probability-(int)probability)
		    		continue;

			length = rand() % (max_length - min_length + 1) + min_length; /* Randomize fragment length between min_length and max_length residues */
			start1  = rand() % (m-length + 1); /* Randomize fragment starting position in query between 0 and (length of protein - length of fragment) */
			start2  = rand() % (n-length + 1); /* Randomize fragment starting position in pdb structure between 0 and (length of protein - length of fragment) */
			loop=0; helix=0; beta=0; ss_score=0; seq_score=0; ramach_score=0;

			/* Compute sequence and secondary structure score for the fragment */
			for (k=start1,k2=start2; k-start1 < length; k++,k2++)
			{
				index1=find(Fasta_Seq[k],Dict2);
				index2=find(DB_Seq[k2],Dict2);
				if((DB_Ramach[start2] >= 7 || DB_Ramach[start2] < 0 ) || (index1 < 0 || index1 >= 20 || index2 < 0 || index2 >=20 ))
				{
					ramach_score=-50;
					ss_score=-50;
					break;
				}
				else
					ramach_score += LOGODS[DB_Ramach[start2]][index1][index2];
		       		ss_score += score_ss(Fasta_SS[k],DB_SS[k2]);

				/* Find predominant SS on the fragment */
				switch(DB_SS[k2])
				{
					case 'H':
						helix++;
						break; 
					case 'E':
						beta++;
						break;
					default:
						loop++;
						break; 
				}
			}
			
			if( ( ramach_score > 6 || ( helix < (length/2+1) && beta < (length/2+1) && ramach_score > 0 )) && (Total_rnd[start1] < MAX_FRAG || ramach_score > Worst_rnd[start1] ) )
			{	


				if( (!( helix < length/2+1 && beta < length/2 + 1 ))  && ss_score < 0 )
					continue; 

				new_frag=malloc(sizeof(fragment));
				if(new_frag == NULL)
				{
					fprintf(stderr,"Malloc error!\n");
					break;
				}
				if(ramach_score > 99) ramach_score = 99; /* This should not happen very often! */

				new_frag->next = Frag_lib_rnd[start1][ramach_score];
         
				sprintf(new_frag->line,"%s\t%c\t%3d\t%3d\t",Header,Chain,start2,start2+length-1);

				for (k=0; k < length ; k++)
					AUX[k]=DB_Seq[start2+k];
				AUX[k]='\0';
				strcat(new_frag->line,AUX);

				if(helix>=length/2+1)
					strcat(new_frag->line,"\tH");
				else { 	if(beta>=length/2+1)
					strcat(new_frag->line,"\tB");
				else { 	if(loop>=length/2+1)
					strcat(new_frag->line,"\tL");
				else 
					strcat(new_frag->line,"\tO"); } }

				sprintf(AUX,"\t%d\t%d\t%d\t%d\t%.2lf\t%d\n",ramach_score,0,length,start1,resolution,ss_score);
				strcat(new_frag->line,AUX);
				Frag_lib_rnd[start1][ramach_score] = new_frag;

				if(Total_rnd[start1]==0) /* First fragment for that position! */
						Worst_rnd[start1]=ramach_score; /* Update the score of the Worst Fragment for that position */	

				if(Total_rnd[start1] < MAX_FRAG)
				{
					Total_rnd[start1]++;	/* Update the number of fragments for that position. */
					if(ramach_score < Worst_rnd[start1]) Worst_rnd[start1]=ramach_score; /* Update the score of the Worst Fragment for that position */
				}
				else 
				{
					/* Remove the worst one */
					if(Frag_lib_rnd[start1][Worst_rnd[start1]] != NULL)
					{
						new_frag = Frag_lib_rnd[start1][Worst_rnd[start1]];
						Frag_lib_rnd[start1][Worst_rnd[start1]] = new_frag->next;	
						free(new_frag);
					}
					if( Frag_lib_rnd[start1][Worst_rnd[start1]] == NULL )
						for( ; Worst_rnd[start1] < 100; Worst_rnd[start1]++)
							if(Frag_lib_rnd[start1][Worst[start1]]!=NULL)
								break;
				}	

			}
		}
		#endif

		/***** EXHAUSTIVE EXTRACTION *****/
		#if EXHAUSTIVE	
		/***** PERFORM THE SEQUENCE ALIGNMENT *****/		

		/* Initialize the Scoring Matrix */
		for(i=0;i<=m;i++) {
			A[i][0]=0; 		A_len[i][0]=0; A_SS[i][0]=0; A_SSconf[i][0]=0; A_ramach[i][0]=0;	}
		for(j=0;j<=n;j++) {	
			A[0][j]=0;		A_len[0][j]=0; A_SS[0][j]=0; A_SSconf[0][j]=0; A_ramach[0][j]=0;	}	

		/* Populate the scoring matrix using sequence score! */

		for(i=1;i<=m;i++)
		{
			for(j=1;j<=n;j++)
    		{
				
   				if(DB_Ramach[j-1] >= 0 && DB_Ramach[j-1] < 7 )
				{
					index1=find(Fasta_Seq[i-1],Dict2);
					index2=find(DB_Seq[j-1],Dict2);
					if(index1 < 0 || index1 >= 20 || index2 < 0 || index2 >=20 )
					 	A_ramach[i][j]=0;
					else
						A_ramach[i][j] = max( A_ramach[i-1][j-1]+LOGODS[DB_Ramach[j-1]][index1][index2] , 0 );
				}
				else
					A_ramach[i][j] = 0;
				
				if( !A_ramach[i][j] )
				{
 					 A[i][j]       = 0;
					 A_len[i][j]   = 0;
					 A_SS[i][j]    = 0;
					 A_SSconf[i][j]= 0.0;

				}
				else
				{
					 A[i][j]	   = 0;
					 A_len[i][j]   = A_len[i-1][j-1]   + 1;
					 A_SS[i][j]    = A_SS[i-1][j-1]    + score_ss(Fasta_SS[i-1],DB_SS[j-1]);
					 A_SSconf[i][j]= A_SSconf[i-1][j-1]+ score_ssconf(Fasta_Conf[i-1],DB_SS[j-1]);
				}
			} 
		}


		/* "Traceback" step of sequence alignment */   
		for(i=m;i>0;i--)
			for(j=n;j>0;j--)
			{

				/* If this fragment has a better score than the worst fragment we extracted so far 
				   for that position and if the length of the fragment is greater than min_length, we extract it!     */
				if(A_len[i][j] >= min_length && A_len[i][j] <= max_length)
				{
					start1 = i - A_len[i][j];
					start2 = j - A_len[i][j];
				}
				else
					continue;

				if( start1 < 0 || start1 >= m || start2 < 0 || start2 >=n ) continue;				
							
				if( (Total[start1] < MAX_FRAG || A_ramach[i][j] > Worst[start1] )  && A_len[i][j] >= min_length && A_SS[i][j] >= 0 && A_len[i][j] <= max_length)
				{	
					loop=0; helix=0; beta=0;
					/* Find predominant SS on the fragment */
					for (k=0; k < A_len[i][j] && k + start2 < n; k++)
					{
						AUX[k]=DB_Seq[start2+k];
						switch(DB_SS[start2+k])
						{
								case 'H':
									helix++;
									break; 
								case 'E':
									beta++;
									break;
								default:
									loop++;
									break; 
						}
					}
					AUX[k]='\0';


					if(A_ramach[i][j] > 99) A_ramach[i][j] = 99; /* This should not happen very often! */

					/* Allocate the fragment */
					new_frag = malloc(sizeof(fragment));
					new_frag->next = Frag_lib[start1][A_ramach[i][j]];
					sprintf(new_frag->line,"%s\t%c\t%3d\t%3d\t",Header,Chain,start2,start2+A_len[i][j]-1);


					strcat(new_frag->line,AUX);
					if(helix>=A_len[i][j]/2+1)
						strcat(new_frag->line,"\tH");
					else { 	if(beta>=A_len[i][j]/2+1)
						strcat(new_frag->line,"\tB");
					else { 	if(loop>=A_len[i][j]/2+1)
						strcat(new_frag->line,"\tL");
					else 
						strcat(new_frag->line,"\tO"); } }

					sprintf(AUX,"\t%d\t%d\t%d\t%d\t%.2lf\t%d\n",A_ramach[i][j],A[i][j],A_len[i][j],start1,resolution, A_SS[i][j]);
					strcat(new_frag->line,AUX);

					Frag_lib[start1][A_ramach[i][j]] = new_frag;
					/* If the number of fragments for that position is less than MAX_FRAG... */

					if(Total[start1]==0) /* First fragment for that position! */
							Worst[start1]=A_ramach[i][j]; /* Update the score of the Worst Fragment for that position */	

					if(Total[start1]<MAX_FRAG)
					{
						Total[start1]++;	/* Update the number of fragments for that position. */
						if(A_ramach[i][j] < Worst[start1]) Worst[start1]=A_ramach[i][j]; /* Update the score of the Worst Fragment for that position */
					}
					else 
					{
							/* Remove the worst one */
							if(Frag_lib[start1][Worst[start1]] != NULL)
							{
								new_frag = Frag_lib[start1][Worst[start1]];
								Frag_lib[start1][Worst[start1]] = Frag_lib[start1][Worst[start1]]->next;	
								free(new_frag);
							}
							if( Frag_lib[start1][Worst[start1]] == NULL )
								for( ; Worst[start1] < 100; Worst[start1]++)
									if(Frag_lib[start1][Worst[start1]]!=NULL)
										break;
					}	
				}
			} 
		#endif			
	}
	fprintf(stderr,"Fragment Extraction finished.\n");

	all=malloc(sizeof(FRAGMENTS));
	for(i=0;i<MaxRes;i++)
	{
		all->res[i]= (RESIDUE *) malloc(sizeof(RESIDUE));
		for(j=0;j<MaxAtom;j++)
			all->res[i]->atom[j]=(ATOM *) malloc(sizeof(ATOM));
	}


	/*** FRAGMENT VALIDATION ****/
	if ( (Atmrec.frag[0]= (FRAGMENTS *) calloc(1,sizeof(FRAGMENTS))) == NULL ) { printf("Error in malloc\n");}
	if ( (Atmrec.frag[1]= (FRAGMENTS *) calloc(1,sizeof(FRAGMENTS))) == NULL ) { printf("Error in malloc\n");}

	for(i=0;i<m-(min_length-1);i++)
	{
		#if EXHAUSTIVE
		for(j=Worst[i];j<100;j++)
		{	
			for(new_frag= Frag_lib[i][j];new_frag!=NULL;new_frag=new_frag->next)
			{

				top=0;
				sscanf(new_frag->line,"%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%lf\t%d",Header,&Chain,&start1,&k,DB_Seq,&c,&ramach_score,&seq_score,&length,&start2,&resolution,&ss_score);
				//printf("%d\n",length);
				//printf("%s\t%c\t%d\t%d\t%s\t%c\n",Header,Chain,start1,length,DB_Seq,c);
				Atmrec.frag[0]->start_res	= start1;
				Atmrec.frag[1]->start_res	= start2;
				sprintf(Atmrec.frag[0]->fname,"%s/%c%c%c%c.pdb",getenv("PDB"),tolower(Header[0]),tolower(Header[1]),tolower(Header[2]),tolower(Header[3]));
				top = coordcol(Atmrec.frag[top]->fname,top,length,Atmrec.frag[top]->start_res,Chain,DB_Seq);

				if(VALIDATE)
				{
					sprintf(Atmrec.frag[1]->fname,"./%s",pdb_file);
					top = coordcol(Atmrec.frag[top]->fname,top,length,Atmrec.frag[top]->start_res,Query_Chain,NULL);
				}
				if(top<1+VALIDATE) /* If, by any reason, was not capable of finding the second fragment! */
					fprintf(stderr,"Could not open the fragment starting at position %d for %s_%c.\n",start1,Header,Chain);
				else
				{
                    /* If validating, calculate the RMSD between fragment and target */
					if(VALIDATE)
					{
						if( Atmrec.frag[0]->natom == Atmrec.frag[1]->natom  &&  Atmrec.frag[1]->natom > 0  &&  Atmrec.frag[0]->natom > 0 && Atmrec.frag[0]->numres > 0 && Atmrec.frag[1]->numres > 0 && Atmrec.frag[0]->numres == Atmrec.frag[1]->numres )
						{
							rmsd =  super4( Atmrec.frag[0], Atmrec.frag[1],length);
							if(rmsd > 10.0 || rmsd < 0.0) continue;
							sprintf(AUX,"%.2lf",rmsd);
							new_frag->line[strlen(new_frag->line)-1]='\t';
							strcat(new_frag->line,AUX);	
						}
					}
					if( coordcol2(all,Header,Chain,DB_Seq,&start_res) ) 
						continue;

                    /* If contacts were provided, check if any contacts occur within the fragment. */
                    if(COEVO)
                    {
                        flag=0;
                        for(ct=0;ct<num_con;ct++)
                        {
                            con1=Contacts[ct][0]-1-start2;
                            con2=Contacts[ct][1]-1-start2;
                            if(con1<0 || con1 >= length || con2<0 || con2 >= length)
                                continue;

                            flag=1;
                            /* If it got here, it means that a contact exists within the fragment */
                            for(c1=0;c1<length;c1++)
                            {
                                found_cb=0;
                                for(c2=0;c2<Atmrec.frag[0]->res[c1]->numatom;c2++)
                                {
                                    if(!strcmp(Atmrec.frag[0]->res[c1]->atom[c2]->atomname,"CB"))
                                    {
                                        XORT[c1] = Atmrec.frag[0]->res[c1]->atom[c2]->x;
                                        YORT[c1] = Atmrec.frag[0]->res[c1]->atom[c2]->y;
                                        ZORT[c1] = Atmrec.frag[0]->res[c1]->atom[c2]->z;
                                        found_cb=1;
                                        break;
                                    }
                                }
                                if(!found_cb)
                                {
                                   for(c2=0;c2<Atmrec.frag[0]->res[c1]->numatom;c2++)
                                   {
                                       if(!strcmp(Atmrec.frag[0]->res[c1]->atom[c2]->atomname,"CA"))
                                       {
                                           XORT[c1] = Atmrec.frag[0]->res[c1]->atom[c2]->x;
                                           YORT[c1] = Atmrec.frag[0]->res[c1]->atom[c2]->y;
                                           ZORT[c1] = Atmrec.frag[0]->res[c1]->atom[c2]->z;
                                           break;
                                       }
                                   }
                                }
                            }
                            dist = sqrt(pow((XORT[con1]-XORT[con2]),2) + pow((YORT[con1]-YORT[con2]),2) + pow((ZORT[con1]-ZORT[con2]),2));
                            if( dist < 8.0)
                                break;
                        }
                    }
                 

                    /* Calculate the torsion angle score */
					new_score=0.0;
					for(k=0, k2=start_res ; k2 < start_res+length ; k++, k2++)
					{
						if(k2<1 || k2 >= all->numres-1)
							continue;
						Phi[k] = Calc_dih( all->res[k2-1] , all->res[k2] , all->res[k2+1], &Psi[k]);
                                                dist1= min(fabs(Phi[k]-Angles[start2+k][0]), fabs(fabs(Phi[k])+fabs(Angles[start2+k][0]) - 360.0) );
                                                dist2= min(fabs(Psi[k]-Angles[start2+k][1]), fabs(fabs(Psi[k])+fabs(Angles[start2+k][1]) - 360.0) );


                                                if(dist1 != dist1 || dist2 != dist2 || fabs(dist1) > 180 || fabs(dist2) > 180)
						{
							new_score=9999.0;
							break;
						}						
						new_score+= dist1+dist2; 
					}
                    printed=0;
					if(!VALIDATE)
					{
							new_frag->line[strlen(new_frag->line)-1]='\t';
                            if(!COEVO || (COEVO && !COEVO_ONLY && ( !flag || dist < 8.0 ) ) || (COEVO && COEVO_ONLY && flag  && dist < 8.0) )
                            {
    							printf("%s%.2lf",new_frag->line,new_score);	
                                printed=1;
                            }
					}
					else
                    {
                        if(!COEVO || (COEVO && !COEVO_ONLY && ( !flag || dist < 8.0 ) ) || (COEVO && COEVO_ONLY && flag && dist < 8.0 ) )
                        {
    						printf("%s\t%.2lf",new_frag->line,new_score);	
                            printed=1;
                        }
                    }
                    if(PRINT_DIST && printed)
                    {
                        if(flag)
                            printf("\t%d\t%d\t%d\t%d\t%.2f\n",Contacts[ct][0],Contacts[ct][1],con1+1,con2+1,dist);
                        else
                            printf("\t-1.0\n");
                    }
                    else
                    {   
                        if(printed)
                            printf("\n");
                    }
				}
			}
		}
		#endif

		#if RANDOM
		for(j=99;j>=0;j--)
		{
			for(new_frag= Frag_lib_rnd[i][j];new_frag!=NULL;new_frag=new_frag->next)
			{
				top=0;
				sscanf(new_frag->line,"%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%lf\t%d",Header,&Chain,&start1,&k,DB_Seq,&c,&ramach_score,&seq_score,&length,&start2,&resolution,&ss_score);
				Atmrec.frag[0]->start_res	= start1;
				Atmrec.frag[1]->start_res	= start2;

				sprintf(Atmrec.frag[0]->fname,"%s/%c%c%c%c.pdb",getenv("PDB"),tolower(Header[0]),tolower(Header[1]),tolower(Header[2]),tolower(Header[3]));
				top = coordcol(Atmrec.frag[top]->fname,top,length,Atmrec.frag[top]->start_res,Chain,DB_Seq);

				if(VALIDATE)
				{				
                    strcpy(Atmrec.frag[1]->fname,pdb_file);
					top = coordcol(Atmrec.frag[top]->fname,top,length,Atmrec.frag[top]->start_res,Query_Chain,NULL);
				}

				if(top<1+VALIDATE) /* If, by any reason, was not capable of finding the second fragment! */
					fprintf(stderr,"Could not open the fragment starting at position %d for %s_%c.\n",start1,Header,Chain);
				else
				{
					if(VALIDATE)
					{
						if( Atmrec.frag[0]->natom == Atmrec.frag[1]->natom  &&  Atmrec.frag[1]->natom > 0  &&  Atmrec.frag[0]->natom > 0 && Atmrec.frag[0]->numres > 0 && Atmrec.frag[1]->numres > 0 && Atmrec.frag[0]->numres == Atmrec.frag[1]->numres )
						{
							rmsd =  super4( Atmrec.frag[0], Atmrec.frag[1],length);
							if(rmsd > 10.0 || rmsd < 0.0) continue;
							sprintf(AUX,"%.2lf",rmsd);
							new_frag->line[strlen(new_frag->line)-1]='\t';
							strcat(new_frag->line,AUX);	
						}					
					}	
					
					if( coordcol2(all,Header,Chain,DB_Seq,&start_res) ) 
						continue;


                    /* If contacts were provided, check if any contacts occur within the fragment. */
                    if(COEVO)
                    {
                        flag=0;
                        for(ct=0;ct<num_con;ct++)
                        {
                            con1=Contacts[ct][0]-1-start2;
                            con2=Contacts[ct][1]-1-start2;

                            if(con1<0 || con1 >= length || con2<0 || con2 >= length)
                                continue;

                            flag=1;
                            /* If it got here, it means that a contact exists within the fragment */
                            for(c1=0;c1<length;c1++)
                            {
                                found_cb=0;
                                for(c2=0;c2<Atmrec.frag[0]->res[c1]->numatom;c2++)
                                {
                                    if(!strcmp(Atmrec.frag[0]->res[c1]->atom[c2]->atomname,"CB"))
                                    {
                                        XORT[c1] = Atmrec.frag[0]->res[c1]->atom[c2]->x;
                                        YORT[c1] = Atmrec.frag[0]->res[c1]->atom[c2]->y;
                                        ZORT[c1] = Atmrec.frag[0]->res[c1]->atom[c2]->z;
                                        found_cb=1;
                                        break;
                                    }
                                }
                                if(!found_cb)
                                {
                                   for(c2=0;c2<Atmrec.frag[0]->res[c1]->numatom;c2++)
                                   {
                                       if(!strcmp(Atmrec.frag[0]->res[c1]->atom[c2]->atomname,"CA"))
                                       {
                                           XORT[c1] = Atmrec.frag[0]->res[c1]->atom[c2]->x;
                                           YORT[c1] = Atmrec.frag[0]->res[c1]->atom[c2]->y;
                                           ZORT[c1] = Atmrec.frag[0]->res[c1]->atom[c2]->z;
                                           break;
                                       }
                                   }
                                }
                            }
                            dist = sqrt(pow((XORT[con1]-XORT[con2]),2) + pow((YORT[con1]-YORT[con2]),2) + pow((ZORT[con1]-ZORT[con2]),2));
                            if( dist < 8.0)
                                break;
                        }
                    }

                    /* Calculate the torsion angle score */
					new_score=0.0;
					for(k=0, k2=start_res ; k2 < start_res+length && k2 < all->numres; k++, k2++)
					{

						if(k2<1)
							continue;
						Phi[k] = Calc_dih( all->res[k2-1] , all->res[k2] , all->res[k2+1], &Psi[k]);
						dist1= min(fabs(Phi[k]-Angles[start2+k][0]), abs(fabs(Phi[k])+fabs(Angles[start2+k][0]) - 360) );
                                                dist2= min(fabs(Psi[k]-Angles[start2+k][1]), abs(fabs(Psi[k])+fabs(Angles[start2+k][1]) - 360) );


						if( dist1 != dist1 || dist2 != dist2 || fabs(dist1) > 180.00 || fabs(dist2) > 180.00)
						{
							new_score=9999.0;
							break;
						}
						new_score += dist1+dist2;
					}

                    printed=0;
					if(!VALIDATE)
					{
							new_frag->line[strlen(new_frag->line)-1]='\t';
                            if(!COEVO || (COEVO && !COEVO_ONLY && ( !flag || dist < 8.0 ) ) || (COEVO && COEVO_ONLY && flag && dist < 8.0 ) )
                            {
                                printed=1;
                                printf("%s%.2lf",new_frag->line,new_score);
                            }
					}
					else
                    {
                        if(!COEVO || (COEVO && !COEVO_ONLY && ( !flag || dist < 8.0 ) ) || (COEVO && COEVO_ONLY && flag && dist < 8.0 ) )
                        {
                            printed=1;
                            printf("%s\t%.2lf",new_frag->line,new_score);
                        }
                    }
                    if(PRINT_DIST && printed)
                    {   
                        if(flag)
                            printf("\t%d\t%d\t%d\t%d\t%.2f\n",Contacts[ct][0],Contacts[ct][1],con1+1,con2+1,dist);
                        else
                            printf("\t-1.0\n");
                    }
                    else
                        if(printed)
                            printf("\n");
				}
			}
		}
		#endif
	}

	free(Atmrec.frag[0]);
	free(Atmrec.frag[1]);

	free(Total);
	free(Total_rnd);
	free(Worst);
	free(Worst_rnd);
	fclose(input_fasta);
	fclose(input_ss);
	fclose(input_pdb);
	fclose(input_phipsi);
	fclose(blossum_file);
	return 0;
}


