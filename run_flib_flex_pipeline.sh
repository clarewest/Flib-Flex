#!/bin/bash

#set -u
#set -e

# SET UP THE PATHS FOR REQUIRED DEPENDENCIES:
export PSIPRED=/data/icarus/west/psipred/
export SPIDER=/data/icarus/west/SPIDER2/misc/
export HHSUITE=/data/icarus/west/hh-suite/bin/
export HHBLITSDB=/data/icarus/west/Databases/pdb70
export BLAST=~/ncbi-blast-2.7.1+/bin
export BLASTDB=/data/icarus/west/Databases/
export BLAST_PDB=/data/icarus/west/Databases/pdbaa
export FLIB=/data/icarus/west/Flib-Flex/
export PDB=/data/icarus/west/PDB/

# ==========================================================================
# DO NOT CHANGE ANYTHING BELOW HERE
# ==========================================================================

usage="\n./$(basename "$0") [options] <NAME>\n\n
Flib: A Fragment Library Generation Software\n\n
Positional arguments:\n
\t<NAME>\t\t\tprefix of the input files\n\n
Optional arguments:\n
\t-h | --help\t\tshow this help text\n
\t-generate_ss\t\trun a local version of PSIPRED\n
\t-generate_spider2\trun a local version of SPIDER2\n
\t-generate_hhr\t\trun a local version of HHblits\n
\t-disable_flib\t\tdon't run a local version of Flib\n
\t-noparsing\t\tdon't parse Flib libraries to SAINT2 format\n
\t-true_ss\t\tspecify to use the native SS (according to DSSP)\n
\t-nohoms\t\tspecify to omit homologs from search\n\n
\t-disable_flex\t\tdon't process lib3000\n\n
\t-flex\t\tspecify to filter out fragments more than FLEX A from validator\n\n
\t-validator_pdb\t\tpdb code of structure for validating against\n\n
\t-segment\t\tlength of region to allow homologues and validate\n\n
\t-terminus\t\tterminus of unknown region\n\n"

# DEFAULT VALUES FOR COMMAND LINE ARGUMENTS
generate_ss=false		# running local version of PSIPRED.
generate_spider2=false		# running local version of SPINE-X.
generate_hhr=false		# running local version of HHBlits.
generate_flib=true		# running local version of Flib.
disable_flex=false
parse_flib=true	        	# parsing Flib libraries to SAINT2 format.
remove_homologs=false		# removing homologs from frag. libraries.
flex=false              # special behaviour for flex region
segment=false         # length for which homologues are allowed 
terminus=false        # terminus of prediction
enrich=false

# CONSUME COMMAND LINE ARGUMENTS
found_positional=false
while [[ $# -gt 0 ]];
do
  key="$1"

  # Check all flags -- basic method but assume positional arguments are no more than one
  if [[ "$key" == -* ]];
  then
    case $key in
      -h|--help)
        echo -e $usage
        exit
        ;;
      -generate_ss)
        generate_ss=true
        ;;
      -generate_spider2)
        generate_spider2=true
        ;;
      -generate_hhr)
        generate_hhr=true
        ;;
      -disable_flib)
        generate_flib=false
        ;;
      -noparsing)
        parse_flib=false
        ;;
      -nohoms)
        remove_homologs=true
        ;;
      -true_ss)
        true_ss=true
        ;;
      -begin)
        begin=$2
        shift
        ;;
      -end)
      end=$2;
      shift
      ;;
    -offset)
      offset=$2;
      shift
      ;;
    -chain)
      CHAIN="$2";
      shift
      ;;
    -flex)
      FLEX=$2;
      shift
      ;;
    -terminus)
      terminus=$2;
      shift
      ;;
    -segment)
      segment=$2;
      shift
      ;;
    -validator_pdb)
      VALIDATORPDB=$2
      shift
      ;;
    -*|--*)
      echo "Unknown option" $key
      exit
      ;;
  esac
elif [[ $found_positional = true ]];
then
  echo "More than one positional argument provided"
  exit
else
  OUTPUT=$key
  found_positional=true
fi

shift
done

# CHECK THAT OUTPUT IS DEFINED
if [[ -z $OUTPUT ]];
then
  echo -e $usage
  exit
fi

##### GENERAL SET UP #####

# Generate SS Prediction using PSIPRED
if [ "$generate_ss" = true ] ; then
  echo "--------------------------------------"
  echo "Generating Secondary Structure Prediction using PSIPRED:"
  $PSIPRED/runpsipredplus ./$OUTPUT.fasta.txt
  echo "Done"
  echo "--------------------------------------"
fi

# Generate SPIDER2 Torsion Angle Prediction
if [ "$generate_spider2" = true ] ; then
  echo "Generating SPIDER2 Torsion Angle Prediction:"
  cat $OUTPUT.fasta.txt > $OUTPUT.seq
  $SPIDER/run_local.sh $OUTPUT.seq 2> $OUTPUT.spd_err
  rm $OUTPUT.seq
  echo "Done"
  echo "--------------------------------------"
fi

# Run HHSearch
if [ "$generate_hhr" = true ] ; then
  echo "Generating HHSearch File:"
  $HHSUITE/hhblits -d $HHBLITSDB -i ./$OUTPUT.fasta.txt -o $OUTPUT.hhr
  echo "Done"
  echo "--------------------------------------"
fi

# Generate list of Homologs ##### FOR THIS SET ONLY HOMOLOGUES WITH SEQ ID >= 95% ARE REMOVED 
if [ "$remove_homologs" = true ] ; then
  echo "Generate list of homologs:"
  $BLAST/blastp -query ./$OUTPUT.fasta.txt -db $BLAST_PDB -evalue 0.05 -outfmt 6  > $OUTPUT.blast
  ### this used to be $4
  awk '{if ($3>=95) print $2}' $OUTPUT.blast | sed -e "s/.*pdb|//g" | cut -c 1-4 | sort | uniq | tr '[:upper:]' '[:lower:]' > $OUTPUT.homol
  echo "Done"
  echo "--------------------------------------"
fi

##### FLIB #####
if [ "$generate_flib" = true ] ; then
  echo "Generating FLIB File:"
  if [ "$disable_flex" = false ] ; then
    LENGTH=$(tail -n1 $OUTPUT.fasta.txt | tr -d '\n' | wc -c)
    if [ $terminus = "C" ]; then
      begin=1
    end=$segment
  elif [ $terminus = "N" ]; then
    begin=$((LENGTH-segment))
  end=$LENGTH
else
  echo "Not a valid terminus"
  exit
fi
echo "Using $VALIDATORPDB as validation from $begin to $end"
    fi
    if [[ -z $begin ]] ; then
      begin=1;
    fi

    if [[ -z $end ]] ; then
    end=$(wc $OUTPUT.fasta.ss | awk '{print $1}');
  fi

  if [[ -z $offset ]] ; then
    offset=0;
  fi

  if [ "$true_ss" = true ] ; then
    echo "Using real secondary structure from $begin to $end"
    $FLIB/bin/Flib -i $OUTPUT -C $OUTPUT.con --true_ss -b $begin -e $end -O $offset -N $CHAIN 2> $OUTPUT.log > $OUTPUT.lib3000
  else
    $FLIB/bin/Flib -i $OUTPUT -C $OUTPUT.con 2> $OUTPUT.log > $OUTPUT.lib3000
  fi

  echo "lib3000 generated:"    

  sort -k 10,10n -k 13,13n $OUTPUT.lib3000 > $OUTPUT.tmp;                   # Sorts LIB3000
  mv $OUTPUT.tmp $OUTPUT.lib3000;                                           
  cp $OUTPUT.lib3000 $OUTPUT.lib3000_ori
fi

if  [ "$disable_flex" = false ] ; then

  ### Validate fragments ###
  echo "Validating fragments:"
  awk -v START=$begin -v STOP=$end '{if (($6>=START)&&($6<=STOP)) print $0}' $VALIDATORPDB.pdb > validator_$OUTPUT.pdb
  cp $OUTPUT.fasta.txt validator_$OUTPUT.fasta.txt
  $FLIB/bin/LibValidator validator_$OUTPUT $OUTPUT.lib3000 $CHAIN > $OUTPUT.lib3000_rmsd 2> $OUTPUT.validator_log   # rmsd of fragments against validator pdb
  cp $OUTPUT.lib3000_rmsd $OUTPUT.lib3000

  ### HOMOLOG REMOVAL only from missing region ####
  echo "Removing homologues"
  if [ "$remove_homologs" = true ] ; then
    cp $OUTPUT.lib3000 $OUTPUT.lib3000_nh;
    for HOMOLOG in $(cat $OUTPUT.homol)
    do
      #			sed -e "/$HOMOLOG/d" "$OUTPUT".lib3000_nh > "$OUTPUT".tmp
      awk -v homolog=$HOMOLOG  '{if (($1!=homolog)||($14>-1.0)) print $0}' "$OUTPUT".lib3000_nh > "$OUTPUT".tmp
      mv $OUTPUT.tmp $OUTPUT.lib3000_nh
    done
    mv $OUTPUT.lib3000_nh $OUTPUT.lib3000
  fi


  ###### lib3000 has been generated and validated to lib3000_rmsd #####
  #       sort -k 10,10n -k 14,14nr $OUTPUT.lib3000_rmsd > $OUTPUT.lib3000_rmsd_a                  # Sort by rmsd (descending)
  #       sort -k 10,10n -k 13,13n  $OUTPUT.lib3000_rmsd > $OUTPUT.lib3000_rmsd_b                  # Sort by predicted torsion angle score
  #       sort -k 10,10n -k 14,14n -k 13,13n $OUTPUT.lib3000_rmsd > $OUTPUT.lib3000_rmsd_c    # Sort by rmsd and predicted torsion angle score     
  #       $FLIB/bin/filter_flex2_lib2 $OUTPUT $OUTPUT.lib3000_rmsd $FLEX $segment		            # Filters out fragments above RMSD threshold and/or sorts by torsion angle score
  #       sed 's/.....$//' $OUTPUT.lib3000_flex$FLEX > $OUTPUT.lib3000		# For now, remove RMSD column and rename in to continue using original scripts
  echo "Processing lib3000"
  python $FLIB/scripts/process_validated.py $OUTPUT.lib3000
  mv $OUTPUT.lib3000_proc $OUTPUT.lib3000
  sort -k 10,10n -k 13,13n $OUTPUT.lib3000 > $OUTPUT.tmp 
  mv $OUTPUT.tmp $OUTPUT.lib3000
  $FLIB/bin/filterlib2 $OUTPUT $OUTPUT.lib3000                             # This will generate LIB20 and LIB500
  mv $OUTPUT.lib20 $OUTPUT.lib20_ori
  mv $OUTPUT.lib500 $OUTPUT.lib500_ori



  if [ "$enrich" = true ]; then  
    ### Fragment Library Enrichment ###
    $FLIB/bin/Flib_Enrich $OUTPUT.lib500_ori $PDB 0.5 0 > $OUTPUT.clib 2> $OUTPUT.error # This will enrich LIB20 with fragments from LIB500: CLIB
    cat $OUTPUT.lib20_ori >> $OUTPUT.lib_tmp2                                     # Merges LIB20 and CLIB
    cat $OUTPUT.clib >> $OUTPUT.lib_tmp2                                      # ...
    sort -k 10,10n $OUTPUT.lib_tmp2 > $OUTPUT.lib_final                       # Sorts LIB20+CLIB : LIB_FINAL
    rm $OUTPUT.lib_tmp2     
  else
    ### could add ~6.5 fragments without RMSD here ###
    cp $OUPUT.lib20_ori $OUTPUT.lib_final 
  fi

  ### Parsing fragments from Threading hits: ###
  python $FLIB/scripts/parse_hhr.py $OUTPUT > $OUTPUT.lib_hhr 2> $OUTPUT.log    # Creates lib from threading hits: LIB_HHR
  sort -k 10,10n -k 11,11nr $OUTPUT.lib_hhr > $OUTPUT.ordered;                    # Sorts LIB_HHR
  $FLIB/bin/parse_hhr $OUTPUT.ordered > $OUTPUT.lib9               # Parses LIB_HHR
  cut -f 1-12 $OUTPUT.lib9 > $OUTPUT.lib9_aux
  mv $OUTPUT.lib9_aux $OUTPUT.lib9
  awk 'BEGIN {OFS = "\t"} ; {print $0,"0.0"}' $OUTPUT.lib9 > $OUTPUT.9_lib
  $FLIB/bin/LibValidator validator_$OUTPUT $OUTPUT.9_lib $CHAIN > $OUTPUT.rmsd_9_lib 2>> $OUTPUT.validator_log   # rmsd of fragments against validator pdb

  ### HOMOLOG REMOVAL ####
  if [ "$remove_homologs" = true ] ; then
    cp $OUTPUT.rmsd_9_lib $OUTPUT.rmsd_9_lib_nh;
    for HOMOLOG in $(cat $OUTPUT.homol)
    do
      #              awk -v homolog=$HOMOLOG -v  '{if (($1!=homolog)||($14>-1.0)) print $0}' "$OUTPUT".lib3000_nh > "$OUTPUT".tmp
      awk -v homolog=$HOMOLOG -v start=$begin -v stop=$end '{if (($1!=homolog)|| ( (($10+$4) < start ) || (10 >= stop))) print $0}' "$OUTPUT".rmsd_9_lib_nh > "$OUTPUT".tmp
      #	       	        sed -e "/$HOMOLOG/d" "$OUTPUT".rmsd_9_lib_nh > "$OUTPUT".tmp
      mv $OUTPUT.tmp $OUTPUT.rmsd_9_lib_nh
    done
    mv $OUTPUT.rmsd_9_lib_nh $OUTPUT.rmsd_9_lib
  fi
  $FLIB/bin/filterlib2 $OUTPUT $OUTPUT.rmsd_9_lib                 # Filter LIB_HHR so it does not contain more than 20 frags. per position.
  cat $OUTPUT.lib_final > $OUTPUT.combined
  awk '$6 == "O" {print }' $OUTPUT.lib20	   >> $OUTPUT.combined
  awk '$6 == "L" {print }' $OUTPUT.lib20     >> $OUTPUT.combined

  sort -k 10,10n -k 13,13nr $OUTPUT.combined > $OUTPUT.lib                        # The resulting library is LIB
  echo "Done"
  echo "--------------------------------------"
fi


##### LIBRARY PARSING #####

# Parse FLIB library into SAINT2 compliant format.
if [ "$parse_flib" = true ] ; then
  echo "Parsing Flib File:"
  python $FLIB/scripts/process_flex_new.py $PDB < $OUTPUT.lib > $OUTPUT.vlib 2> $OUTPUT.log
  mv $OUTPUT.vlib $OUTPUT.vlib_flex"$FLEX"_seg$segment.$terminus
  echo "Done"
  echo "--------------------------------------"
fi



