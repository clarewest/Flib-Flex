## Usage: generate_validator.sh TARGET VALIDATORID TERMINUS SEGMENT
SCRIPTS=~/Project/Scripts/

OUTPUT=$1
VALIDATOR=validator_$OUTPUT
VALIDATORPDB=$2
terminus=$3
segment=$4             # Length of known region
CHAIN=A
FLAG=$5

if [ "$FLAG" == "template" ];
then VALIDATOR=template_$OUTPUT
fi

### Setting start and end residues of validator (inclusive)###
LENGTH=$(tail -n1 $OUTPUT.fasta.txt | tr -d '\n' | wc -c)
if [ $terminus = "C" ]; then
  begin=1
  end=$segment
elif [ $terminus = "N" ]; then
  begin=$((LENGTH-segment+1))
  end=$LENGTH
else
  echo "Not a valid terminus"
  exit
fi

$SCRIPTS/get_domain_from_pdb.py $VALIDATORPDB.pdb $CHAIN $begin $end $VALIDATOR.pdb
