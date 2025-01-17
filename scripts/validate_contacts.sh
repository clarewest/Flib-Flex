TARGET=$1;
VALIDATOR=$2
CHAIN_TARGET=A;

export QA=/data/icarus/not-backed-up/west/QualityAssessment/

# Clare edit: _stage1 will be empty if there are no contacts even if this has already run
if [ ! -s $TARGET.metapsicov_stage1_messages ]
then
    if [ ! -s $TARGET.metapsicov.stage1 ]
    then
        echo "Contact file for $TARGET not found. Moving onto the next target" 1>&2;
        exit;
    fi
    # Compute the true contact map and the PDB sequence:
    $QA/bin/getcontactsSparse3.py $VALIDATOR.pdb $CHAIN_TARGET 8.0 2> $TARGET.log
    # Contacts are usually calculated using the fasta sequence of the target ($TARGET.metapsicov.stage1).
    # This line corrects the contacts so the numbering matches the sequence on the PDB (and not on the fasta)
    $QA/bin/convert $VALIDATOR.proxy_fasta $TARGET.aln $TARGET.metapsicov.stage1 $VALIDATOR.proxy_map > $TARGET.metapsicov_stage1 2> $TARGET.metapsicov_stage1_messages
fi
