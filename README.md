----------------------------------------
HOW TO RUN FLIB-FLEX 
----------------------------------------
FlibFlex generates a fragment library for protein structure prediction using SAINT2, where part of
the structure is known. 

A Validator PDB structure containing the known region is provided or generated. For the
corresponding region of the fragment library, these torsion angles and secondary structure (DSSP) are used
for fragment generation, and homologues are not excluded. Fragments for this region are validated
against the structure, and fragments with higher RMSD than the "Flex" cutoff are discared (a minumum
of 20 fragments per position are included). 


---------------------------------------
FILE PREPARATION
--------------------------------------

```bash
export FLIB=/data/icarus/west/Flib-Flex/
export SCRIPTS=~/Project/Scripts/
```

Note that PDB files need to be cleaned up before use. This means indexed from 1, single chain (A),
standard amino acids and with multiple occupancies removed. Some scripts that can help with this:

```bash
$SCRIPTS/rechain.py $TARGET.pdb
for PDB in $(cat list.txt); do mv $PDB.pdb orig_$PDB.pdb; $SCRIPTS/pdb-tools/pdb_delocc.py orig_$PDB.pdb > $PDB.pdb; done
```

To generate the validator where a terminal end is missing:

```bash
bash $FLIB/scripts/generate_validator.sh TARGET VALIDATORID TERMINUS SEGMENT
``` 

Where segment is the length of the known region.

For each target, you will need:

$TARGET.fasta.txt     : fasta sequence of target
$TARGET.fasta
$TARGET.ss8           : Predicted secondary structure (the 8-category output of DeepCNF)
$TARGET.spd3          : Predicted torsion angles (the output of Spider3)

To convert the DeepCNF output to a form readable by Flib, use the following script:

```bash
for TARGET in $(cat lisst.txt); do $FLIB/scripts/convert_ss8.sh $TARGET ; done
```

This will remove the header and summarise the 8 columns into 3 columns. The actual values are currently
not used in Flib, only the category for each residue. 

$TARGET.fasta.ss      : Predicted secondary structure in correct format

-----------------------------------------
PREDICTED CONTACTS
-----------------------------------------
These will need to be prepared separately.

$TARGET.fasta         : fasta sequence of target

```bash
./run_metapsicov $TARGET.fasta >> $TARGET.log
(or use parallelise_contacts.sh)
```

This will generate, among other intermediate files:
$TARGET.metapsicov.stage1
$TARGET.metapsicov.stage2

For protein structure prediction we use the output of stage1.
We next need to validate the contacts against the known region of the protein, removing known false
positives and keeping all predictions from the missing regions.

First use:

```bash
for PDB in $(cat list.txt); do bash /data/icarus/west/Flib-Flex/scripts/validate_contacts.sh $PDB model_$PDB ; done
```

This uses scripts in QualityAssessment to align the target fasta (used to generate the predictions)
and the structure (used to validate the predictions) to validate the prediction. This will output:

$TARGET.metapsicov_stage1             : the prediction file with 1 or 0 for each contact
$TARGET.metapsicov_stage1_messages    : alignment information, contact count, TP/FP rates
$TARGET.log                           : BioPython messages

```bash
python $FLIB/scripts/get_vcon.py $TARGET $MISSING $TERMINUS
```

Note that model_$TARGET.fasta must exist ( cp $PDB.fasta model_$TARGET.fasta)

This parses the predicted contacts, removing all known incorrect contacts and keeping all
predictions from the missing region. This outputs:

$TARGET.vcon         : validated predicted contact list

---------------------
RUNNING FLIB-FLEX
---------------------
Use parallelise_pipeline.sh to run 

```bash
$FLIB/run_flib_flex_pipeline.sh $TARGET -nohoms -true_ss -generate_hhr -generate_validator
-validator_pdb $VALIDATOR -chain A -segment $SEG -terminus $TERMINUS -flex 1.5
```

- generate_hhr        : run HHblits to get homologues
- nohoms              : Removes >=95% sequence identical homologues from the Missing region
- generate_validator  : generate the validator PDB structure based on -segment and -terminus options
- validator_pdb       : ID of structure used for validator generation
- true_ss             : Uses the DSSP from the validator structure 
- flex                : Angstrom cutoff to remove bad validated fragments from the library
 







  



