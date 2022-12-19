#! /bin/csh -f

## Preprocessing script for converting MPMs into surface recon-able
## volumes that can be fed through FS processing
## Adapted from Fred Dick's PreProcScriptMPM (for single subjects)
##
## NOTE: AFNI and freesurfer need to be on path
## NOTE: R1 and PD volume prefixes are hardcoded as are the truncation thresholds
##       and the scaling factor. Change if needed for your data.
##
## Usage
echo $#argv
if ($#argv < 2) then
	echo ""
	echo ""
	echo "###################################################################"
	echo ""
	echo " $0 Usage: DataDir SubjectList"
	echo ""
	echo "     Each argument is separated by a space"
	echo "     DataDir is the absolute path to hi-res MPMs"
	echo "     SubjectList is a list of subject names to recon (e.g., (sub-01 sub-02 sub-03))"
	echo ""
	echo "###################################################################"
	echo ""
	goto done
endif


## grab input args
set datadir=$1
shift argv
set subjlist=$argv

## check freesurfer setup
 if ! $?FREESURFER_HOME then
 	echo "You need to setup freesurfer environment and run again"
 	goto done
 endif

## check freesurfer dir exists!
 if (! -d $FREESURFER_HOME) then
 	echo "I can't find MGH freesurfer, please make sure it exists"
 	goto done
 endif

## check if SUBJECTS_DIR exists
 if ! -d $SUBJECTS_DIR then
 	echo "I can't find SUBJECTS_DIR!"
 	goto done
 endif

foreach sub ($argv)
## check if datadir exists
	echo "===Computing subject $sub"
 	if (! -d $datadir/$sub/) then
 	echo "$datadir does not exist, please make sure datadir is correct"
 	goto done
 endif

## check if R1 volume exists
cd $datadir/$sub/Results
set r1vol = sub*R1.nii
 	if (! -e $datadir/$sub/Results/$r1vol) then
 	echo "$r1vol doesn't exist, please check name"
 	goto done
 endif

## check if PD volume exists
set pdvol = sub*PD.nii
 	if (! -e $datadir/$sub/Results/$pdvol) then
 	echo "$pdvol doesn't exist, please check name"
 	goto done
 endif

## confirm name assignment
if ($#argv == 2) then
	echo "Subject is $sub, DataDir is $datadir"
	echo "R1-Vol is $r1vol, PD-Vol is $pdvol"
endif

## Move to data directory
cd $datadir/$sub/Results

## To create the appropriate PD for input to mri_synthesize
##  remove all of the negative values and scale

## Here, $pdvol is the quantitative PD calculated by the hMRI Toolbox

3dcalc -a $pdvol -prefix $sub-PD-vol-noneg-scale.nii -expr '(a*step(a))*100'

## R1 volume ($r1vol) is created by the hMRI toolbox,
## needs to be  T1 (e.g., 1/R1) and also truncated,
## with no negative values, and scaled

3dcalc -a ./$r1vol -prefix ./$sub-T1.nii -expr '1000/a'
3dcalc -a ./$sub-T1.nii -prefix ./$sub-T1-mask.nii -expr 'within(a,0,8000)'
3dcalc -a ./$sub-T1-mask.nii -b ./$sub-T1.nii -prefix ./$sub-T1-trunc.nii -expr 'a*b'

## Output of mri_synthesize help is:

## usage: mri_synthesize [options] <TR> <alpha (deg)>
## <TE> <T1 volume> <PD volume> <output volume>

## The -w switch will use a fixed weighting in order
## to generate an output volume with optimal gray/white contrast

## It seems to need those three parameters, but
## ignores them....

mri_synthesize -w 20 30 2.5 ./$sub-T1-trunc.nii \
./$sub-PD-vol-noneg-scale.nii ./$sub-synth.nii

## Finally eliminate extreme values and scale

## This should get signal values into correct ballpark (80 < wm < 120)
## NOTE - the scaling factor may have to be a variable
## I am setting it so that it is easy to find.
## NOTE - larger scale number means lower vals as is denominator

set scale = 3

3dcalc -a ./$sub-synth.nii -prefix ./$sub-synth-trunc-scale.nii.gz \
-expr '(a*(within(a,0,820))/'$scale')'

## Labels to goto to get out of script
done:
	exit 0
error:
	exit 1



