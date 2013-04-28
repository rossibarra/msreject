 
#$-S /bin/bash
#$-cwd
#$-N palu.island
#$-t 1-15
#$ -v LD_LIBRARY_PATH=/opt/lib:/usr/lib:/usr/local/lib:/opt/gridengine/lib/:/opt/gridengine/lib/lx26-x86/:/usr/include

PATH=/opt/gridengine/bin/lx26-x86:/opt/gridengine/bin/lx26-x86:/usr/java/jdk1.5.0_07/bin:/opt/globus/bin:/opt/globus/sbin:/opt/gridengine/bin/lx26-x86:/opt/gridengine/bin/lx26-x86:/usr/java/jdk1.5.0_07/bin:/opt/gridengine/bin/lx26-x86:/usr/kerberos/bin:/opt/gridengine/bin/lx26-x86:/usr/java/jdk1.5.0_07/bin:/usr/local/bin:/usr/bin:/bin:/usr/X11R6/bin:/opt/Bio/ncbi/bin:/opt/Bio/mpiblast/bin/:/opt/Bio/hmmer/bin:/opt/Bio/Emboss/bin:/opt/Bio/clustalw/bin:/opt/Bio/t_coffee/bin:/opt/Bio/phylip/exe:/opt/Bio/mrbayes:/opt/Bio/fasta:/opt/Bio/glimmer/bin://opt/Bio/glimmer/scripts:/opt/Bio/gromacs/bin:/opt/chromium/bin/Linux:/opt/eclipse:/opt/ganglia/bin:/opt/ganglia/sbin:/opt/lam/gnu/bin:/opt/maven/bin:/opt/openmpi/bin/:/usr/share/pvm3/lib:/usr/share/pvm3/lib/LINUX:/usr/share/pvm3/bin/LINUX:/opt/rocks/bin:/opt/rocks/sbin:/titus/users/jross/.sage/bin:/opt/Bio/ncbi/bin:/opt/Bio/mpiblast/bin/:/opt/Bio/hmmer/bin:/opt/Bio/Emboss/bin:/opt/Bio/clustalw/bin:/opt/Bio/t_coffee/bin:/opt/Bio/phylip/exe:/opt/Bio/mrbayes:/opt/Bio/fasta:/opt/Bio/glimmer/bin://opt/Bio/glimmer/scripts:/opt/Bio/gromacs/bin:/opt/chromium/bin/Linux:/opt/eclipse:/opt/ganglia/bin:/opt/ganglia/sbin:/opt/lam/gnu/bin:/opt/maven/bin:/opt/openmpi/bin/:/usr/share/pvm3/lib:/usr/share/pvm3/lib/LINUX:/usr/share/pvm3/bin/LINUX:/opt/rocks/bin:/opt/rocks/sbin:/titus/users/jross/.sage/bin:/opt/Bio/ncbi/bin:/opt/Bio/mpiblast/bin/:/opt/Bio/hmmer/bin:/opt/Bio/Emboss/bin:/opt/Bio/clustalw/bin:/opt/Bio/t_coffee/bin:/opt/Bio/phylip/exe:/opt/Bio/mrbayes:/opt/Bio/fasta:/opt/Bio/glimmer/bin://opt/Bio/glimmer/scripts:/opt/Bio/gromacs/bin:/opt/chromium/bin/Linux:/opt/eclipse:/opt/ganglia/bin:/opt/ganglia/sbin:/opt/lam/gnu/bin:/opt/maven/bin:/opt/openmpi/bin/:/usr/share/pvm3/lib:/usr/share/pvm3/lib/LINUX:/usr/share/pvm3/bin/LINUX:/opt/rocks/bin:/opt/rocks/sbin:/titus/users/jross/.sage/bin:/titus/users/jross/bin
LD_LIBRARY_PATH=/opt/gridengine/lib/lx26-x86:/titus/users/jross/.sage/lib:/opt/gridengine/lib/lx26-x86:/opt/globus/lib:/opt/gridengine/lib/lx26-x86:/titus/users/jross/.sage/lib:/opt/gridengine/lib/lx26-x86:/opt/gridengine/lib/lx26-x86:/titus/users/jross/.sage/lib:/opt/gridengine/lib/lx26-x86:/opt/lam/gnu/lib:/opt/lam/gnu/lib:/opt/lam/gnu/lib:/usr/local/lib
scratch=/state/partition1

loci=26
taxon=palu
model=island
cutoff=0.30
nruns_multilocus=2000
singletons=yes 
verbose=-V
errorcheck=no
meanonly=
acceptall=
nmlog=
taulog=
n2log=
nalog=

#get model command
k=$SGE_TASK_ID
if [ $model == "iso" ]; then
		model_command="-ma x tbs tbs x -n 2 tbs -ej tbs 2 1 -eN tbs tbs" 
		modelnum=3
	#	echo "meantheta	meanrho	meannm1	meannm2	X	X	X	X	X	X	n2	tau	X	na" >  $scratch/$taxon.$model.$k; 
	elif [ $model == "island" ]; then
		model_command="-ma x tbs tbs x -n 2 tbs -ej tbs 2 1 -eN tbs tbs"
		modelnum=2
	#	echo "meantheta	meanrho	meannm1	meannm2	X	X	X	X	X	X	nm1	nm2	n2	tau	X	na" >  $scratch/$taxon.$model.$k; 
	elif [ $model == "sym" ]; then
		model_command="-ma x 0 0 x -n 2 tbs -ema tbs 2 x tbs tbs x -ej tbs 2 1 -eN tbs tbs" 
		modelnum=4
	#	echo "meantheta	meanrho	meannm1	meannm2	X	X	X	X	X	X	n2	tausmall	nm1	nm2	tau	X	na" >  $scratch/$taxon.$model.$k; 
	elif [ $model == "allo" ]; then
		model_command="-ma x tbs tbs x -n 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs"
		modelnum=1
	#	echo "meantheta	meanrho	meannm1	meannm2	X	X	X	X	X	X	nm1	nm2	n2	tausmall	tau	X	na" >  $scratch/$taxon.$model.$k; 
	else echo "model not recognized"; kill 0;
fi

#funxions
function checkfile(){
	if [ ! -s "$1" ]; then
		echo "FATAL ERROR: $1 does not exist or is empty"
		kill 0; 
	fi
}

# checkfiles and mv files to scratch
checkfile $taxon".in"
cp $taxon".in" $scratch/$taxon.$model".in."$k
echo "model = $model" >> $scratch/$taxon.$model".in."$k
echo "bpfile = $scratch/$taxon.$model.bpfile.$k" >> $scratch/$taxon.$model".in."$k
echo "nreps =  $nruns_multilocus" >> $scratch/$taxon.$model".in."$k
prior_in=$taxon.$model".in."$k
checkfile $taxon.bpfile
cp $taxon.bpfile $scratch/$taxon.$model.bpfile.$k
if [ $singletons == "yes" ] 
	then datafile="data.$taxon"; ffilter=""
	else datafile="data.$taxon.q"; ffilter="-F 2";
fi
checkfile $datafile
cp $datafile $scratch/$datafile.$k

#seed prior_in and generate prior
r0=$(($RANDOM*$k))
echo "seed = $r0" >> $scratch/$prior_in
echo "prior = $scratch/prior.$taxon.$model.$k" >> $scratch/$prior_in
priorgen $nmlog $taulog $n2log $nalog --arg-file $scratch/$prior_in
checkfile $scratch/prior.$taxon.$model.$k
	
#seed ms, generate samples, reject	
r1=$(($RANDOM*$k))
r2=$(($RANDOM*$k))
r3=$(($RANDOM*$k))
if [ $errorcheck == "yes" ]
then
	cat $scratch/$prior_in >> error
	cp  $scratch/prior.$taxon.$model.$k prior.error
	echo "msnsam tbs $(($nruns_multilocus*$loci)) -t tbs -r tbs tbs -I 2 tbs tbs $model_command $ffilter -seed $r1 $r2 $r3 < $scratch/prior.$taxon.$model.$k" >> error
	echo "rejection_zea -f $scratch/$taxon.$model.$k -d $scratch/$datafile.$k -l $loci -m $modelnum -c $cutoff $verbose $meanonly" >> error
fi
msnsam tbs $(($nruns_multilocus*$loci)) -t tbs -r tbs tbs -I 2 tbs tbs $model_command $ffilter -seed $r1 $r2 $r3 < $scratch/prior.$taxon.$model.$k |
rejection_zea -f $scratch/$taxon.$model.$k -d $scratch/$datafile.$k -l $loci -c $cutoff $verbose $meanonly $acceptall > $scratch/$taxon.$model.out.$k

#cat to main outfile remove extra files from run, go to next run
rm -f $scratch/$prior_in $scratch/prior.$taxon.$model.$k  $scratch/$datafile.$k $scratch/$taxon.$model.bpfile.$k 
mv $scratch/$taxon.$model.$k ~/projects/maize/abc/$taxon/$taxon.$model.$k
mv $scratch/$taxon.$model.out.$k ~/projects/maize/abc/$taxon/$taxon.$model.out.$k
