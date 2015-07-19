$baseDir = "/afs/cern.ch/work/m/mlucchin/TB_timing/jobs";
$jobDir = $baseDir;


# create job file for each configuration set

$NJOBS = 100;

for($iJob = 0; $iJob < $NJOBS; ++$iJob)
{
#    $jobDir = $baseDir."/jobdir_".$iJob;    
#    system("mkdir ".jobDir);
    $jobFile = $jobDir."/job_".$iJob.".sh\n" ;
    system("cat ".$baseDir."/job_template.sh | sed -e s%ID%".$iJob."%g > ".$jobFile);    
}



# 
$lanciaFile = "./lancia_jobs.sh";
open(LANCIAFILE, ">", $lanciaFile);

for($iJob = 0; $iJob < $NJOBS; ++$iJob)
{
#	print LANCIAFILE "mkdir ".$baseDir."/jobdir_".$iJob." \n";
	print LANCIAFILE "bsub -cwd /afs/cern.ch/work/m/mlucchin/TB_timing/jobs/output -q 2nd /afs/cern.ch/work/m/mlucchin/TB_timing/jobs/job_".$iJob.".sh \n"; 
}
