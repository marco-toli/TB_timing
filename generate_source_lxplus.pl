#$baseDir = "/afs/cern.ch/work/m/mlucchin/TB_timing/jobs";
$baseDir = ".";
$jobDir = $baseDir;


# create job file for each configuration set

$NJOBS = 10;

for($iJob = 0; $iJob < $NJOBS; ++$iJob)
{
    $jobFile = $jobDir."/job_".$iJob.".sh\n" ;
    system("cat ".$baseDir."/job_template.sh | sed -e s%ID%".$iJob."%g > ".$jobFile);    
}



# 
$lanciaFile = "./lancia_jobs.sh";
open(LANCIAFILE, ">", $lanciaFile);

for($iJob = 0; $iJob < $NJOBS; ++$iJob)
{
	print LANCIAFILE "bsub -cwd /afs/cern.ch/work/m/mlucchin/TB_timing/jobs/ -q 1nh /afs/cern.ch/work/m/mlucchin/TB_timing/jobs/job_".$iJob.".sh \n"; 
}
