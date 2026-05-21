#  Running Local Lazypipe on Puhti
### For VZRU members

This manual will guide you through installing and running a local version of Lazypipe v3.1 on *Puhti CSC HPC*

In this module you will learn to:

* set up working environment on *CSC Puhti*
* install local version of *Lazypipe* to your application directory
* analyze individual samples with *Lazypipe*
* analyze sample batches with *Lazypipe* array jobs

**Prerequisites:**

* account on CSC Puhti
* no experience with Unix command line or NGS analysis is required 

For more information please refer to these guides:

* [**Getting started with Puhti**](https://docs.csc.fi/support/tutorials/puhti_quick/)
* [**Lazypipe User Guides**](https://bitbucket.org/plyusnin/lazypipe/wiki)

### Table of Content
1. [Installing Lazypipe](#InstallingLazypipe)
    * [Login to Puhti Web-interface](#LoginPuhti)
    * [Setting up directories](#SettingUpDirectories)
    * [Cloning repository](#CloningRepository)
    * [Installing reference databases](#InstallingReferenceDatabases)
2. [Running Lazypipe](#RunningLazypipe)
    * [Example 1: analysing sample data on command line](#Example1)
    * [Example 2: analysing sample data using sbatch script](#Example2)
    * [Example 3: analyse batches of fastq files with array-jobs](#Example3)
3. [Citing Lazypipe](#Citing)
4. [Contact](#Contact)

<a id="InstallingLazypipe"></a>
## Installing

<a id="LoginPuhti"></a>
### Login to CSC Puhti server

Both MacOS and Windows users can access Puhti via Puhti web-interface. We recommend this option for all users that are new to Unix/CSC working environment:

* Login to Puhti web-interface by following the link: [**Puhti web-interface**](https://www.puhti.csc.fi/)
* From the main Dashboard click on "Login node shell" to open the terminal

<a id="SettingUpDirectories"></a>
### Setting up directories

After you have logged in to Puhti continue working in your terminal.

Start by listing projects you have access to:

```
csc-workspaces
```

This tutorial assumes you have access to project *COVID19 Analysis of SARS-CoV-2 sequences* (*project_2002989*). In this case all *reference databases* and *background filters* are pre-installed and you only need to create directory for *Lazypipe data*. You can use any directory on *scratch*, for example:

    mkdir -p /scratch/project_2002989/$USER/lazypipe
    
Now add environment variables that will point to the location of *Lazypipe* *databases*, *NCBI taxonomy*, *Background filters* (i.e. *Host genomes*), your *data directory*, and *Lazypipe application directgory*. From Puhti Dashboard open *Home Directory*. Select *Show Dotfiles*, then open for editing file named `.bashrc` (configuration files for the unix terminal). In the `.bashrc` add the following lines:

    export databases=/scratch/project_2002989/lazydata/databases
    export hostgenomes=/scratch/project_2002989/lazydata/hostgen
    export taxonomy=/scratch/project_2002989/lazydata/taxonomy
    export mydata=/scratch/project_2002989/$USER/lazypipe
    export lazypipe=/projappl/project_2002989/$USER/lazypipe
    
In the terminal load `.bashrc` settings and create `$mydata` and `$lazypipe` directories:

    source $HOME/.bashrc
    mkdir -p $mydata $lazypipe
    
<a id="CloningRepository"></a>
### Cloning the repository

Clone *Lazypipe* repository to `$lazypipe` directory you created in the previous step:

    git clone https://plyusnin@bitbucket.org/plyusnin/lazypipe.git  $lazypipe

Load module dependencies and check that *Lazypipe* is installed by printing online help:

    cd $lazypipe
    module load r-env-singularity biokit lazypipe pigz
    perl lazypipe.pl -h

<a id="InstallingReferenceDatabases"></a>
### Installing Reference Databases

If you are working under project *COVID19 Analysis of SARS-CoV-2 sequences* (*project_2002989*) all required databases are pre-installed. If you are working under different project or need to install additional databases please see https://bitbucket.org/plyusnin/lazypipe/wiki/UserGuide.v3.1#InstallingReferenceDatabases. 


<a id="RunningLazypipe"></a>
## Running Lazypipe

<a id="Example1"></a>
### Example 1: analysing sample data on command line

Here we will run *Virus Discovery* pipelines with sample data (Mink feces, Illumina PE):

Run the main pipeline that will *preprocess* reads, *remove human and mink background*, *assemble*, *realign reads* to assembly, *annotate* with Minimap2 > BLASTN > BLASTP on *RefSeq viruses*, *print reports* and pack results to a \*.tar.gz* archive. Sample library is small so you can run all steps on Puhti *login node*:

    cd $lazypipe
    perl lazypipe.pl -1 data/samples/M15small_R1.fastq -p main --hostgen Neovison_vison,Homo_sapiens --ann1 minimap.refseq.vi,blastn.refseq.vi,blastp.uniref100.vi --res $mydata -s M15test -v
    
Note that we used background filters called `Neovison_vison` and `Homo_sapiens`. Available background filters are listed by name in `config.yaml` in section `host.databases`. Also note which reference databases were used: `minimap.refseq.vi`, `blastn.refseq.vi` and `blastp.uniref100.vi`. Available databases are listed in `config.yaml` in section `ann.databases`.
    
Once finished check your results:

    cd $mydata/M15test
    ls -l
    
    
<a id="Example2"></a> 
### Example 2: analysing sample data using a sbatch job script

In this example we will create a sample *sbatch job* that will run *Virus-Discovery-Chain1* on sample data.

Copy template *sbatch job script* to your working directory:

    cd $lazypipe
    cp /scratch/project_2002989/lazydata/templates/template.job.bash $lazypipe/example.job.bash
    
Then from Puhti Dashboard navigate to your `$lazypipe` directory and open `example.job.bash` for editing.

Set names for main logs:

    #SBATCH --error=logs/job.name.log
    #SBATCH --output=logs/job.name.log

Add paths for forward and reverse reads. If your reads are named following `*_R1.fastq` `*_R2.fastq` convention, the pipe will guess reverse reads from forward reads.

    r1=data/samples/M15small_R1.fastq
    r2=data/samples/M15small_R2.fastq
    
Add path for results and step-wise logs or keep default.

    res=$mydata/results
    logs=$res/logs

Then save your file and submit to the *sbatch queue*

    sbatch example.job.bash
    
Check the state of your job

    sacct
    
If errors occure, check the error log `logs/job.name.log`.

When completed check results from `$mydata/results/M15small`.


<a id="Example3"></a>
### Example 3: analyse batches of fastq files with array-jobs

Copy template *array job script* to your working directory:

    cd $lazypipe
    cp /scratch/project_2002989/lazydata/templates/template.array.job.bash $lazypipe/example.ajob.bash
    
Create a list of input files, one file per line. Here, we will list all FASTQ files in sample-directory, but you can list all FASTQ files in your input data directory. We will assume that your filenames follow *_R1.fastq*/*_R2.fastq* convention, so we will list only forward reads. You can optionally add a *sample name* for each input FASTQ library that will be used to name output directories.

     ls -1 $lazypipe/data/samples/*_R1.fastq.gz 1> array.files
    
Then from Puhti Dashboard navigate to your `$lazypipe` directory and open `example.ajob.bash` for editing.

Set array indices to match the number of files you have in your `array.files`. For example if you have 10 files set this to:

    #SBATCH --array=1-10

Add names for the main error logs. The `%a` placeholder stands for the index of the array-job that is executed. For example log for job-1 will be printed to `logs/job.name.1.log`.

    #SBATCH --error=logs/job.name.%a.log
    #SBATCH --output=logs/job.name.%a.log

Add paths for forward and reverse reads. If your reads are named following `*_R1.fastq` `*_R2.fastq` convention, the pipe will guess reverse reads from forward reads.

    r1=$(sed -n ${SLURM_ARRAY_TASK_ID}p array.files | cut -f1)

If you have samples names in the second column of `array.files` uncommend `sample=..` line and add `-s $sample` to the pipeline call:

    sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p array.files | cut -f2)
    perl lazypipe.pl .. -s $sample 

Add path for results and step-wise logs or keep default.

    res=$mydata/results
    logs=$res/logs

Then save your file and submit to the *sbatch queue*

    sbatch example.ajob.bash
    
Check the state of your job

    sacct
    
If errors occure, check the error log `logs/job.name.1.log`.

When completed check results from `$mydata/results/mysample`.


#### End notes

For more information see [**Lazypipe User Guides**](https://bitbucket.org/plyusnin/lazypipe/wiki)
