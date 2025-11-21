# Predicting Protein Structures with AlphaFold3 on the CHTC GPU Capacity

A two-phase workflow: alignment generation → structure prediction

AlphaFold3 (AF3) is a next-generation biomolecular structure prediction system capable of modeling multichain protein complexes, DNA/RNA interactions, ligands, and modified residues (e.g., 2′-O-methylated RNAs or glycosylated amino acids). This guide provides a step-by-step workflow for running AF3 on the CHTC HTC ecosystem, including how to structure your project, prepare input JSON files, run many inference jobs with HTCondor, and store your results using OSDF.

AlphaFold3 workloads in high-throughput environments are best organized into two separate job types:
* **Step 1: Generating the Alignments (Data-Only Pipeline)**
  * MMseqs2 search
  * HMMER search
  * Template retrieval (if enabled)
  * MSA + templates packaged into an AF3 “fold_input.json” 
* **Step 2: Predicting the Structure (Inference Pipeline)**
  * Load the precomputed fold_input.json 
  * Run the AF3 diffusion model 
  * Produce PDB models, ranking, trajectories, and metrics


You will learn how to:

* **Understand the overall computational architecture of AlphaFold3 on CHTC**, including how the data-generation and inference stages map to CPU and GPU resources. 
* **Design, organize, and manage large-scale AF3 workloads**, including preparing inputs, structuring job directories, and generating automated job manifests. 
* **Leverage CHTC’s GPU capacity for high-throughput structure prediction**, including selecting appropriate resources based on input complexity. 
* **Use containers, staged databases, and HTCondor data-transfer** mechanisms to build reproducible, portable, and scalable AF3 workflows. 
* **Submit and monitor hundreds to thousands of AF3 jobs**, using standard HTCondor patterns and best practices for reliable execution on distributed compute sites.

All of these steps run across hundreds (or thousands) of jobs using the HTCondor workload manager and Apptainer containers to execute your software reliably and reproducibly at scale. The tutorial uses realistic genomics data and emphasizes performance, reproducibility, and portability. You will work with real data and see how high-throughput computing (HTC) can accelerate your workflows.

![Workflow_Diagram.png](.images/Workflow_Diagram.png)

> [!NOTE]
> If you are new to running jobs on CHTC, complete the CHTC ["Roadmap to getting started
"](https://chtc.cs.wisc.edu/uw-research-computing/htc-roadmap/) and our ["Practice: Submit HTC Jobs using HTCondor"](https://chtc.cs.wisc.edu/uw-research-computing/htcondor-job-submission) guide before starting this tutorial.

**Let’s get started!**
- [Tutorial Setup](#tutorial-setup)
   * [Assumptions](#assumptions)
   * [Prerequisites](#prerequisites)
- [Understanding the AlphaFold3 Workflow](#understanding-the-alphafold3-workflow)
   * [The CPU-Only Pipeline: Generating Alignments (Step 1)](#the-cpu-only-pipeline-generating-alignments-step-1)
      + [What the Data Pipeline Does](#what-the-data-pipeline-does)
      + [Why This Runs on CPU Machines at CHTC](#why-this-runs-on-cpu-machines-at-chtc)
- [The GPU-Accelerated Pipeline: Structural Prediction (Step 2)](#the-gpu-accelerated-pipeline-structural-prediction-step-2)
      + [What the Inference Pipeline Does](#what-the-inference-pipeline-does)
      + [Why This Runs on GPU Machines at CHTC](#why-this-runs-on-gpu-machines-at-chtc)
- [Running AlphaFold3 on CHTC](#running-alphafold3-on-chtc)
   * [Set Up Your Software Environment](#set-up-your-software-environment)
   * [Data Wrangling and Preparing AlphaFold3 Inputs](#data-wrangling-and-preparing-alphafold3-inputs)
      + [Setting Up AlphaFold3 Input JSONs and Job Directories](#setting-up-alphafold3-input-jsons-and-job-directories)
      + [Preparing Your _List of (AlphaFold) Jobs_](#preparing-your-list-of-alphafold-jobs)
   * [Submit Your AlphaFold3 Jobs - CPU-Intensive Alignment Generation (Step 1)](#submit-your-alphafold3-jobs---cpu-intensive-alignment-generation-step-1)
      + [AlphaFold3 Databases Availability on CHTC](#alphafold3-databases-availability-on-chtc)
   * [Visualize Your AlphaFold3 Results](#visualize-your-alphafold3-results)
- [Next Steps](#next-steps)
   * [Software](#software)
   * [Data](#data)
      + [Key AF3 data components include:](#key-af3-data-components-include)
   * [GPUs](#gpus)
      + [Key AF3 GPU considerations:](#key-af3-gpu-considerations)
- [Getting Help](#getting-help)

<!-- TOC end -->


## Tutorial Setup

### Assumptions

This tutorial assumes that you:

* Have basic command-line experience (e.g., navigating directories, using bash, editing text files)
* Have a working CHTC account and can log into an Access Point (e.g., ap2001/2002.uw.osg-htc.org)
* Are familiar with HTCondor job submission, including writing simple `.sub` files and tracking job status with `condor_q`
* Have access to a machine with a GPU-enabled execution environment (provided automatically via CHTC's GPU Open Capacity)
* Have sufficient disk quota and file permissions in your CHTC `/home` and `/staging` directories

### Prerequisites

1. [ ] A CHTC HTC account. If you do not have one, request access at the [CHTC Account Request Page](https://chtc.cs.wisc.edu/uw-research-computing/form.html).
2. [ ] Basic familiarity with HTCondor job submission. If you are new to HTCondor, complete the CHTC ["Roadmap to getting started
"](https://chtc.cs.wisc.edu/uw-research-computing/htc-roadmap/) and read the ["Practice: Submit HTC Jobs using HTCondor"](https://chtc.cs.wisc.edu/uw-research-computing/htcondor-job-submission).
3. [ ] AlphaFold3 Model Weights. Request the AF3 model weights from the [DeepMind AlphaFold Team](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md#obtaining-model-parameters).

> [!WARNING]
> Requesting AlphaFold3 model weights requires agreeing to DeepMind's terms of service. Ensure you comply with all licensing and usage restrictions when using AF3 for research. This tutorial does not distribute AF3 model weights. Requesting the weights can take up to several weeks. Ensure you have them before starting the tutorial.

Log into your CHTC account:

    ```bash
    ssh user.name@ap##.uw.osg-htc.org
    ```

To obtain a copy of the tutorial files, you can:

* Clone the repository:

  ```bash
  git clone https://github.com/dmora127/tutorial-CHTC-AF3.git
  cd tutorial-CHTC-AF3/
  bash tutorial-setup.sh <username>
  ```
  _This script creates the directory structure in your home directory `/home/<user.name>/tutorial-CHTC-AF3/` and CHTC directory `/staging/<netID>/tutorial-CHTC-AF3/`, along with several subdirectories used in this tutorial._

* Or download the toy dataset using Pelican **CURRENT UNAVAILBLE - WORK IN PROGRESS**:

  ```bash
  pelican object get [PENDING PELICAN LINK] ./
  ```

## Understanding the AlphaFold3 Workflow

AlphaFold3 (AF3) uses a two-stage workflow that separates alignment generation from structure prediction. This separation allows AF3 to run efficiently in high-throughput environments like CHTC, where jobs may need to scale across hundreds or thousands of independent sequences.

At a high level, AF3 works by transforming biological sequences into a structured representation (MSAs, templates, and features) and then using a diffusion-based deep learning model to iteratively refine the predicted 3D structure. The data pipeline produces all the information AF3 needs to understand evolutionary context, residue co-variation, structural homology, and template similarity. The inference pipeline uses that information, together with model weights, to construct atomic-level predictions and associated confidence metrics.

On CHTC, these two stages map naturally to the heterogeneous compute environment:

 - **Stage 1 (Data Pipeline)** runs on CPU-only execution points and is typically the bottleneck for large batches of sequences.

 - **Stage 2 (Inference Pipeline)** runs on GPU-enabled execution points in the GPU Lab and GPU Open Capacity.

Separating these steps allows researchers to mix-and-match resource allocation strategies, reuse MSA/template features for multiple inference runs, and test alternative model parameters without re-running alignment-heavy jobs.

### The CPU-Only Pipeline: Generating Alignments (Step 1)

The first stage of AlphaFold3 prepares all input features needed for structure prediction. This step is entirely CPU-driven and dominated by database searches and feature construction.

#### What the Data Pipeline Does

AF3’s data pipeline performs the following major tasks:

   1. MMseqs2 Search<br>
   Uses MMseqs2 to search large protein sequence databases (e.g., UniRef, MGnify) to build MSAs.<br>
   _This provides evolutionary depth and signals conserved interactions._

   2. HMMER Search<br>
   Runs profile HMM searches to detect more remote homology and increase sensitivity for sequences with few known relatives.

   3. Template Retrieval (Optional)<br>
   Searches structural databases (e.g., PDB mmCIF) for homologous structures.<br>
   Templates provide complementary spatial constraints and can dramatically improve prediction accuracy when good matches exist.

   4. Feature Construction and Packaging <br>
   AF3 organizes MSAs, templates, and metadata into a standardized directory and packages them into:
       ```bash
       <job>.data_pipeline.tar.gz
       ```

#### Why This Runs on CPU Machines at CHTC

* Database searches are CPU-intensive and parallelize extremely well.
* Storage needs are substantial (~750 GB of AF3 databases), so CHTC provides pre-staged databases on select nodes.
* Running these jobs on GPU nodes would unnecessarily waste GPU time, increase GPU queue pressure, and reduce throughput.

This stage can be run for dozens, hundreds, or thousands of sequences simultaneously, making it ideal for HTC workloads.

## The GPU-Accelerated Pipeline: Structural Prediction (Step 2)

Once the data pipeline has produced MSAs and templates, AF3’s second stage uses this information to generate atomic-resolution structural models.

#### What the Inference Pipeline Does

The inference stage loads:

* precomputed MSAs and template features (from Step 1)
* model weights (user-supplied)
* AF3 container environment

It then performs:

1. Model Initialization & Tokenization <br>
The system converts all inputs into internal tokens for attention-based processing.<br>
Token counts relate directly to GPU memory demands.

2. Diffusion-Based Structure Prediction<br>
AF3 iteratively denoises and refines a 3D structure using diffusion and attention mechanisms.<br>
This is the most GPU-intensive part of the workflow.

3. Per-Seed Sampling & Ranking<br>
AF3 generates multiple predictions per job (depending on settings) and outputs:
   * ranked structures
   * confidence scores (e.g., pLDDT, PAE-like metrics)
   * intermediate trajectory snapshots (optional)

4. Output Packaging<br>
Results are tarballed as:
    ```bash
    <job>.inference_pipeline.tar.gz
    ```
#### Why This Runs on GPU Machines at CHTC

The inference pipeline requires:

* large matrix multiplications on GPUs
* high GPU memory capacity for complexes and nucleic acid assemblies
* JAX/XLA runtime optimizations available only on GPU-enabled nodes

* CHTC’s GPU Lab and GPU Open Capacity provide the necessary range of GPUs and allow users to scale inference jobs across many machines in parallel.

## Running AlphaFold3 on CHTC

### Set Up Your Software Environment
CHTC maintains a shared Apptainer container for AlphaFold3, which we **highlyrecommend most researchers use on CHTC's systems**. However, if you wish to build your own AlphaFold3 container (for example, to include custom models or software versions), follow the steps below to create your own Apptainer container image.

1. On your local machine, clone the AlphaFold3 repository:

    ```bash
    git clone https://github.com/google-deepmind/alphafold3.git
    ```
2. Navigate to the `alphafold3/` directory:

    ```bash
    cd alphafold3/
    ```
   
3. Build a docker image using the provided `Dockerfile`:

    ```bash
    docker build -t alphafold3:latest ./docker/
    ```
   
4. Push the docker image to a container registry (e.g., Docker Hub, Google Container Registry):

    ```bash
    docker tag alphafold3:latest <your-dockerhub-username>/alphafold3:latest
    docker push <your-dockerhub-username>/alphafold3:latest
    ```
   
5. On your CHTC Access Point, pull the docker image and convert it to an Apptainer image:

    ```bash
    apptainer build alphafold3.sif docker://<your-dockerhub-username>/alphafold3:latest
   ```
   
6. Verify that the Apptainer image was created successfully:

    ```bash
    apptainer exec alphafold3.sif python3 -c "import alphafold3; print(alphafold3.__version__)"
   ```

### Data Wrangling and Preparing AlphaFold3 Inputs

In this section, we will prepare the input files and job directories needed to run AlphaFold3 on CHTC. This includes creating input JSON files for each protein sequence and organizing them into job-specific directories.

#### Setting Up AlphaFold3 Input JSONs and Job Directories
Alphafold3 requires input JSON files that contain the input query sequence(s) along with their corresponding metadata, such as chain IDs, sequence names, and molecule type. Additionally, if you are using precomputed MSAs and templates, these JSON files should also reference the paths to them. For this tutorial, we will create individual JSON files for each protein sequence. Each AlphaFold3 job will have its own directory containing the input JSON file and any associated supporting files, such as

```bash
./AF3_Jobs/Job3_ProteinZ/
├── data_inputs                            # data_inputs directory for the data pipeline
│   └── fold_input.json              # input JSON for data pipeline
└── inference_inputs                       # inference_inputs directory for the GPU inference pipeline
```

As your jobs progress through the two phases of the AlphaFold3 workflow, the output files will be organized as follows:

```bash
./AF3_Jobs/Job3_ProteinZ/
├── data_inputs                            # data_inputs directory for the data pipeline
│   └── fold_input.json              # input JSON for data pipeline
├── inference_inputs                       # inference_inputs directory for the GPU inference pipeline
│   └── node9.data_pipeline.tar.gz   # tar.gz output from data pipeline (input for inference pipeline)
├── data_pipeline_job.err                  # stderr from data pipeline runs
├── data_pipeline_job.out                  # stdout from data pipeline runs
├── inference_pipeline_job.err             # stderr from inference pipeline runs
├── inference_pipeline_job.out             # stdout from inference pipeline runs
└── node9.inference_pipeline.tar.gz        # tar.gz output from inference pipeline (final results)
```

We need to setup a job directory for each AlphaFold3 job we plan to run. Each job directory will contain the necessary input JSON files and will be used to store the output files generated during the two phases of the AlphaFold3 workflow. We can use a simple bash script to automate the creation of these job directories and input JSON files.

1. Setup a CSV manifest containing your protein sequences in FASTA format in the `data/protein_sequences/` directory. Each FASTA file should contain a single protein sequence. The CSV should follow this format: 

    ```bash
    job_name,molecule_type,chain_id,sequence
    ProteinA,protein,A,MKTAYIAKQRQIS
    ProteinB,protein,A,GAVLILALLAVF
    ```
    
    If you plan to model multiple molecules, simply add more columns to the CSV manifest. For example, if you have two proteins and one DNA molecule, your CSV might look like this:
    
    ```bash
    job_name,molecule_type,chain_id,sequence
    ProteinA,protein,A,MKTAYIAKQRQIS,protein,B,MKTAYIAKQRQIS
    ProteinB,protein,A,GAVLILALLAVF,dna,B,GCGTACGTAGCTAGC
    ```
    
    You can also model multimeric complexes by including multiple chain_ids in the CSV manifest. Each chain_id must be wrapped in double-quotes `(")` and the full chain_id set must be wrapped in single quotes `(')`. For example, if you have a trimeric protein complex, your CSV might look like this:
    
    ```bash
    job_name,molecule_type,chain_id,sequence
    ProteinComplex,protein,'"A", "B", "C"',MKTAYIAKQRQIS
    ```

2. Run the helper script in `scripts/setup_af3_jobs.sh` to read the CSV manifest and create the necessary job directories and input JSON files for each protein sequence.:

    ```bash
   ./scripts/setup_af3_jobs.sh --manifest data/protein_sequences/manifest.csv --output_dir AF3_Jobs/
    ```

3. Verify that the job directories and input JSON files have been created correctly:

    ```bash
   tree AF3_Jobs/
   ```

   You should see a list of job directories corresponding to each protein sequence in your CSV manifest.

#### Preparing Your _List of (AlphaFold) Jobs_

Often you will hear the term "list of jobs" when working with HTCondor. A "list of jobs" is a simple way to specify multiple jobs to be run using a single HTCondor submit file. If you are familar with HTCondor, you may have used the `queue` command in your submit files to specify multiple jobs. HTCondor uses this `queue` statement to process multiple job submissions in one single submit file. Each job is defined by a set of arguments passed to the executable. This constitutes a "list of jobs".

There are many ways to create a "list of jobs" directly in your HTCondor submit files. For example, you can use the `queue` command with a `from` clause to read job arguments from a file or use the `queue` command with a `matching files` clause to generate jobs based on files in a directory. However, when working with large numbers of jobs, it is often more convenient to use a separate file to list the jobs you want to run. This file is referred to as a "list of jobs".

 In this tutorial, we will create a "list of jobs" file that contains the names of each AlphaFold3 job directory we plan to run. This file will be used in our HTCondor submit files to specify which jobs to execute.

1. Generate the "list of jobs" file by listing the job directories in your `AF3_Jobs/` directory:

    ```bash
    ls AF3_Jobs/ > list_of_af3_jobs.txt
    ```

1. Examine the contents of your "list of jobs" file:

    ```bash
    cat list_of_af3_jobs.txt
    ```

> [!NOTE]  
> Each line in `list_of_af3_jobs.txt` corresponds to a single AlphaFold3 job that will be executed using HTCondor. You can modify this file to add or remove jobs as needed. This file, functionally, is a one-column comma-separated values (CSV) file where each line represents a job name. You could add additional columns to this file if you wanted to pass more variables to your HTCondor jobs.

### Submit Your AlphaFold3 Jobs - CPU-Intensive Alignment Generation (Step 1)

The data-pipeline stage prepares all alignments, templates, and features needed for AF3 prediction. These CPU-only jobs run on CHTC’s standard compute nodes and can be scaled to many sequences at once. In the steps below, you’ll submit one data-pipeline job per sequence, producing the feature tarballs required for the GPU inference stage.

[INSERT GRAPHIC OF ALPHAFOLD3 WORKFLOW WITH DATA PIPELINE HIGHLIGHTED]

#### AlphaFold3 Databases Availability on CHTC

CHTC maintains a full, pre-extracted copy of the AlphaFold3 reference databases on a subset of CPU execute points. When your data-pipeline jobs match to one of these machines, they can use the local /alphafold3 directory directly, avoiding the costly transfer and extraction of several hundred gigabytes of database files. This dramatically reduces startup time, disk requirements, and overall job runtime. If a job lands on a machine without pre-staged databases, the script automatically falls back to unpacking the databases in the job’s scratch space, ensuring that every job can run regardless of where it matches.

You can target these pre-staged database nodes specifically by adding the following requirement to your submit file:

```bash
requirements = (HasAlphafold3 == true)
```

1. Change to your `tutorial-CHTC-AF3/` directory:
    ```bash
    cd ~/tutorial-CHTC-AF3/
   ```

2. Review your Data Pipeline executable script `scripts/data_pipleine.sh`. Generally, no changes will be necessary. You should, however, review the following information and options below as your AF3 jobs may require additional non-default options.

    #### Overview: AlphaFold3 Data Pipeline Executable (Data-Only Stage)

    This script implements the **data-generation portion** of the AlphaFold3 workflow — the part that creates alignments, MSA features, and template features required for structure prediction. It **does not** run structure inference step. It is designed for execution inside HTCondor jobs on CHTC systems and supports both containerized and pre-staged database workflows.
    
    The script handles:

   - Preparing a clean per-job working directory  
   - Detecting database availability (local `/alphafold3` vs. copying from staging)  
   - Copying and/or binding input JSON files  
   - Running the AlphaFold3 **data pipeline only**  
   - Tarballing the resulting `data_pipeline` outputs for transfer back to the Access Point  
   - Cleaning up the job sandbox to minimize disk use

    #### 1. Checking for Pre-Staged AF3 Databases

    CHTC hosts a copy of the AlphaFold3 databases locally on certain machines. Machines with these databases advertise this resource availability using the `HasAlphafold3` HTCondor MachineAd. The script is able to run on machines with/without these databases pre-loaded. The script inspects `.machine.ad` to see whether the matched machine advertises the availbility of AF3 databases:
    
    ```
    HasAlphafold3 = true
    ```
    
    The script also checks that the `/alphafold3/` path is not empty. If the MachineAd is set to `true` **and** `/alphafold3/` is populated:
    
    ✔ Uses the pre-extracted AF3 databases  
    ✘ Otherwise falls back to decompressing databases locally.

    The databases hosted by CHTC include the full AF3 default database suite:

    - PDB mmCIF  
    - MGnify Clusters 
    - UniRef90  
    - UniProt (Full UniProt Database) 
    - BFD First Non-Consensus Sequences
    - RNAcentral  
    - NT RNA clusters  
    - Rfam
   
    Databases may be updated periodically by CHTC staff. If you require a custom database set (e.g., reduced-size databases for testing), you can modify the script to use your own database tarballs. Contact the CHTC Research Computing Facilitation team for assistance.

    #### 2. Command-Line Options
    
    | Flag | Meaning |
    |------|---------|
    | `--work_dir_ext <name>` | Name appended to working directory (`work.<name>`) |
    | `--tmpdir <path>` | Override scratch location |
    | `--verbose` / `--silent` | Adjust verbosity |
    | `--no_copy` | Don’t copy containers/databases locally |
    | `--container <image>` | Path to Apptainer image |
    | `--random_sleep_minutes <N>` | Randomized delay to prevent thundering herd |
    | `--smalldb` | Use a reduced-size database set |
    | `--extracted_database_path <path>` | Use provided database directory |

    #### 3. Working Directory Setup
    
    Creates:
    
    ```
    work.<ext>/
        af_input/
        af_output/
        models/
        public_databases/
        tmp/
    ```
    
    Moves all `*.json` into `af_input/`.
    
    #### 4. Running the AlphaFold3 Data Pipeline
    
    Inside the container:
    
    ```
    python run_alphafold.py --run_data_pipeline=true --run_inference=false
    ```
    
    Generates all alignment and feature files under `af_output/<job_name>/`.
     
    #### 5. Packaging Results
    
    Each output directory is archived:
    
    ```
    <target>.data_pipeline.tar.gz
    ```
    
    These are returned to the submit host.

3. Create your submit file `data_pipeline.sub`. The submit file below works out-of-the-box (if you've setup your directories as specified in section [Setting Up AlphaFold3 Input JSONs and Job Directories]()). You can specify additional parameters for the exectuable in the `arguments` attribute as needed. Refer to the [AlphaFold3 Data Pipeline Executable - Command-Line Options]() section above for available options.

    ```bash
    # CHTC maintained container for AlphaFold3 as of January 2025
    #container_image = osdf:///ospool/uw-shared/OSG-Staff/public/alphafold3/alphafold3.minimal.22Jan2025_v1.sif
    container_image = file:///staging/groups/glbrc_alphafold/af3/alphafold3.minimal.22Jan2025.sif
    
    executable = scripts/data_pipeline.sh
    
    log = ../logs/data_pipeline.log
    output = data_pipeline_$(Cluster)_$(Process).out
    error  = data_pipeline_$(Cluster)_$(Process).err
    
    initialdir = $(directory)
    # transfer all files in the data_inputs directory
    transfer_input_files = data_inputs/
   
   # transfer output files back to the submit node
   transfer_output_files = $(directory).data_pipeline.tar.gz
   transfer_output_remaps = "$(directory).data_pipeline.tar.gz=inference_inputs/$(directory).data_pipeline.tar.gz"
    
    should_transfer_files = YES
    when_to_transfer_output = ON_EXIT
    
    # We need this to transfer the databases to the execute node
    Requirements = (Target.HasCHTCStaging == true) && (TARGET.HasAlphafold3 == true)
    
    if defined USE_SMALL_DB
      # testing requirements
      request_memory = 8GB
      request_disk = 16GB
      request_cpus = 4
      arguments = --smalldb --work_dir_ext $(Cluster)_$(Process) --verbose
    else
      # full requirements
      request_memory = 8GB
      # Request less disk if matched machine already has AF3 DB preloaded (650GB savings)
      request_disk = 700000 - ( (TARGET.HasAlphaFold3?: 1) * 650000)
      request_cpus = 8
      arguments = --work_dir_ext $(Cluster)_$(Proc) 
    endif
    
    #queue directory matching job*
    queue directory from list_of_af3_jobs.txt
   ```

This submit file will read the contents of `list_of_af3_jobs.txt`, iterate through each line, and assign the value of each line to the variable `$(directory)`. This allows you to programmatically submit _N_ jobs, where _N_ equals the number of AlphaFold3 job directories you previously created. Each job processes one AlphaFold3 job directory and uses the CHTC-maintained AlphaFold3 container image, which is transferred to the Execution Point (EP) by HTCondor.

The submit files will attempt to match to machines that are advertising the `HasAlphafold3` resource, which ensures that the necessary AlphaFold3 databases are transferred to the EP for each job. You can remove the `&& (TARGET.HasAlphafold3 == true)` clause from the `Requirements` attribute if you want to run jobs on machines with and without pre-staged databases. Removing this clause will increase the number of available machines for your jobs, leading to faster _start_ times for your jobs, but may also increase job runtime due to long database transfer times.

The script will also check the matched machine's MachineAd, after the job has matched, to see if the `HasAlphafold3` attribute is set to `true`. If it is, the submit file will request significantly less disk space, as the databases are already present on the machine. If the attribute is not set to `true`, the submit file will request more disk space to accommodate the database transfer. Lower disk requests can lead to increased number of running jobs, as the scheduler has more flexibility in matching jobs to machines.

    > [!TIP]  
    > AlphaFold3 jobs can be resource-intensive, especially when using the very conserved query sequences. Conserved sequenced can generally **very deep alignments** which will require signifantly more memory. If you encounter out-of-memory errors during job execution, consider increasing the `request_memory` attribute in your submit file. You can also utilize the `retry_request_memory = <memory/expression>` command in your submit file to request a retry if the job holds for an out-of-memory error. For more information on how to use `retry_request_memory`, visit our [Request variable memory](https://chtc.cs.wisc.edu/uw-research-computing/variable-memory#use-retry_request_memory) documentation page.

```bash

4. Submit your data-pipeline jobs:

    ```
   condor_submit scripts/data_pipeline.sub
   ```

    > [!TIP]  
    > You can test the data pipeline using reduced-size databases by defining the `USE_SMALL_DB=1` variable when submitting your jobs. This is useful for debugging and testing purposes, as it reduces resource requirements and speeds up job execution. To use the small database set, submit your jobs with the following command: 
       > ```
       > condor_submit USE_SMALL_DB=1 data_pipeline.sub
       > ```

5. Track your job progress:

    ```
   condor_watch_q
   ```

### Submit Your AlphaFold3 Jobs - GPU-Accelerated Structural Prediction (Step 2)

Once the data-pipeline jobs have finished generating alignments and features, the next stage is to run the AlphaFold3 inference pipeline on GPU-enabled execute points. This stage loads the model weights, expands the feature tarball from Step 1, and performs the diffusion-based structure prediction. Because inference is GPU-intensive, these jobs run on CHTC’s GPU Lab and GPU Open Capacity, and can be scaled across many GPUs in parallel. In the steps below, you will set up a submit file that launches one inference job per sequence, each producing a final structure tarball ready for download and visualization.

This stage **does not** require the full AlphaFold3 databases, only the model weights and the feature tarballs produced in Step 1. As a result, you can run these jobs on a wider range of GPU Execute Points without worrying about database availability, including GPU EPs outside of CHTC on the OSPool. To learn more about using additional capacity beyond CHTC, visit our guide on [Scale Beyond Local HTC Capacity](https://chtc.cs.wisc.edu/uw-research-computing/scaling-htc). 

[INSERT GRAPHIC OF ALPHAFOLD3 WORKFLOW WITH INFERENCE PIPELINE HIGHLIGHTED]

1. Change to your `tutorial-CHTC-AF3/` directory:
    ```bash
    cd ~/tutorial-CHTC-AF3/
   ```

2. Review your Data Pipeline executable script `scripts/inference_pipleine.sh`. Generally, no changes will be necessary. You should, however, review the following information and options below as your AF3 jobs may require additional non-default options.

    #### Overview: AlphaFold3 Inference Pipeline Executable

    This script implements the **structure inference portion** of the AlphaFold3 workflow. The inference pipeline loads model weights, runs the deep learning model on previously generated MSA/template features, and produces predicted 3D structures. It **does not** generate alignments or template features. This script is designed for execution inside HTCondor jobs and supports both containerized workflows and user-specified model weight files.

    The script handles:

   - Preparing a clean per-job working directory on the GPU execute point
   - Expanding `.data_pipeline.tar.gz` results from the previous data pipeline stage 
   - Detecting and loading AlphaFold3 model parameters (**User supplied**) 
   - Running AlphaFold3 **inference only**  
   - Packaging predicted structures and metadata  
   - Cleaning up the job sandbox to minimize disk use

    #### 1. Model Parameters and Containers

    Unlike the data pipeline, the inference stage does **not need the full AlphaFold3 databases**. It only requires:

    - A valid AlphaFold3 **model parameter file** (`.bin.zst`)
    - A valid AlphaFold3 **container image** (if running outside the container)

    The script allows users to specify a model file using:

    ```
    --model_param_file <name_of_model_weights_file.zst>
    ```

    If the model weights are compressed (`.zst`), the script automatically decompresses them into the working directory. 

    #### 2. Command-Line Options

    | Flag | Meaning |
    |------|---------|
    | `--work_dir_ext <name>` | Name appended to working directory (`work.<name>`) |
    | `--verbose` / `--silent` | Adjust verbosity |
    | `--no_copy` | Do not copy container or model parameters locally |
    | `--container <image>` | Path to Apptainer/Singularity image |
    | `--model_param_file <path>` | Path to AlphaFold3 model weights |
    | `--user_specified_alphafold_options` | Additional arguments passed directly to `run_alphafold.py` |
    | `--enable_unified_memory` | Enable unified GPU memory mode for large jobs |

    #### 3. Working Directory Setup

    Creates:

    ```
    work.<ext>/
        af_input/
        af_output/
        models/
        public_databases/
    ```

    The inference stage expects the output from the data pipeline in the form of:

    ```
    *.data_pipeline.tar.gz
    ```

    These archives are extracted into `af_input/`, reconstructing the full feature directory structure required by AlphaFold3.

    #### 4. Preparing Model Weights

    Valid model weights may be provided in either uncompressed or `.zst` form. The script handles both:

    - If compressed:  
      ```
      zstd --decompress > models/af3.bin
      ```
    - If uncompressed:  
      ```
      cp af3.bin models/
      ```

    #### 5. GPU Capability Handling

    The script inspects the GPU compute capability using:

    ```
    python -c "import jax; print(jax.local_devices(backend='gpu')[0].compute_capability)"
    ```

    For GPUs with compute capability **7.x**, it automatically:

    - Forces flash attention to `xla`
    - Sets required XLA flags to avoid runtime errors

    Users may also enable unified memory mode if working with especially large complexes.

    #### 6. Running the AlphaFold3 Inference Stage

    Inside the container, the script executes:

    ```
    python run_alphafold.py  --run_data_pipeline=false --run_inference=true --input_dir=/root/af_input --model_dir=/root/models --output_dir=/root/af_output
    ```

    Users may append custom AF3 flags through:

    ```
    --user_specified_alphafold_options "<options>"
    ```
   
    This allows for advanced customization of the inference run. It can be used to modify parameters such as:
    - Enabling different size model compilation bucket sizes (`--buckets 5132, 5280, 5342`)
    - Altering or Disabling JAX GPU Memory Preallocation (useful for larger complexes over 8k tokens, see our section on [Handling Large Complexes]()
    - Other AlphaFold3-specific options
    For more information on available AlphaFold3 options, refer to the [AlphaFold3 GitHub Repository](https://github.com/google-deepmind/alphafold3/)

    Output includes:

    - Predicted 3D structures  
    - Ranking scores  
    - Seed/sample-level predictions  

    #### 7. Packaging Results

    Each prediction directory under `af_output/` is archived:

    ```
    <target>.inference_pipeline.tar.gz
    ```

    These tarballs are returned to the submit host for downstream use.

    #### 8. Cleanup

    After packaging the results, the script removes:

    - `work.<ext>/`  
    - Shell history files (`.bash_history`, `.lesshst`, etc.)

    This ensures minimal disk usage on execute nodes and prepares the job sandbox for automatic cleanup by HTCondor.

3. Create your submit file `inference_pipeline.sub`. You will need to edit the `MODEL_WEIGHT_PATH` **and** `gpus_minimum_memory`. You can specify additional parameters for the exectuable in the `arguments` attribute as needed. Refer to the [AlphaFold3 Data Pipeline Executable - Command-Line Options]() section above for available options.

    ```bash
    container_image = file:///staging/groups/glbrc_alphafold/af3/alphafold3.minimal.22Jan2025.sif
    
    executable = inference_pipeline.sh
    
    environment = "myjobdir=$(directory)"
    
    MODEL_WEIGHTS_PATH = /staging/damorales4/af3/weights/af3.bin.zst
    
    log = ../logs/inference_pipeline.log
    output = inference_pipeline_$(Cluster)_$(Process).out
    error  = inference_pipeline_$(Cluster)_$(Process).err
    
    initialdir = $(directory)
    # transfer all files in the inference_inputs directory
    transfer_input_files = inference_inputs/, $(MODEL_WEIGHTS_PATH)
    
    should_transfer_files = YES
    when_to_transfer_output = ON_EXIT
    
    request_memory = 8GB
    # need space for the container (3GB) as well
    request_disk = 10GB
    request_cpus = 4
    request_gpus = 1
    
    # we should be able to run short jobs on CUDA_CAPABILITY=7.x but need
    # other environment variables and options to be set. This is done automatically
    # in inference_pipeline.sh
    gpus_minimum_memory = 0
    
    # short jobs 4-6 hours so it is okay to use is_resumable
    +GPUJobLength = "short"
    +WantGPULab = true
    +is_resumable = true
    
    # Use --user-specified-alphafold-options to pass any extra options to AlphaFold3, such as
    # arguments = --model_param_file af3.bin.zst --work_dir_ext $(Cluster)_$(Process) --user-specified-alphafold-options "--buckets 5982"
    arguments = --model_param_file af3.bin.zst --work_dir_ext $(Cluster)_$(Process)
    
    queue directory from list_of_af3_jobs.txt   
   ```

This submit file will read the contents of `list_of_af3_jobs.txt`, iterate through each line, and assign the value of each line to the variable `$(directory)`. This allows you to programmatically submit _N_ jobs, where _N_ equals the number of AlphaFold3 job directories you previously created. Each job processes one AlphaFold3 job directory and uses the CHTC-maintained AlphaFold3 container image, which is transferred to the Execution Point (EP) by HTCondor.

The GPU-accelerated inference jobs do not require the full AlphaFold3 databases, so there is no need to check for pre-staged databases on the execute nodes. However, it does require the model weights to be transferred to the EP. The submit file above transfers the model weights specified in the `MODEL_WEIGHTS_PATH` variable to the EP for each job.

Beyond the model weights, your AlphaFold3 jobs may require different resource requests depending on the size and complexity of the proteins you are modeling. Typically, this mainly involves adjusting the `request_disk` and `gpus_minimum_memory` attributes in your submit file. Estimating the required GPU memory can be challenging, as it depends on factors such as sequence length, number of chains, and model parameters. As a starting point, you can refer to the following general guidelines:

| Tokens     | Estimated GPU Memory Requirement |
|------------|----------------------------------|
| Up to 1200 | 8-10 GB                          |
| 1200-1850  | 15-20 GB                         |
| 2000-3000  | 35-40 GB                         |
| Over 3000  | 70+ GB                           |

You can generally estimate number of tokens as approximately equivalent to 1.2x the sequence length. For example, a protein with 1000 amino acids would have roughly 1200 tokens. If you are modeling multimeric complexes, you will need to account for the combined sequence lengths of all chains. For structures that include nucleic acids (DNA/RNA), each nucleotide is considered one token. Post-translational modifications (PTMs) and additional ligands may also increase token counts, depending on the specific modifications involved. 

For very large complexes exceeding 10k tokens, you may need to enable unified memory mode using the `--enable_unified_memory` flag in the executable script. This allows AlphaFold3 to utilize system RAM in addition to GPU memory, which can help accommodate larger models. However, it may also lead to **significantly slower performance** due to increased data transfer times between system RAM and GPU memory. You will also need to ensure that the execute node has sufficient system RAM to support unified memory mode by increasing the `request_memory` attribute in your submit file. We recommend reaching out to the CHTC Research Computing Facilitation team for assistance with very large complexes.

    > [!TIP]  
    > CHTC has a number of smaller RTX series GPUs (e.g., RTX 4000, RTX 5000) with only 8-16GB of GPU memory. If your job requires less than 10GB of GPU memory, you can set `gpus_minimum_memory = 10000` in your submit file to allow your jobs to match to these smaller GPUs. This can help increase the number of available machines for your jobs, leading to faster start times.
    > 
    > If your jobs require more GPU memory than is available on these smaller GPUs, you can adjust the `gpus_minimum_memory` attribute in your submit file to request machines with larger GPUs (e.g., A100, V100). For example, to request machines with at least 32GB of GPU memory, you can set `gpus_minimum_memory = 32000` in your submit file. This will help ensure that your jobs are matched to machines with sufficient GPU memory to run successfully.

```bash

4. Submit your data-pipeline jobs:

    ```
   condor_submit scripts/inference_pipeline.sub
   ```

5. Track your job progress:

    ```
   condor_watch_q
   ```


### Visualize Your AlphaFold3 Results

AlphaFold3 generates a variety of output files, including predicted 3D structures in PDB format, ranking scores, and seed/sample-level predictions. These files are packaged into tarballs named `<job_name>.inference_pipeline.tar.gz` and are returned to the submit host for downstream analysis. You will need to download these results from the Access Point to your local machine for visualization and further analysis.

1. Once your inference pipeline jobs have completed successfully, `exit` the Access Point SSH session and return to your local machine.:

    ```bash
    [bbdager@ap2002 tutorial-CHTC-AF3]$ exit
    logout
    Shared connection to ap2002.chtc.wisc.edu closed.
    Bucky@MyLaptop ~ % 
   ```

    **_Note:_** Notice the chance in the prompt, indicating you are back on your local laptop. 

2. For each job directory, download and extract the `<job_name>.inference_pipeline.tar.gz` file to your local machine:

    ```bash
    scp <netID>@ap2002.chtc.wisc.edu:~/tutorial-CHTC-AF3/AF3_Jobs/Job1_ProteinA/Job1_ProteinA.inference_pipeline.tar.gz ./
    tar -xzvf Job1_ProteinA.inference_pipeline.tar.gz
    ```
   
3. After extracting the tarballs, you will find the predicted structures and associated metadata in each job directory under `af_output/<job_name>/`. You can visualize the predicted 3D structures using molecular visualization software such as PyMOL, Chimera, or VMD. To visualize a predicted structure using PyMOL, you can use the following command:

    ```bash
    pymol AF3_Jobs/Job1_ProteinA/af_output/Job1_ProteinA/predicted_structure.pdb
    ```

![pymol.png](.images/pymol.png)

## Next Steps

Now that you've successfully run the full AlphaFold3 two-stage workflow on the CHTC GPU capacity, you’re ready to extend this workflow to your own research projects, larger datasets, and more complex biomolecular systems. Below are recommended next steps to build on this tutorial.

🧬 Apply the Workflow to Your Own Data
* Replace the tutorial sequences with your own proteins, RNA molecules, complexes, or mixed multimers.
* Use the provided CSV manifest system to generate hundreds or thousands of AF3 job directories automatically.
* Experiment with:
  * Multimeric protein assemblies 
  * Protein–RNA or protein–DNA complexes 
  * Ligands or modified residues (e.g., 2′-O-methyl RNA, PTMs)
* Re-use the data pipeline outputs to run multiple inference configurations (different seeds, different AF3 options) without recomputing MSAs.

🚀 Run Larger Analyses
* Once you’re comfortable with the basics, try:
  * Large protein–RNA complexes (>10k tokens)
  * Multi-seed inference strategies 
  * Modeling structural evolution by comparing AF3 predictions across species 
* Integrating AF3 predictions into:
  * Molecular dynamics (MD) or Docking simulations 
  * Variant impact analyses 
  * Evolutionary analyses 
* Using AF3 for structural annotation of genomes, e.g., predicting full gene families with DAGMan workflows.

🧑‍💻 Get Help or Collaborate
* Reach out to [chtc@cs.wisc.edu](mailto:chtc@cs.wisc.edu) for one-on-one help with scaling your research.
* Attend office hours or training sessions—see the [CHTC Help Page](https://chtc.cs.wisc.edu/uw-research-computing/get-help.html) for details.

### Software

In this tutorial, we created several *starter* Apptainer containers, including tools like: Dorado, SAMtools, Minimap, and Sniffles2. These containers can serve as a *jumping-off* for you if you need to install additional software for your workflows. 

Our recommendation for most users is to use "Apptainer" containers for deploying their software.
For instructions on how to build an Apptainer container, see our guide [Using Apptainer/Singularity Containers](https://chtc.cs.wisc.edu/uw-research-computing/apptainer-htc#main).
If you are familiar with Docker, or want to learn how to use Docker, see our guide [Using Docker Containers](https://chtc.cs.wisc.edu/uw-research-computing/docker-jobs#main).

This information can also be found in our guide [Using Software on CHTC](https://chtc.cs.wisc.edu/uw-research-computing/software-overview-htc#main).

### Data

AlphaFold3 jobs involve large and complex datasets, especially during the data pipeline stage. Understanding how data moves through the HTC system—and how to store it efficiently—is essential for scaling AF3 workloads.

#### Key AF3 data components include:
* AF3 reference databases (~750 GB unpacked in total)
  * Stored on selected CHTC nodes and accessed via `(TARGET.HasAlphafold3 == true)` requirements
  * Avoids expensive per-job transfers of MGnify, UniRef, UniProt, PDB mmCIF, Rfam, RNAcentral, etc.
* Model weights (user-supplied)
  * Typically ~2–4 GB depending on compression
  * Stored in your `/staging` directory for fast transfer to GPU execute nodes
* Data pipeline outputs
  * `.data_pipeline.tar.gz` for each job
  * Typically 100–300 MB depending on alignment depth
* Inference outputs
  * `.inference_pipeline.tar.gz` per job
  * Includes ranked PDBs, metrics, seeds, confidence scores

For guides on how data movement works on the HTC system, see our [Data Staging and Transfer to Jobs](https://chtc.cs.wisc.edu/uw-research-computing/htc-job-file-transfer) guides.

### GPUs

AlphaFold3 inference is GPU-intensive, and selecting the right GPU resources is crucial for reliability and efficiency. CHTC provides a broad range of available GPUs, from smaller RTX-series cards (8–16 GB) to large-memory accelerators such as A100s (40–80 GB).

#### Key AF3 GPU considerations:

Token count drives GPU memory needs
  * Tokens ≈ ~1.2 × total sequence length across all chains
  * RNA/DNA bases count as tokens
  * Larger complexes → more GPU memory
Typical GPU memory guidelines:
  * | Tokens     | Estimated GPU Memory Requirement |
    |------------|----------------------------------|
    | Up to 1200 | 8-10 GB                          |
    | 1200-1850  | 15-20 GB                         |
    | 2000-3000  | 35-40 GB                         |
    | Over 3000  | 70+ GB                           |
* Unified memory support
  * Can rescue very large complexes (>10k tokens)
  * Significantly slower
  * Increase `request_memory` if enabling `--enable_unified_memory` in your executable arguments

If you would like to learn more about our GPU capacity, please visit our [GPU Guide on CHTC Documentation Portal](https://chtc.cs.wisc.edu/uw-research-computing/gpu-jobs).

## Getting Help

The CHTC Research Computing Facilitators are here to help researchers using the CHTC resources for their research. We provide a broad swath of research facilitation services, including:

* **Web guides**: [CHTC Guides](https://chtc.cs.wisc.edu/uw-research-computing/htc/guides.html) - instructions and how-tos for using the CHTC cluster.
* **Email support**: get help within 1-2 business days by emailing [chtc@cs.wisc.edu](mailto:chtc@cs.wisc.edu).
* **Virtual office hours**: live discussions with facilitators - see the [Email, Office Hours, and 1-1 Meetings](https://chtc.cs.wisc.edu/uw-research-computing/get-help.html) page for current schedule.
* **One-on-one meetings**: dedicated meetings to help new users, groups get started on the system; email [chtc@cs.wisc.edu](mailto:chtc@cs.wisc.edu) to request a meeting.

This information, and more, is provided in our [Get Help](https://chtc.cs.wisc.edu/uw-research-computing/get-help.html) page.
