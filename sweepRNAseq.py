# sweepRNAseq.py

"""
sweepRNAseq Pipeline

An object-oriented Python pipeline designed to remove DNA and rRNA contamination from RNA-seq datasets.

The pipeline performs the following steps:

1. Runs FastQC for initial quality control and contamination check.
2. Removes rRNA and genomic contamination using SortMeRNA
3. Runs FastQC again on the cleaned reads to check for contamination.
4. Aligns reads to the reference genome using STAR to detect potential contamination.
5. Verifies cleanup by re-running FastQC on filtered reads.

The program assumes the tools fastqc, STAR, and sortmerna are installed in the "sweepRNAseq-env" conda environment.

Usage example: python sweepRNAseq.py sample1.fastq sample2.fastq annotation.gtf reference-genome.fa rrna-database.fa

"""

import subprocess
from pathlib import Path # lets us safely handle file paths
import sys   # allows the program to run terminal commands
from typing import Union, List


class FastQC:
    """
    Small helper class to run FastQC from Python instead of typing
    long commands in the terminal.

    This class is meant to be used inside a larger RNA-seq pipeline.
    It takes a list of FASTQ/FASTQ.GZ files, checks that they exist,
    creates an output directory for the reports, builds the FastQC
    command, and then runs it with subprocess.
    """

    def __init__(
        self,
        inputFiles: Union[str, List[str]],
        outputDir: Union[str, Path],
        fastqcExecutable: str = "fastqc",
        threads: int = 2
    ) -> None:
        """
        Initialize a FastQC job.

        Parameters are

        inputFiles :
            Either a single FASTQ file (string / Path) or a list of FASTQ files.
            These are typically raw or trimmed reads, in .fastq or .fastq.gz format.

        outputDir :
            Directory where FastQC should write the HTML and ZIP reports.

        fastqcExecutable :
            Name of the FastQC executable. By default this is just "fastqc"
            assuming it is on the user's PATH.

        threads :
            Number of CPU threads to use. More threads generally means faster runs.
        """

        # Normalize inputFiles so that we always store a list of Path objects.
        if isinstance(inputFiles, (str, Path)):
            # Single file case so it wraps into a list for consistent handling downstream.
            self.inputFiles = [Path(inputFiles)]
        else:
            # Multiple files case to convert each item into a Path object.
            self.inputFiles = [Path(f) for f in inputFiles]

        # Store output directory as a Path object for convenient filesystem operations.
        self.outputDir = Path(outputDir)

        # Name/path of the FastQC executable.
        self.fastqcExecutable = fastqcExecutable

        # Number of threads to pass to FastQC.
        self.threads = threads

    def CheckInputs(self) -> None:
        """
        Verify that all input FASTQ files actually exist on disk.

        This guards against confusing FastQC errors later on.
        If any file is missing, we raise a clear FileNotFoundError
        listing the problematic files.
        """

        # Collect any input files that do not exist.
        missing = [str(f) for f in self.inputFiles if not f.is_file()]

        # If we found at least one missing file, fail early with a helpful message.
        if missing:
            raise FileNotFoundError(
                f"These FASTQ files were not found: {', '.join(missing)}"
            )

    def PrepareOutputDir(self) -> None:
        """
        Ensures that the output directory exists before running FastQC.

        If the directory does not exist,
        it is created. If it already exists, nothing happens.
        """

        # mkdir with parents=True recursively creates any missing parent directories.
        # exist_ok=True avoids raising an error if the directory is already there.
        self.outputDir.mkdir(parents=True, exist_ok=True)

    def BuildCommand(self) -> List[str]:
        """
        Construct the FastQC command as a list of strings.
        """

        # Base command with executable, thread count, and output directory.
        cmd = [
            self.fastqcExecutable,           # usually "fastqc"
            "--threads", str(self.threads),  # convert to string for subprocess
            "--outdir", str(self.outputDir)  # where FastQC should write results
        ]

        # Append all FASTQ file paths to the command.
        # FastQC can process multiple files in a single call.
        cmd.extend(str(f) for f in self.inputFiles)

        return cmd

    def Run(self) -> None:
        """
        Execute FastQC using the stored settings.

        Pipeline of steps:
        1. Check that all input files exist.
        2. Create the output directory if needed.
        3. Build the FastQC command.
        4. Run the command with subprocess.run().
        5. Print a short status message so the user knows where outputs are.

        Any FastQC failure (non zero exit code) will raise a CalledProcessError
        because we use check=True.
        """

        # 1. Fail early if any input is missing.
        self.CheckInputs()

        # 2. Make sure the output directory is ready.
        self.PrepareOutputDir()

        # 3. Build the FastQC command for the given inputs and settings.
        cmd = self.BuildCommand()

        # Show the exact command being run is useful for debugging or reproducing.
        print(f"[FastQC] Running: {' '.join(cmd)}")

        # 4. Call FastQC via subprocess. If FastQC fails, this will raise an error.
        subprocess.run(cmd, check=True)

        # 5. Confirm completion and print the absolute path of the output directory.
        print(f"[FastQC] Done! Reports saved in: {self.outputDir.resolve()}")


class SortMeRNA:
  '''
  A class for running SortMeRNA on FASTQ files.
  SortMeRNA is a tool used for filtering ribosomal RNA (rRNA) from RNA-seq datasets.

  This class builds a SortMeRNA index from a reference rRNA FASTA file and filters rRNA reads from single-end or paired-end FASTQ datasets.

  This class is written from SortMeRNA version 4.3.7 which automatically builds the index the first time you run it with --ref.
  '''
  def __init__(
      self,
      refFasta: Path | str,
      indexPrefix: Path | str,
      executable: str = "sortmerna",
      indexExecutable: str = "indexdb_rna",
      threads: int = 10,
  ):
      '''
      Initialize a SortMeRNA helper object.
            refFasta:
                Path to the rRNA reference FASTA file used by SortMeRNA.
            indexPrefix:
                Path prefix where index files would be stored. With 4.3.7 this
                is mostly kept for consistency and future flexibility.
            executable:
                Name or path of the main SortMeRNA executable (default "sortmerna").
            indexExecutable:
                Name or path of the index-building executable (default "indexdb_rna").
                Not actively used in 4.3.7, but stored for backwards compatibility.
            threads:
                Number of CPU threads SortMeRNA should use (default 10).
      '''
      # Turn the reference FASTA into a Path object so it is easy to work with.
      self.refFasta = Path(refFasta)
      # Path prefix for index files. Even though indexing is automatic in
      # this version, we keep this attribute to match older usage patterns.
      self.indexPrefix = Path(indexPrefix)
      # Store the executable names used when we build the SortMeRNA commands.
      self.executable = executable
      self.indexExecutable = indexExecutable
      # Store the threads as a string, because subprocess.run expects
      # arguments as strings, not integers.
      self.threads = str(threads)

  def buildIndex(self) -> None:
      """Build SortMeRNA index from the rRNA reference FASTA."""
      # SortMeRNA 4.3.7 builds indexes automatically, so this is a no-op
      # The index will be created when sortmerna runs with --ref
      print(f"[SortMeRNA] Index will be built automatically on first use")
      # We keep the function body as 'pass' because there is no action to take.
      pass

  def filterRRNA(
      self,
      reads1: Path,
      reads2: Path | None = None,
      alignedPrefix: Path | str = "rrna_reads",
      otherPrefix: Path | str = "clean_reads",
  ) -> tuple[Path, Path]:

      """
        Run SortMeRNA to separate rRNA and non-rRNA reads.

        Args:
            reads1:
                Path to the first read file (R1) in FASTQ or FASTQ.GZ format.
            reads2:
                Optional path to the second read file (R2) for paired-end data.
                If None, SortMeRNA is run in single-end mode.
            alignedPrefix:
                Prefix for the rRNA-aligned output files (rRNA reads).
            otherPrefix:
                Prefix for the non-rRNA ("clean") output files.

        Returns:
            A tuple (rRNAOut, cleanedFASTQ) where:
                rRNAOut:
                    Path to the FASTQ (or FASTQ.GZ) file containing rRNA reads.
                cleanedFASTQ:
                    Path to the FASTQ (or FASTQ.GZ) file containing
                    non-rRNA (cleaned) reads.

        Raises:
            FileNotFoundError:
                If the expected cleaned output file cannot be found after
                SortMeRNA finishes.
      """
      # Convert the output prefixes to Path objects so we can easily
      # work with parent directories, stems, etc.
      alignedPrefix = Path(alignedPrefix)
      otherPrefix = Path(otherPrefix)

      # SortMeRNA 4.3.7 uses --ref with just the FASTA file (index built automatically)
      # Use unique workdir per file to avoid conflicts
      sortmernaWorkDir = alignedPrefix.parent / f"workdir_{Path(reads1).stem}"
      # Make the directory if it does not exist yet. exist_ok=True means:
      # "If it already exists, do not crash, just continue.""
      sortmernaWorkDir.mkdir(exist_ok=True)

      # Build the base SortMeRNA command as a list of strings.
      # We use a list instead of a single string so subprocess.run can safely
      # pass each element as a separate argument.
      cmd = [
          self.executable,
          "--ref",
          str(self.refFasta),
          "--reads",
          str(reads1),
          "--threads",
          self.threads,
          "--aligned",
          str(alignedPrefix),
          "--other",
          str(otherPrefix),
          "--fastx",
          "--log",
          "--workdir",
          str(sortmernaWorkDir),
        ]

      if reads2 is not None:
          # add second reads and paired flag
          cmd.insert(cmd.index("--threads"), str(reads2))
          cmd.insert(cmd.index("--threads"), "--reads")
          cmd.append("--paired_in")

      print("Running: sortmerna", " ".join(cmd[1:]))
      subprocess.run(cmd, check=True)

      # Find output files - SortMeRNA 4.3.7 creates: inputname_cleaned.fq.gz
      reads1Name = Path(reads1).name.replace(".gz", "")
      outputDir = otherPrefix.parent
      cleanedFASTQ = next((f for f in outputDir.glob(f"{reads1Name}_cleaned.fq*")), None)
      if not cleanedFASTQ:
        cleanedFASTQ = next((f for f in outputDir.glob("*_cleaned.fq*")), None)
      if not cleanedFASTQ:
        raise FileNotFoundError(f"Could not find SortMeRNA cleaned output in {outputDir}")

      alignedDir = alignedPrefix.parent
      rRNAOut = next((f for f in alignedDir.glob(f"{alignedPrefix.name}.*") if f.suffix in [".fq", ".fastq"] or f.suffixes[-1] == ".gz"), alignedPrefix.with_suffix(".fastq"))

      return rRNAOut, cleanedFASTQ

class STAR:
  """
  STAR is a tool for aligning RNA-seq reads to a reference genome.

  Inputs:

  - FASTQ files (can be a single file or a list of files).
  - Reference genome (FASTA file)
  - Output directory

  Outputs:
  - Aligned reads (SAM file).
  - Log file (detailed information about the reads and alignment process).
  - Output directory (contain the SAM file and the log file).
  """

  def __init__(self, inputFASTQ, referenceGenome, outputDirectory, annotationFile=None, threads=10):
    """
    Initializes the STAR object with the genome directory and the number of threads to use.

    Parameters:
    - genomeDirectory: the directory containing the genome files.
    - threads: the number of threads to use for the alignment process (default: 10 CPU cores).

    The genome directory should contain the following files:
    - reference genome (FASTA file)
    - annotation file (GTF file)
    - genome index (STAR index)
    """
    if isinstance(inputFASTQ, str):
      self.inputFASTQ = [inputFASTQ]
    else:
      self.inputFASTQ = inputFASTQ
    self.referenceGenome = Path(referenceGenome)
    self.outputDirectory = Path(outputDirectory)
    self.annotationFile = Path(annotationFile) if annotationFile else None
    self.threads = threads
    self.genomeIndexDirectory = self.outputDirectory / "genomeIndex"
    self.mappedDirectory = self.outputDirectory / "mapped"

  def generateGenomeIndex(self):
    """
    Generates the genome index for the reference genome.

    Inputs:
    - Reference genome (FASTA file)
    - Annotation file (GTF file) - optional

    Outputs:
    - Genome index (STAR index) - stored in the output directory

    """
    if (self.genomeIndexDirectory / "Genome").exists():
      return
    self.genomeIndexDirectory.mkdir(parents=True, exist_ok=True)
    indexArguments = [
      "STAR",
      "--runMode", "genomeGenerate",
      "--runThreadN", str(self.threads), # the number of threads to use for the genome index generation
      "--genomeDir", str(self.genomeIndexDirectory), # the directory to store the genome index
      "--genomeFastaFiles", str(self.referenceGenome), # the reference genome to use for the genome index
    ]
    if self.annotationFile:
      # Add the annotation file and the sjdbOverhang argument to the index arguments
      # The sjdbGTFfile argument is the annotation file to use for the genome index
      # The sjdbOverhang argument is set to 100, which is the default value recommended from the STAR manual
      indexArguments.extend([
        "--sjdbGTFfile", str(self.annotationFile),
        "--sjdbOverhang", "100"
      ])
    subprocess.run(indexArguments, check=True)

  def alignReads(self):
    self.mappedDirectory.mkdir(parents=True, exist_ok=True)
    for fastqFile in self.inputFASTQ:
      fastqPath = Path(fastqFile) #
      inputName = fastqPath.stem
      outputPrefix = self.mappedDirectory / f"{inputName}_star"
      # STAR alignment arguments from the STAR manual (parameters set to default values recommended from the manual)
      alignArguments = [
        "STAR",
        "--runThreadN", str(self.threads),
        "--genomeDir", str(self.genomeIndexDirectory),
        "--readFilesIn", str(fastqFile),
        "--outFileNamePrefix", f"{outputPrefix}_",
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outputFilterType", "Normal",
        "--outFilterMultimapScoreRange", "1",
        "--outReadsUnmapped", "Fastx",
        "--outFilterMultimapNmax", "10",
        "--outFilterMismatchNmax", "10",
        "--outFilterMismatchNoverLmax", "0.3",
        "--outFilterMismatchNoverReadLmax", "1.0",
        "--outFilterScoreMin", "0",
        "--outFilterScoreMinOverLread", "0.66",
        "--outFilterMatchNmin", "0",
        "--outFilterIntronMotifs", "None",
        "--outFilterIntronStrands", "RemoveInconsistentStrands",
        "--alignEndsType", "EndToEnd",
      ]
      subprocess.run(alignArguments, check=True)

  def run(self):
    """
    Runs the STAR alignment process.

    Parameters:
    - fastqFile: the FASTQ file to align.
    - outputDirectory: the directory to store the aligned reads.

      Inputs:
      - FASTQ file (the reads to align)
      - Output directory (the directory to store the aligned reads)

    Outputs:
    - Aligned reads (SAM file)
    - Log file (detailed information about the reads and alignment process).
    """
    self.outputDirectory.mkdir(parents=True, exist_ok=True)
    self.generateGenomeIndex()
    self.alignReads()


class Pipeline:
  """
  The pipeline class coordinates the execution of the FastQC, STAR, and SortMeRNA steps sequentially.
  """
  def __init__(
    self,
    inputFASTQ: List[str],
    referenceGenome: str,
    rRNADatabase: str,
    outputDir: str = "output",
    annotationFile: str = None,
    threads: int = 10
  ):
    """
    Initializes the pipeline with input FASTQ files, reference genome, and rRNA database.

    Parameters:
    - inputFASTQ: List of FASTQ file paths
    - referenceGenome: Path to reference genome FASTA file
    - rRNADatabase: Path to rRNA reference FASTA file
    - outputDir: Output directory for all results (default: "output")
    - annotationFile: Optional GTF annotation file for STAR
    - threads: Number of threads to use (default: 4)
    """
    if len(inputFASTQ) == 0:
      raise ValueError("Error: No FASTQ files provided.")

    self.fastqFiles = [Path(filePath) for filePath in inputFASTQ]
    self.referenceGenome = Path(referenceGenome)
    self.rRNADatabase = Path(rRNADatabase)
    self.outputDir = Path(outputDir)
    self.annotationFile = Path(annotationFile) if annotationFile else None
    self.threads = threads

    # Create subdirectories for different outputs
    self.fastqcDir = self.outputDir / "fastqc"
    self.sortmernaDir = self.outputDir / "sortmerna"
    self.starDir = self.outputDir / "star"
    self.cleanedDir = self.outputDir / "cleaned"

  def run(self):
    """
    Runs the complete pipeline:
    1. FastQC on raw reads
    2. SortMeRNA to remove rRNA
    3. FastQC on cleaned reads
    4. STAR alignment on cleaned reads
    5. FastQC on aligned reads
    """
    # Initialize SortMeRNA
    rRNAIndexPrefix = self.sortmernaDir / "rrna_index"
    sortmerna = SortMeRNA(
      refFasta=self.rRNADatabase,
      indexPrefix=rRNAIndexPrefix,
      threads=self.threads
    )

    # SortMeRNA 4.3.7 builds indexes automatically, no need to pre-build
    self.sortmernaDir.mkdir(parents=True, exist_ok=True)

    for fastqFile in self.fastqFiles:
      print(f"\n[Pipeline] Processing: {fastqFile.name}")

      # Step 1: FastQC on raw reads
      print("[Pipeline] Step 1: Running FastQC on raw reads...")
      fastqcRaw = FastQC(
        inputFiles=[fastqFile],
        outputDir=self.fastqcDir / "raw" / fastqFile.stem,
        threads=2  # FastQC threading is for parallel files, 2 is optimal
      )
      fastqcRaw.Run()

      # Step 2: SortMeRNA to remove rRNA
      print("[Pipeline] Step 2: Removing rRNA contamination with SortMeRNA...")
      self.cleanedDir.mkdir(parents=True, exist_ok=True)
      alignedPrefix = self.sortmernaDir / f"{fastqFile.stem}_rrna"
      otherPrefix = self.cleanedDir / f"{fastqFile.stem}_cleaned"

      _, cleanedFASTQ = sortmerna.filterRRNA(
        reads1=fastqFile,
        alignedPrefix=alignedPrefix,
        otherPrefix=otherPrefix
      )

      # Step 3: FastQC on cleaned reads
      print("[Pipeline] Step 3: Running FastQC on cleaned reads...")
      fastqcCleaned = FastQC(
        inputFiles=[cleanedFASTQ],
        outputDir=self.fastqcDir / "cleaned" / fastqFile.stem,
        threads=2  # FastQC threading is for parallel files, 2 is optimal
      )
      fastqcCleaned.Run()

      # Step 4: STAR alignment on cleaned reads
      print("[Pipeline] Step 4: Aligning cleaned reads with STAR...")
      star = STAR(
        inputFASTQ=[cleanedFASTQ],
        referenceGenome=self.referenceGenome,
        outputDirectory=self.starDir / fastqFile.stem,
        annotationFile=self.annotationFile,
        threads=self.threads
      )
      star.run()

      # Step 5: FastQC on aligned reads
      print("[Pipeline] Step 5: Running FastQC on aligned reads...")
      # Find the BAM file from STAR output
      starOutputDir = self.starDir / fastqFile.stem / "mapped"
      bamFile = list(starOutputDir.glob("*Aligned.sortedByCoord.out.bam"))
      if bamFile:
        fastqcAligned = FastQC(
          inputFiles=[bamFile[0]],
          outputDir=self.fastqcDir / "aligned" / fastqFile.stem,
          threads=2  # FastQC threading is for parallel files, 2 is optimal
        )
        fastqcAligned.Run()
      else:
        print(f"[Pipeline] Warning: Could not find aligned BAM file in {starOutputDir}")

      print(f"[Pipeline] Completed processing: {fastqFile.name}")

if __name__ == "__main__":
  if len(sys.argv) < 4:
    print("Usage: python sweepRNAseq.py <fastq_file1> [fastq_file2 ...] [annotation_file] <reference_genome> <rrna_database>")
    print("Example: python sweepRNAseq.py sample1.fastq sample2.fastq annotation.gtf reference-genome.fa rrna-database.fa")
    sys.exit(1)

  # Parse command line arguments
  # Last 2-3 arguments are reference files, rest are FASTQ files
  fastqFiles = sys.argv[1:-2]  # All but last 2 are FASTQ files
  referenceGenome = sys.argv[-2]  # Second to last is reference genome
  rRNADatabase = sys.argv[-1]  # Last is rRNA database

  # Check if annotation file was provided (optional)
  annotationFile = None
  if len(sys.argv) > 4 and sys.argv[-3].endswith(('.gtf', '.gff')):
    annotationFile = sys.argv[-3]
    fastqFiles = sys.argv[1:-3]  # Adjust FASTQ files if annotation provided

  pipeline = Pipeline(
    inputFASTQ=fastqFiles,
    referenceGenome=referenceGenome,
    rRNADatabase=rRNADatabase,
    annotationFile=annotationFile
  )
  pipeline.run()
  print("\n[Pipeline] Pipeline completed successfully!")