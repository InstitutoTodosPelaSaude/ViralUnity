"""Configuration file generation for ViralUnity pipelines."""

import os
from typing import Any, Dict, List, Optional, Union

import yaml

from viralunity.constants import ConfigKeys, DataType, ResourceDefaults
from viralunity.exceptions import ConfigurationError


class ConfigGenerator:
    """Generates YAML configuration files for Snakemake workflows."""

    # Section names used as comment headers in the output YAML
    SECTION_PARAMETERS = "parameters"
    SECTION_DATABASES = "databases"
    SECTION_RESOURCES = "resources"
    SECTION_REFERENCE_ASSEMBLY = "reference assembly"

    def __init__(self, config_path: str):
        """Initialize config generator.

        Args:
            config_path: Path where the config file will be written
        """
        self.config_path = config_path
        self.config: Dict[str, Any] = {}
        # Track which section each key belongs to
        self._sections: Dict[str, str] = {}

    def _set(self, key: str, value: Any, section: str) -> None:
        """Set a config key and tag it to a section.

        Args:
            key: Configuration key name
            value: Configuration value
            section: Section name (parameters, databases, resources)
        """
        self.config[key] = value
        self._sections[key] = section

    def add_samples(self, samples: Dict[str, List[str]], data_type: str) -> None:
        """Add samples to configuration.

        Args:
            samples: Dictionary mapping sample names to file paths
            data_type: Type of sequencing data (illumina or nanopore)
        """
        formatted_samples = {}
        for sample_name, file_paths in samples.items():
            key = f"sample-{sample_name}"
            if data_type == DataType.ILLUMINA:
                if len(file_paths) != 2:
                    raise ConfigurationError(
                        f"Illumina sample {sample_name} must have 2 files, "
                        f"found {len(file_paths)}"
                    )
                formatted_samples[key] = f"{file_paths[0]} {file_paths[1]}"
            else:
                if len(file_paths) != 1:
                    raise ConfigurationError(
                        f"Nanopore sample {sample_name} must have 1 file, "
                        f"found {len(file_paths)}"
                    )
                formatted_samples[key] = file_paths[0]

        self._set(ConfigKeys.SAMPLES, formatted_samples, self.SECTION_PARAMETERS)
        self._set(ConfigKeys.DATA, data_type, self.SECTION_PARAMETERS)

    def add_output(self, output_dir: str, run_name: str) -> None:
        """Add output directory to configuration.

        Args:
            output_dir: Base output directory
            run_name: Name of the run
        """
        self._set(
            ConfigKeys.OUTPUT,
            os.path.join(output_dir, run_name, ""),
            self.SECTION_PARAMETERS,
        )

    def add_threads(self, threads: int) -> None:
        """Add thread count to configuration.

        Args:
            threads: Number of threads
        """
        self._set(ConfigKeys.THREADS, threads, self.SECTION_PARAMETERS)

    def add_consensus_nanopore_settings(
        self,
        minimum_read_length: int,
        af_threshold: float,
        chunk_size: int,
        clair3_model: str,
        variant_quality: int,
        variant_depth: int,
        minimum_map_quality: int,
    ) -> None:
        """Add Nanopore consensus-specific settings to configuration.

        Args:
            minimum_read_length: Minimum read length threshold
            af_threshold: Allele frequency threshold to call a variant into consensus
            chunk_size: Size of chunks to process [clair3]
            clair3_model: Model to use for variant calling [clair3]
            variant_quality: Minimum variant quality to call a variant into consensus [clair3]
            variant_depth: Minimum alt allele depth to call a variant into consensus [clair3]
            minimum_map_quality: Minimum map quality to call a variant into consensus [clair3]
        """
        P = self.SECTION_PARAMETERS
        self._set(ConfigKeys.MINIMUM_LENGTH, minimum_read_length, P)
        self._set(ConfigKeys.AF_THRESHOLD, af_threshold, P)
        self._set(ConfigKeys.CHUNK_SIZE, chunk_size, P)
        self._set(ConfigKeys.CLAIR3_MODEL, clair3_model, P)
        self._set(ConfigKeys.VARIANT_QUALITY, variant_quality, P)
        self._set(ConfigKeys.VARIANT_DEPTH, variant_depth, P)
        self._set(ConfigKeys.MINIMUM_MAP_QUALITY, minimum_map_quality, P)

    def add_illumina_settings(
        self,
        adapters: str,
        minimum_read_length: int,
        trim_head: Optional[int] = None,
        trim_tail: Optional[int] = None,
        cut_front_mean_quality: int = 20,
        cut_tail_mean_quality: int = 20,
        cut_right_window_size: int = 4,
        cut_right_mean_quality: int = 20,
        af_threshold: float = 0.5,
        af_isnv_threshold: float = 0.05,
        run_isnv: bool = False,
    ) -> None:
        """Add Illumina-specific settings to configuration (fastp QC).

        Args:
            adapters: Path to adapters file or "NA" for auto-detection
            minimum_read_length: Minimum read length threshold
            trim_head: Bases to trim from 5'
            trim_tail: Bases to trim from 3'
            cut_front_mean_quality: fastp cut_front mean quality threshold
            cut_tail_mean_quality: fastp cut_tail mean quality threshold
            cut_right_window_size: fastp cut_right window size
            cut_right_mean_quality: fastp cut_right mean quality threshold
            af_threshold: Allele frequency threshold to call a variant into consensus
            af_isnv_threshold: Minimum allele frequency threshold to call a variant into iSNV analysis
            run_isnv: Whether to run iSNV analysis
        """
        P = self.SECTION_PARAMETERS
        self._set(ConfigKeys.ADAPTERS, adapters, P)
        self._set(ConfigKeys.MINIMUM_LENGTH, minimum_read_length, P)
        self._set(ConfigKeys.TRIM_HEAD, trim_head if trim_head is not None else 0, P)
        self._set(ConfigKeys.TRIM_TAIL, trim_tail if trim_tail is not None else 0, P)
        self._set(ConfigKeys.CUT_FRONT_MEAN_QUALITY, cut_front_mean_quality, P)
        self._set(ConfigKeys.CUT_TAIL_MEAN_QUALITY, cut_tail_mean_quality, P)
        self._set(ConfigKeys.CUT_RIGHT_WINDOW_SIZE, cut_right_window_size, P)
        self._set(ConfigKeys.CUT_RIGHT_MEAN_QUALITY, cut_right_mean_quality, P)
        self._set(ConfigKeys.AF_THRESHOLD, af_threshold, P)
        self._set(ConfigKeys.AF_ISNV_THRESHOLD, af_isnv_threshold, P)
        self._set(ConfigKeys.RUN_ISNV, run_isnv, P)

    def add_consensus_settings(
        self,
        reference: Union[str, Dict[str, str]],
        primer_scheme: str,
        minimum_coverage: int,
        minimap2_consensus_align_flags: str = "-a --sam-hit-only --secondary=no --score-N=0",
    ) -> None:
        """Add consensus-specific settings to configuration.

        Args:
            reference: Path to reference genome (str) or dict mapping
            segment names to paths for segmented viruses
            primer_scheme: Path to primer scheme file or "NA"
            minimum_coverage: Minimum coverage for consensus
            minimap2_consensus_align_flags: Flags passed to minimap2 when
                re-aligning the per-sample consensus FASTA back to the
                reference for the final multiple-sequence alignment. The
                default keeps the historical behaviour.
        """
        P = self.SECTION_PARAMETERS
        self._set(ConfigKeys.REFERENCE, reference, P)
        self._set(ConfigKeys.SCHEME, primer_scheme, P)
        self._set(ConfigKeys.MINIMUM_DEPTH, minimum_coverage, P)
        self._set(
            ConfigKeys.MINIMAP2_CONSENSUS_ALIGN_FLAGS,
            minimap2_consensus_align_flags,
            P,
        )

    def add_metagenomics_settings(
        self,
        kraken2_database: str,
        krona_database: str,
        remove_human_reads: bool,
        remove_unclassified_reads: bool,
        host_reference: str = "NA",
        deacon_index: str = "NA",
        taxdump: str = "NA",
        run_denovo_assembly: bool = False,
        run_kraken2_reads: bool = True,
        run_kraken2_contigs: bool = True,
        run_diamond_reads: bool = False,
        run_diamond_contigs: bool = False,
        taxids: str = "NA",
        diamond_database: str = "NA",
        diamond_sensitivity: str = "sensitive",
        evalue: float = 0.001,
        bleed_fraction: float = 0.005,
        negative_controls: Optional[List[str]] = None,
        minimum_hit_group: int = 4,
        diamond_max_target_seqs: int = 1,
        kraken2_extra_flags: str = "--report-minimizer-data",
        compute_rpkm: bool = False,
        enrichment_pseudocount: float = 1.0,
        z_score_threshold: float = 3.0,
        log2_ratio_threshold: float = 1.0,
    ) -> None:
        """Add metagenomics-specific settings to configuration.

        Args:
            kraken2_database: Path to Kraken2 database
            krona_database: Path to Krona database
            remove_human_reads: Whether to remove human reads
            remove_unclassified_reads: Whether to remove unclassified reads
            host_reference: Path to host genome FASTA for dehosting (or "NA")
            deacon_index: Path to Deacon minimizer index for host depletion (or "NA"). If set, used instead of host_reference for dehosting.
            taxdump: Path to NCBI taxdump dir (nodes.dmp, names.dmp)
            run_denovo_assembly: Whether to run MEGAHIT assembly
            run_kraken2_reads: Whether to run Kraken2 on reads
            run_kraken2_contigs: Whether to run Kraken2 on contigs (if assembly)
            run_diamond_reads: Whether to run DIAMOND on reads
            run_diamond_contigs: Whether to run DIAMOND on contigs (if assembly)
            taxids: NCBI assembly summary (protein2taxid) for Diamond taxonomy (or "NA")
            diamond_database: Path to Diamond DB (protein FASTA)
            diamond_sensitivity: Diamond sensitivity (e.g. sensitive, mid-sensitive)
            evalue: E-value threshold for Diamond
            bleed_fraction: Max-RPM bleed filter fraction
            negative_controls: Sample IDs to use as negative controls
            minimum_hit_group: Kraken2 --minimum-hit-group (default: 4)
            diamond_max_target_seqs: DIAMOND --max-target-seqs value (default 1).
            kraken2_extra_flags: Extra flags appended to every kraken2
                invocation alongside ``--threads`` and ``--minimum-hit-group``.
                Defaults to ``--report-minimizer-data``; set to ``""`` to drop it.
            compute_rpkm: Whether to compute RPKM (requires viral_genomes and
                viral_taxids to be set).  Derived automatically from
                ``viral_genomes != "NA"`` in the CLI layer.
            enrichment_pseudocount: Pseudocount added to both numerator and
                denominator when computing fold-enrichment and log2-ratio.
            z_score_threshold: Minimum z-score to consider a hit above
                background (used when ≥2 negative controls are present).
            log2_ratio_threshold: Minimum log2-ratio to consider a hit above
                background (used when exactly 1 negative control is present, or
                when z-score is undefined due to zero control variance).
        """
        P = self.SECTION_PARAMETERS
        D = self.SECTION_DATABASES
        # Database paths
        self._set(ConfigKeys.KRAKEN2_DATABASE, kraken2_database, D)
        self._set(ConfigKeys.KRONA_DATABASE, krona_database, D)
        self._set(ConfigKeys.HOST_REFERENCE, host_reference, D)
        self._set(ConfigKeys.DEACON_INDEX, deacon_index, D)
        self._set(ConfigKeys.TAXDUMP, taxdump, D)
        self._set(ConfigKeys.TAXIDS, taxids, D)
        self._set(ConfigKeys.DIAMOND_DATABASE, diamond_database, D)
        # Pipeline parameters
        self._set(ConfigKeys.REMOVE_HUMAN_READS, remove_human_reads, P)
        self._set(ConfigKeys.REMOVE_UNCLASSIFIED_READS, remove_unclassified_reads, P)
        self._set(ConfigKeys.RUN_DENOVO_ASSEMBLY, run_denovo_assembly, P)
        self._set(ConfigKeys.RUN_KRAKEN2_READS, run_kraken2_reads, P)
        self._set(ConfigKeys.RUN_KRAKEN2_CONTIGS, run_kraken2_contigs, P)
        self._set(ConfigKeys.RUN_DIAMOND_READS, run_diamond_reads, P)
        self._set(ConfigKeys.RUN_DIAMOND_CONTIGS, run_diamond_contigs, P)
        self._set(ConfigKeys.DIAMOND_SENSITIVITY, diamond_sensitivity, P)
        self._set(ConfigKeys.EVALUE, evalue, P)
        self._set(ConfigKeys.BLEED_FRACTION, bleed_fraction, P)
        self._set(ConfigKeys.NEGATIVE_CONTROLS, negative_controls or [], P)
        self._set(ConfigKeys.MINIMUM_HIT_GROUP, minimum_hit_group, P)
        self._set(ConfigKeys.DIAMOND_MAX_TARGET_SEQS, diamond_max_target_seqs, P)
        self._set(ConfigKeys.KRAKEN2_EXTRA_FLAGS, kraken2_extra_flags, P)
        self._set(ConfigKeys.COMPUTE_RPKM, compute_rpkm, P)
        self._set(ConfigKeys.ENRICHMENT_PSEUDOCOUNT, enrichment_pseudocount, P)
        self._set(ConfigKeys.Z_SCORE_THRESHOLD, z_score_threshold, P)
        self._set(ConfigKeys.LOG2_RATIO_THRESHOLD, log2_ratio_threshold, P)

    def add_reference_assembly_settings(
        self,
        run_reference_assembly: bool = False,
        method: Optional[str] = None,
        source: Optional[str] = None,
        reads_count: int = 100,
        contigs_count: int = 1,
        families: str = "Coronaviridae,Orthomyxoviridae,Flaviviridae,Herpesviridae,Papillomaviridae,Paramyxoviridae,Adenoviridae",
        reference_selection_strategy: str = "taxid",
        blast_qcov: int = 80,
        blast_pident: int = 80,
        viral_genomes: str = "databases/virus_genomes/viral.genomes.fasta",
        viral_taxids: str = "databases/virus_genomes/genome2taxid.tsv",
    ) -> None:
        """Add settings for reference assembly on metagenomics hits."""
        section = self.SECTION_REFERENCE_ASSEMBLY
        self._set(ConfigKeys.RUN_REFERENCE_ASSEMBLY, run_reference_assembly, section)
        self._set(ConfigKeys.REF_ASSEMBLY_METHOD, method, section)
        self._set(ConfigKeys.REF_ASSEMBLY_SOURCE, source, section)
        self._set(ConfigKeys.REF_ASSEMBLY_READS_COUNT, reads_count, section)
        self._set(ConfigKeys.REF_ASSEMBLY_CONTIGS_COUNT, contigs_count, section)
        self._set(
            ConfigKeys.REF_ASSEMBLY_FAMILIES,
            [f.strip() for f in families.split(",")],
            section,
        )
        self._set(ConfigKeys.REF_SELECTION_STRATEGY, reference_selection_strategy, section)
        self._set(ConfigKeys.REF_BLAST_QCOV, blast_qcov, section)
        self._set(ConfigKeys.REF_BLAST_PIDENT, blast_pident, section)

        db_section = self.SECTION_DATABASES
        self._set(ConfigKeys.VIRAL_GENOMES, viral_genomes, db_section)
        self._set(ConfigKeys.VIRAL_TAXIDS, viral_taxids, db_section)

    def add_nanopore_settings(
        self,
        run_polish_racon: bool = False,
        run_polish_medaka: bool = False,
        medaka_model: Optional[str] = None,
        clair3_model: Optional[str] = None,
    ) -> None:
        """Add Nanopore-specific settings (polishing: Racon, Medaka; clair3 model).

        Args:
            run_polish_racon: Whether to run Racon polishing after MEGAHIT.
            run_polish_medaka: Whether to run Medaka polishing (after Racon if both enabled).
            medaka_model: Medaka model name (e.g. r941_min_high_g360). Optional; Medaka uses default if not set.
            clair3_model: Clair3 model name used by the reference-assembly
                consensus path. Optional; the rule falls back to its own
                default if not set.
        """
        P = self.SECTION_PARAMETERS
        self._set(ConfigKeys.RUN_POLISH_RACON, run_polish_racon, P)
        self._set(ConfigKeys.RUN_POLISH_MEDAKA, run_polish_medaka, P)
        if medaka_model is not None:
            self._set(ConfigKeys.MEDAKA_MODEL, medaka_model, P)
        if clair3_model is not None:
            self._set(ConfigKeys.CLAIR3_MODEL, clair3_model, P)

    def add_workflow_path(self, workflow_path: str) -> None:
        """Add workflow path to configuration.

        Args:
            workflow_path: Path to the workflow directory
        """
        self._set("workflow_path", workflow_path, self.SECTION_PARAMETERS)

    def add_resource_settings(self, args: Dict[str, Any], rule_names: list) -> None:
        """Add per-rule resource settings (CPUs and RAM) to configuration.

        For each rule name in rule_names, writes ``<rule>_cpus`` and
        ``<rule>_ram`` keys to the config dict.  Values are taken from
        *args* if present; otherwise the defaults from
        ``ResourceDefaults`` are used.

        Args:
            args: Dictionary of pipeline arguments (from the CLI).
            rule_names: List of Snakemake rule name strings that should
                receive resource entries.
        """
        R = self.SECTION_RESOURCES
        for rule in rule_names:
            cpus_key = f"{rule}_cpus"
            ram_key = f"{rule}_ram"
            self._set(cpus_key, args.get(cpus_key, ResourceDefaults.DEFAULT_CPUS), R)
            self._set(ram_key, args.get(ram_key, ResourceDefaults.DEFAULT_RAM), R)

    @classmethod
    def write_skeleton(
        cls,
        pipeline: str,
        data_type: str,
        config_path: str,
        placeholder_dir: str,
    ) -> str:
        """Write a placeholder YAML config that lets Snakemake parse a workflow.

        Used by ``viralunity setup`` to materialize per-rule conda envs via
        ``snakemake(..., conda_create_envs_only=True)``. Even in
        envs-only mode Snakemake walks the DAG and verifies that rule
        inputs exist, so the caller is expected to populate
        ``placeholder_dir`` with empty files for the FASTQs, references,
        and database paths emitted here. ``viralunity_setup_cli`` does
        this for the package; tests use the same convention.

        Args:
            pipeline: ``"consensus"`` or ``"metagenomics"``.
            data_type: ``"illumina"`` or ``"nanopore"``.
            config_path: Where to write the YAML.
            placeholder_dir: Directory that contains the placeholder input
                files referenced by the generated config.

        Returns:
            ``config_path``.
        """
        gen = cls(config_path)
        root = placeholder_dir.rstrip("/")

        if data_type == DataType.ILLUMINA:
            placeholder_samples = {
                "skel": [
                    f"{root}/reads/skel_R1.fastq.gz",
                    f"{root}/reads/skel_R2.fastq.gz",
                ],
            }
        else:
            placeholder_samples = {"skel": [f"{root}/reads/skel.fastq.gz"]}

        gen.add_samples(placeholder_samples, data_type)
        gen.add_output(f"{root}/output", "skeleton_run")
        gen.add_threads(1)

        # Every optional feature flag is enabled below so Snakemake's DAG
        # walk reaches every rule and ``conda_create_envs_only=True``
        # materializes every per-rule env that pipeline could possibly need.
        # If a flag is left off here, its env is silently skipped by
        # ``viralunity setup`` and falls back to dynamic env creation at
        # the user's first real run — exactly the failure mode setup
        # exists to prevent.
        if pipeline == "consensus":
            gen.add_consensus_settings(
                reference=f"{root}/references/skel.reference.fasta",
                primer_scheme="NA",
                minimum_coverage=20,
            )
            if data_type == DataType.ILLUMINA:
                # run_isnv=True pulls the detect_isnv (LoFreq) rule into
                # the DAG, which uses envs/consensus.yaml.
                gen.add_illumina_settings(
                    adapters="NA",
                    minimum_read_length=50,
                    run_isnv=True,
                )
            else:
                gen.add_consensus_nanopore_settings(
                    minimum_read_length=50,
                    af_threshold=0.51,
                    chunk_size=10000,
                    clair3_model="r1041_e82_400bps_sup_v500",
                    variant_quality=20,
                    variant_depth=10,
                    minimum_map_quality=30,
                )
            gen.add_workflow_path(".")
        elif pipeline == "metagenomics":
            if data_type == DataType.ILLUMINA:
                gen.add_illumina_settings(adapters="NA", minimum_read_length=50)
            gen.add_metagenomics_settings(
                kraken2_database=f"{root}/kraken2/",
                krona_database=f"{root}/krona/taxonomy/",
                remove_human_reads=False,
                remove_unclassified_reads=False,
                taxdump=f"{root}/taxdump/",
                taxids=f"{root}/diamond/protein2taxid.tsv",
                diamond_database=f"{root}/diamond/nr.dmnd",
                # All four classifier toggles on => kraken2 + diamond
                # rules on reads *and* contigs are in the DAG; pulls
                # taxonomy.yaml.
                run_denovo_assembly=True,
                run_kraken2_reads=True,
                run_kraken2_contigs=True,
                run_diamond_reads=True,
                run_diamond_contigs=True,
            )
            # Reference-assembly checkpoint pulls genome_selection.yaml.
            gen.add_reference_assembly_settings(
                run_reference_assembly=True,
                method="kraken2",
                source="reads",
                viral_genomes=f"{root}/virus_genomes/viral.genomes.fasta",
                viral_taxids=f"{root}/virus_genomes/genome2taxid.tsv",
            )
            if data_type == DataType.NANOPORE:
                # Polish flags pull racon (consensus.yaml) and medaka.yaml.
                gen.add_nanopore_settings(
                    run_polish_racon=True,
                    run_polish_medaka=True,
                )
        else:
            raise ValueError(
                f"Unknown pipeline: {pipeline!r} (expected 'consensus' or 'metagenomics')"
            )

        gen.save()
        return config_path

    # Empty files (relative to ``placeholder_dir``) that
    # ``ConfigGenerator.write_skeleton`` expects to exist on disk for
    # ``snakemake --conda-create-envs-only`` to succeed. Kept here so
    # the skeleton and its placeholder-file expectations stay in lock-step.
    SKELETON_PLACEHOLDERS: Dict[str, Dict[str, List[str]]] = {
        "consensus": {
            "illumina": [
                "reads/skel_R1.fastq.gz",
                "reads/skel_R2.fastq.gz",
                "references/skel.reference.fasta",
            ],
            "nanopore": [
                "reads/skel.fastq.gz",
                "references/skel.reference.fasta",
            ],
        },
        "metagenomics": {
            "illumina": [
                "reads/skel_R1.fastq.gz",
                "reads/skel_R2.fastq.gz",
                "kraken2/hash.k2d",
                "krona/taxonomy/taxonomy.tab",
                "taxdump/nodes.dmp",
                "taxdump/names.dmp",
                "diamond/nr.dmnd",
                "diamond/protein2taxid.tsv",
                "virus_genomes/viral.genomes.fasta",
                "virus_genomes/genome2taxid.tsv",
            ],
            "nanopore": [
                "reads/skel.fastq.gz",
                "kraken2/hash.k2d",
                "krona/taxonomy/taxonomy.tab",
                "taxdump/nodes.dmp",
                "taxdump/names.dmp",
                "diamond/nr.dmnd",
                "diamond/protein2taxid.tsv",
                "virus_genomes/viral.genomes.fasta",
                "virus_genomes/genome2taxid.tsv",
            ],
        },
    }

    def save(self) -> None:
        """Save configuration to YAML file with section comment headers.

        Keys are grouped into sections (# parameters, # databases,
        # resources) based on the tags assigned by ``_set()``.  Within
        each section the insertion order of keys is preserved.

        Raises:
            ConfigurationError: If config directory cannot be created
        """
        config_dir = os.path.dirname(self.config_path)
        if config_dir:  # Only create directory if path contains a directory component
            os.makedirs(config_dir, exist_ok=True)

        # Group keys by section, preserving insertion order
        section_order = [
            self.SECTION_PARAMETERS,
            self.SECTION_REFERENCE_ASSEMBLY,
            self.SECTION_DATABASES,
            self.SECTION_RESOURCES,
        ]
        grouped: Dict[str, Dict[str, Any]] = {s: {} for s in section_order}

        for key, value in self.config.items():
            section = self._sections.get(key, self.SECTION_PARAMETERS)
            grouped[section][key] = value

        try:
            with open(self.config_path, "w") as f:
                first = True
                for section in section_order:
                    items = grouped[section]
                    if not items:
                        continue
                    if not first:
                        f.write("\n")
                    f.write(f"# {section}\n")
                    yaml.dump(
                        items,
                        f,
                        default_flow_style=False,
                        sort_keys=False,
                    )
                    first = False
        except (OSError, IOError) as e:
            raise ConfigurationError(
                f"Failed to write config file to {self.config_path}: {e}"
            ) from e
