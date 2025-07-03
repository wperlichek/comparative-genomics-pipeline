import asyncio
import logging
import sys
import subprocess
from .client import UniProtClient, NCBIClient, EBIClient, pdp_client
from .client.clinvar_client import ClinVarClient
from .config import logging_config, path_config
from .util import file_util
from .service import biopython_service
from pathlib import Path

logger = logging.getLogger(__name__)


async def collect_orthologous_sequences(uni_prot_client, ncbi_client) -> bool:
    """Collect orthologous sequences from UniProt and NCBI APIs.
    
    Returns:
        bool: True if successful, False if critical errors occurred
    """
    try:
        genes_to_proteins = file_util.open_file_return_as_json(
            f"{path_config.DATA_INPUT_DIR}/genes_to_proteins.json"
        )
        if not genes_to_proteins:
            logger.error("No genes found in configuration file")
            return False
        
        # Validate configuration structure
        if not file_util.validate_genes_config(genes_to_proteins):
            logger.error("Configuration validation failed")
            return False
            
        logger.info(f"Processing {len(genes_to_proteins)} genes for ortholog collection")
        
        for gene_name, ortholog_list in genes_to_proteins.items():
            logger.info(f"Collecting orthologs for {gene_name}")
            all_orthologs_as_fasta = ""
            successful_retrievals = 0
            
            for i, ortholog in enumerate(ortholog_list):
                fasta = ""
                species = ortholog.get("species", "unknown")
                
                try:
                    if "uniprot_id" in ortholog and ortholog["uniprot_id"] != "":
                        logger.debug(f"Fetching UniProt sequence for {species}: {ortholog['uniprot_id']}")
                        fasta = await uni_prot_client.fetch_protein_fasta_sequence_by_accession_id(
                            ortholog["uniprot_id"]
                        )
                    
                    if not fasta and "entrez_protein_id" in ortholog and ortholog["entrez_protein_id"] != "":
                        logger.debug(f"Fetching NCBI sequence for {species}: {ortholog['entrez_protein_id']}")
                        fasta = await ncbi_client.fetch_protein_fasta_by_entrez_id(
                            ortholog["entrez_protein_id"]
                        )
                    
                    if fasta:
                        all_orthologs_as_fasta += fasta
                        successful_retrievals += 1
                        logger.debug(f"Successfully retrieved sequence for {species}")
                    else:
                        logger.warning(f"No sequence retrieved for {species} (ortholog {i+1})")
                        
                except Exception as e:
                    logger.error(f"Failed to retrieve sequence for {species}: {e}")
                    continue
            
            if successful_retrievals == 0:
                logger.error(f"No sequences retrieved for {gene_name}. Skipping.")
                continue
            elif successful_retrievals < len(ortholog_list):
                logger.warning(f"Only {successful_retrievals}/{len(ortholog_list)} sequences retrieved for {gene_name}")
            
            # Save sequences to file
            try:
                orthologs_dir = path_config.ORTHOLOGS_OUTPUT_DIR
                orthologs_dir.mkdir(parents=True, exist_ok=True)
                file_path = orthologs_dir / f"{gene_name}.fasta"
                
                with open(file_path, "w") as fasta_file:
                    fasta_file.write(all_orthologs_as_fasta)
                
                logger.info(f"Saved {successful_retrievals} sequences for {gene_name} to {file_path}")
                
            except OSError as e:
                logger.error(f"Failed to save sequences for {gene_name}: {e}")
                return False
        
        return True
        
    except Exception as e:
        logger.error(f"Critical error in ortholog collection: {e}")
        return False


async def align_sequences_msa(ebi_client) -> bool:
    """Perform multiple sequence alignment using EBI Clustal Omega.
    
    Returns:
        bool: True if at least one alignment succeeded, False if all failed
    """
    try:
        orthologs_dir = path_config.ORTHOLOGS_OUTPUT_DIR
        msa_dir = path_config.MSA_OUTPUT_DIR
        msa_dir.mkdir(parents=True, exist_ok=True)
        
        fasta_files = list(orthologs_dir.glob("*.fasta"))
        if not fasta_files:
            logger.error("No FASTA files found for alignment")
            return False
        
        logger.info(f"Starting MSA for {len(fasta_files)} files")
        successful_alignments = 0
        
        for file_path in fasta_files:
            try:
                gene_name = file_path.stem
                logger.info(f"Aligning sequences for {gene_name}")
                
                ortholog_fasta = file_util.open_file_return_as_str(file_path)
                if not ortholog_fasta or not ortholog_fasta.strip():
                    logger.error(f"Empty or unreadable FASTA file for {gene_name}")
                    continue
                
                job_id = await ebi_client.submit_job(ortholog_fasta)
                if not job_id:
                    logger.error(f"Failed to submit MSA job for {gene_name}")
                    continue
                
                logger.debug(f"MSA job submitted for {gene_name}: {job_id}")
                
                # Poll for job completion with timeout
                max_polls = 60  # 5 minutes timeout (5 seconds * 60)
                poll_count = 0
                status = ""
                
                while status not in ("FINISHED", "ERROR", "FAILURE") and poll_count < max_polls:
                    try:
                        status = await ebi_client.check_status(job_id)
                        logger.debug(f"MSA job {job_id} status: {status} (poll {poll_count+1}/{max_polls})")
                        
                        if status == "FINISHED":
                            result = await ebi_client.get_result(job_id, "fa")
                            if result:
                                msa_path = msa_dir / f"{gene_name}.fasta"
                                with open(msa_path, "w") as msa_file:
                                    msa_file.write(result)
                                logger.info(f"MSA completed for {gene_name}: {msa_path}")
                                successful_alignments += 1
                            else:
                                logger.error(f"Empty result returned for MSA job {job_id}")
                            break
                            
                        elif status in ("ERROR", "FAILURE"):
                            logger.error(f"MSA job {job_id} failed with status: {status}")
                            break
                            
                        else:
                            await asyncio.sleep(5)
                            poll_count += 1
                            
                    except Exception as e:
                        logger.error(f"Error polling MSA job {job_id}: {e}")
                        break
                
                if poll_count >= max_polls:
                    logger.error(f"MSA job {job_id} timed out after {max_polls * 5} seconds")
                    
            except Exception as e:
                logger.error(f"Failed to process MSA for {gene_name}: {e}")
                continue
        
        if successful_alignments == 0:
            logger.error("No MSA alignments completed successfully")
            return False
        
        logger.info(f"MSA completed: {successful_alignments}/{len(fasta_files)} alignments successful")
        return True
        
    except Exception as e:
        logger.error(f"Critical error in MSA alignment: {e}")
        return False


async def generate_phylogenetic_trees(ebi_client) -> bool:
    """Generate phylogenetic trees using EBI services.
    
    Returns:
        bool: True if at least one tree generated successfully, False if all failed
    """
    try:
        orthologs_dir = path_config.ORTHOLOGS_OUTPUT_DIR
        trees_dir = path_config.TREES_OUTPUT_DIR
        trees_dir.mkdir(parents=True, exist_ok=True)
        
        fasta_files = list(orthologs_dir.glob("*.fasta"))
        if not fasta_files:
            logger.error("No FASTA files found for phylogenetic tree generation")
            return False
        
        logger.info(f"Starting phylogenetic tree generation for {len(fasta_files)} files")
        successful_trees = 0
        
        for file_path in fasta_files:
            try:
                gene_name = file_path.stem
                logger.info(f"Generating phylogenetic tree for {gene_name}")
                
                with open(file_path, "r") as f:
                    fasta_str = f.read()
                
                if not fasta_str.strip():
                    logger.error(f"Empty FASTA file for {gene_name}")
                    continue
                
                job_id = await ebi_client.submit_job(fasta_str)
                if not job_id:
                    logger.error(f"Failed to submit tree generation job for {gene_name}")
                    continue
                
                logger.debug(f"Tree generation job submitted for {gene_name}: {job_id}")
                
                # Wait for job with timeout and polling
                max_polls = 120  # 10 minutes timeout (5 seconds * 120)
                poll_count = 0
                initial_wait = 15  # Initial wait before first poll
                
                await asyncio.sleep(initial_wait)
                
                while poll_count < max_polls:
                    try:
                        status = await ebi_client.check_status(job_id)
                        logger.debug(f"Tree job {job_id} status: {status} (poll {poll_count+1}/{max_polls})")
                        
                        if status == "FINISHED":
                            tree = await ebi_client.get_phylogenetic_tree(job_id)
                            if tree and tree.strip():
                                tree_path = trees_dir / f"{gene_name}.nwk"
                                with open(tree_path, "w") as tree_file:
                                    tree_file.write(tree)
                                logger.info(f"Phylogenetic tree saved for {gene_name}: {tree_path}")
                                successful_trees += 1
                            else:
                                logger.error(f"Empty or invalid tree returned for {gene_name}")
                            break
                            
                        elif status in ("ERROR", "FAILURE"):
                            logger.error(f"Tree generation job {job_id} failed with status: {status}")
                            break
                            
                        elif status in ("RUNNING", "PENDING", "QUEUED"):
                            await asyncio.sleep(5)
                            poll_count += 1
                            
                        else:
                            logger.warning(f"Unknown job status for {job_id}: {status}")
                            await asyncio.sleep(5)
                            poll_count += 1
                            
                    except Exception as e:
                        logger.error(f"Error polling tree generation job {job_id}: {e}")
                        break
                
                if poll_count >= max_polls:
                    logger.error(f"Tree generation job {job_id} timed out after {max_polls * 5 + initial_wait} seconds")
                    
            except Exception as e:
                logger.error(f"Failed to process tree generation for {gene_name}: {e}")
                continue
        
        if successful_trees == 0:
            logger.error("No phylogenetic trees generated successfully")
            return False
        
        logger.info(f"Tree generation completed: {successful_trees}/{len(fasta_files)} trees successful")
        return True
        
    except Exception as e:
        logger.error(f"Critical error in phylogenetic tree generation: {e}")
        return False


async def async_main() -> int:
    """Main async pipeline execution.
    
    Returns:
        int: Exit code (0 for success, 1 for failure)
    """
    uni_prot_client = None
    ncbi_client = None
    ebi_client = None
    pdb_client = None
    
    try:
        # Initialize clients
        logger.info("Initializing API clients...")
        uni_prot_client = UniProtClient()
        ncbi_client = NCBIClient()
        ebi_client = EBIClient()
        pdb_client = pdp_client.PDBClient()
        
        # Step 1: Collect orthologous sequences
        logger.info("Step 1: Collecting orthologous sequences...")
        if not await collect_orthologous_sequences(uni_prot_client, ncbi_client):
            logger.error("Failed to collect orthologous sequences. Aborting pipeline.")
            return 1
        
        # Step 2: Perform multiple sequence alignment
        logger.info("Step 2: Performing multiple sequence alignment...")
        if not await align_sequences_msa(ebi_client):
            logger.error("Failed to perform MSA. Aborting pipeline.")
            return 1
        
        # Step 3: Generate phylogenetic trees
        logger.info("Step 3: Generating phylogenetic trees...")
        if not await generate_phylogenetic_trees(ebi_client):
            logger.error("Failed to generate phylogenetic trees. Aborting pipeline.")
            return 1
        
        # Step 4: Fetch ClinVar variants
        logger.info("Step 4: Fetching ClinVar variants...")
        clinvar_client = ClinVarClient()
        genes_to_proteins = file_util.open_file_return_as_json(
            f"{path_config.DATA_INPUT_DIR}/genes_to_proteins.json"
        )
        
        if genes_to_proteins:
            for gene_name in genes_to_proteins.keys():
                try:
                    variants = await clinvar_client.get_gene_variants(gene_name, limit=100)
                    if variants:
                        # Save variants to CSV
                        import pandas as pd
                        df = pd.DataFrame(variants)
                        variants_dir = path_config.VARIANTS_OUTPUT_DIR
                        variants_dir.mkdir(parents=True, exist_ok=True)
                        df.to_csv(variants_dir / f"{gene_name}_clinvar_variants.csv", index=False)
                        logger.info(f"Saved {len(variants)} ClinVar variants for {gene_name}")
                except Exception as e:
                    logger.error(f"Failed to fetch ClinVar variants for {gene_name}: {e}")
        
        # Step 5: Fetch UniProt variants
        logger.info("Step 5: Fetching UniProt variants...")
        try:
            variant_fetch_errors = 0
            for gene_name, ortholog_list in genes_to_proteins.items():
                try:
                    # Find canonical human sequence
                    canonical = next((
                        o for o in ortholog_list 
                        if o.get("species", "").lower() in ["homo sapiens", "human"] 
                        and o.get("uniprot_id")
                    ), None)
                    
                    if not canonical:
                        canonical = next((o for o in ortholog_list if o.get("uniprot_id")), None)
                    
                    if canonical and canonical.get("uniprot_id"):
                        accession = canonical["uniprot_id"]
                        output_dir = str(path_config.VARIANTS_OUTPUT_DIR)
                        logger.info(f"Fetching variants for {gene_name} ({accession})...")
                        
                        await uni_prot_client.fetch_protein_variants_by_accession_id(accession, output_dir, gene_name)
                        
                    else:
                        logger.warning(f"No UniProt ID found for {gene_name}. Skipping variant fetch.")
                        
                except Exception as e:
                    logger.error(f"Failed to fetch variants for {gene_name}: {e}")
                    variant_fetch_errors += 1
                    continue
            
            if variant_fetch_errors > 0:
                logger.warning(f"Variant fetch completed with {variant_fetch_errors} errors")
                
        except Exception as e:
            logger.error(f"Critical error in variant fetching: {e}")
            return 1
        
        # Step 5: Compute conservation scores
        logger.info("Step 5: Computing conservation scores...")
        try:
            biopython_service.compute_conservation_for_all_msas()
        except Exception as e:
            logger.error(f"Failed to compute conservation scores: {e}")
            return 1
        
        # Step 6: Generate visualizations
        logger.info("Step 6: Generating visualizations...")
        try:
            logger.info("Generating traditional plots...")
            biopython_service.plot_all_conservation_scores()
            biopython_service.visualize_and_save_trees()
            
            logger.info("Generating scientific publication-quality plots...")
            biopython_service.plot_all_conservation_scientific()
            biopython_service.visualize_trees_scientific()
            
        except Exception as e:
            logger.error(f"Failed to generate basic visualizations: {e}")
            return 1
        
        # Step 7: Generate variant plots
        logger.info("Step 7: Generating variant analysis plots...")
        try:
            genes_to_proteins = file_util.open_file_return_as_json(
                f"{path_config.DATA_INPUT_DIR}/genes_to_proteins.json"
            )
            
            if not genes_to_proteins:
                logger.error("Failed to load genes configuration for variant plotting")
                return 1
            
            variant_plot_errors = 0
            for gene_name, ortholog_list in genes_to_proteins.items():
                try:
                    canonical = next((
                        o for o in ortholog_list 
                        if o.get("species", "").lower() in ["homo sapiens", "human"] 
                        and o.get("uniprot_id")
                    ), None)
                    
                    if not canonical:
                        canonical = next((o for o in ortholog_list if o.get("uniprot_id")), None)
                    
                    if canonical and canonical.get("uniprot_id"):
                        accession = canonical["uniprot_id"]
                        conservation_csv = str(Path(path_config.CONSERVATION_OUTPUT_DIR) / f"{gene_name}_conservation.csv")
                        variants_csv = str(Path(path_config.VARIANTS_OUTPUT_DIR) / f"{gene_name}_{accession}_variants.csv")
                        output_dir = str(path_config.VARIANTS_OUTPUT_DIR)
                        
                        # Check if required files exist
                        if not Path(conservation_csv).exists():
                            logger.warning(f"Conservation file not found for {gene_name}: {conservation_csv}")
                            continue
                        
                        if not Path(variants_csv).exists():
                            logger.warning(f"Variants file not found for {gene_name}: {variants_csv}")
                            continue
                        
                        logger.info(f"Generating variant plots for {gene_name}...")
                        biopython_service.plot_variants_on_conservation(conservation_csv, variants_csv, output_dir)
                        biopython_service.plot_variants_scientific(conservation_csv, variants_csv, output_dir)
                        
                    else:
                        logger.warning(f"No UniProt ID found for {gene_name}. Skipping variant plots.")
                        
                except Exception as e:
                    logger.error(f"Failed to generate variant plots for {gene_name}: {e}")
                    variant_plot_errors += 1
                    continue
            
            if variant_plot_errors > 0:
                logger.warning(f"Variant plot generation completed with {variant_plot_errors} errors")
                
        except Exception as e:
            logger.error(f"Critical error in variant plot generation: {e}")
            return 1
        
        # Step 8: Generate ClinVar comparison plot
        logger.info("Step 8: Generating ClinVar comparison plot...")
        try:
            from .visualization.scientific_plots import ClinVarPlotter
            clinvar_plotter = ClinVarPlotter()
            
            scn1a_clinvar = path_config.VARIANTS_OUTPUT_DIR / "SCN1A_clinvar_variants.csv"
            depdc5_clinvar = path_config.VARIANTS_OUTPUT_DIR / "DEPDC5_clinvar_variants.csv"
            
            if scn1a_clinvar.exists() or depdc5_clinvar.exists():
                clinvar_plotter.plot_clinvar_variants(scn1a_clinvar, depdc5_clinvar, path_config.VARIANTS_OUTPUT_DIR)
                logger.info("Generated ClinVar variants comparison plot")
            else:
                logger.warning("No ClinVar variant files found for plotting")
                
        except Exception as e:
            logger.error(f"Failed to generate ClinVar comparison plot: {e}")
            return 1
        
        # Step 9: Generate protein domain visualizations
        logger.info("Step 9: Generating protein domain visualizations...")
        try:
            from .visualization.domain_visualization import visualize_protein_domains
            from .visualization.conservation_domain_plot import create_conservation_domain_plot
            
            genes_to_proteins = file_util.open_file_return_as_json(
                f"{path_config.DATA_INPUT_DIR}/genes_to_proteins.json"
            )
            
            if genes_to_proteins:
                for gene_name, ortholog_list in genes_to_proteins.items():
                    # Focus on SCN1A P35498 as requested
                    if gene_name == "SCN1A":
                        canonical = next((
                            o for o in ortholog_list 
                            if o.get("species", "").lower() in ["homo sapiens", "human"] 
                            and o.get("uniprot_id") == "P35498"
                        ), None)
                        
                        if canonical and canonical.get("uniprot_id"):
                            accession = canonical["uniprot_id"]
                            logger.info(f"Fetching domain data for {gene_name} {accession}")
                            
                            # Fetch domain data from UniProt
                            domain_data = pdb_client.fetch_uniprot_domains(accession)
                            if domain_data:
                                # Save domain data
                                pdb_client.save_domain_data(domain_data)
                                
                                # Generate standard domain visualization
                                plot_path = visualize_protein_domains(domain_data)
                                if plot_path:
                                    logger.info(f"Generated domain visualization: {plot_path}")
                                else:
                                    logger.warning(f"Failed to generate domain visualization for {accession}")
                                
                                # Generate conservation-domain overlay plot
                                try:
                                    conservation_plot_path = create_conservation_domain_plot(gene_name, accession)
                                    logger.info(f"Generated conservation-domain plot: {conservation_plot_path}")
                                except Exception as e:
                                    logger.warning(f"Failed to generate conservation-domain plot for {gene_name}: {e}")
                                    
                            else:
                                logger.warning(f"Failed to fetch domain data for {accession}")
                        else:
                            logger.warning(f"P35498 not found in {gene_name} orthologs")
            else:
                logger.warning("No genes configuration found for domain visualization")
                
        except Exception as e:
            logger.error(f"Failed to generate domain visualizations: {e}")
            # Don't return 1 here - domain visualization is not critical
        
        logger.info("Pipeline completed successfully!")
        return 0
        
    except Exception as e:
        logger.error(f"Critical pipeline error: {e}")
        return 1
        
    finally:
        # Clean up async clients
        logger.info("Cleaning up resources...")
        cleanup_errors = []
        
        try:
            if ebi_client:
                await ebi_client.close()
        except Exception as e:
            cleanup_errors.append(f"EBI client cleanup failed: {e}")
        
        try:
            if uni_prot_client:
                await uni_prot_client.close()
        except Exception as e:
            cleanup_errors.append(f"UniProt client cleanup failed: {e}")
        
        try:
            if ncbi_client:
                await ncbi_client.close()
        except Exception as e:
            cleanup_errors.append(f"NCBI client cleanup failed: {e}")
        
        try:
            if pdb_client:
                pdb_client.close()
        except Exception as e:
            cleanup_errors.append(f"PDB client cleanup failed: {e}")
        
        if cleanup_errors:
            for error in cleanup_errors:
                logger.warning(error)


def run_tests() -> bool:
    """Run tests before starting the pipeline.
    
    Returns:
        bool: True if all tests pass, False otherwise
    """
    try:
        print("Running tests before starting pipeline...")
        result = subprocess.run(
            ["python", "-m", "pytest", "-v", "--tb=short"],
            capture_output=True,
            text=True,
            timeout=120  # 2 minute timeout
        )
        
        if result.returncode == 0:
            print("✅ All tests passed!")
            return True
        else:
            print("❌ Tests failed!")
            print("STDOUT:", result.stdout[-1000:])  # Last 1000 chars
            print("STDERR:", result.stderr[-1000:])  # Last 1000 chars
            return False
            
    except subprocess.TimeoutExpired:
        print("❌ Tests timed out after 2 minutes")
        return False
    except Exception as e:
        print(f"❌ Error running tests: {e}")
        return False


def main() -> int:
    """Main entry point for the comparative genomics pipeline.
    
    Returns:
        int: Exit code (0 for success, 1 for failure)
    """
    try:
        logging_config.setup_logging()
        
        # Run tests first
        if not run_tests():
            logger.error("Tests failed - aborting pipeline execution")
            print("\n❌ Pipeline aborted due to test failures")
            print("Fix the failing tests before running the pipeline")
            return 1
        
        logger.info("Starting comparative genomics pipeline...")
        
        # Validate critical paths exist
        input_config = Path(f"{path_config.DATA_INPUT_DIR}/genes_to_proteins.json")
        if not input_config.exists():
            logger.error(f"Configuration file not found: {input_config}")
            return 1
        
        exit_code = asyncio.run(async_main())
        
        if exit_code == 0:
            logger.info("Pipeline completed successfully!")
        else:
            logger.error(f"Pipeline failed with exit code {exit_code}")
        
        return exit_code
        
    except KeyboardInterrupt:
        logger.info("Pipeline interrupted by user")
        return 1
    except Exception as e:
        logger.error(f"Fatal error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
