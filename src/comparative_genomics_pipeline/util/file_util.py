import json
import logging
from pathlib import Path
from typing import Optional, Dict, Any, Union
from ..config import path_config

logger = logging.getLogger(__name__)


def open_file_return_as_json(file: Union[str, Path]) -> Optional[Dict[str, Any]]:
    """Open and parse a JSON file with comprehensive error handling.
    
    Args:
        file: Path to the JSON file (string or Path object)
        
    Returns:
        Dict containing JSON data if successful, None if failed
    """
    try:
        file_path = Path(file)
        
        if not file_path.exists():
            logger.error(f"JSON file not found: {file_path}")
            return None
        
        if not file_path.is_file():
            logger.error(f"Path is not a file: {file_path}")
            return None
        
        if file_path.stat().st_size == 0:
            logger.error(f"JSON file is empty: {file_path}")
            return None
        
        with open(file_path, "r", encoding="utf-8") as f:
            try:
                contents = json.load(f)
                logger.debug(f"Successfully loaded JSON from {file_path}")
                return contents
            except json.JSONDecodeError as e:
                logger.error(f"Invalid JSON in file {file_path}: {e}")
                return None
            except UnicodeDecodeError as e:
                logger.error(f"Encoding error reading file {file_path}: {e}")
                return None
                
    except PermissionError:
        logger.error(f"Permission denied accessing file: {file}")
        return None
    except OSError as e:
        logger.error(f"OS error accessing file {file}: {e}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error loading JSON from {file}: {e}")
        return None


def open_file_return_as_str(file: Union[str, Path]) -> Optional[str]:
    """Open and read a file as string with comprehensive error handling.
    
    Args:
        file: Path to the file (string or Path object)
        
    Returns:
        File contents as string if successful, None if failed
    """
    try:
        file_path = Path(file)
        
        if not file_path.exists():
            logger.error(f"File not found: {file_path}")
            return None
        
        if not file_path.is_file():
            logger.error(f"Path is not a file: {file_path}")
            return None
        
        with open(file_path, "r", encoding="utf-8") as f:
            try:
                contents = f.read()
                logger.debug(f"Successfully read {len(contents)} characters from {file_path}")
                return contents
            except UnicodeDecodeError as e:
                logger.error(f"Encoding error reading file {file_path}: {e}")
                return None
                
    except PermissionError:
        logger.error(f"Permission denied accessing file: {file}")
        return None
    except OSError as e:
        logger.error(f"OS error accessing file {file}: {e}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error reading file {file}: {e}")
        return None


def validate_genes_config(config_data: Dict[str, Any]) -> bool:
    """Validate the structure of genes_to_proteins.json configuration.
    
    Args:
        config_data: Parsed JSON configuration data
        
    Returns:
        True if valid, False if validation errors found
    """
    if not isinstance(config_data, dict):
        logger.error("Configuration must be a JSON object")
        return False
    
    if not config_data:
        logger.error("Configuration is empty")
        return False
    
    valid = True
    
    for gene_name, ortholog_list in config_data.items():
        if not isinstance(gene_name, str) or not gene_name.strip():
            logger.error(f"Invalid gene name: {gene_name}")
            valid = False
            continue
        
        if not isinstance(ortholog_list, list):
            logger.error(f"Gene {gene_name}: ortholog list must be an array")
            valid = False
            continue
        
        if not ortholog_list:
            logger.error(f"Gene {gene_name}: ortholog list is empty")
            valid = False
            continue
        
        for i, ortholog in enumerate(ortholog_list):
            if not isinstance(ortholog, dict):
                logger.error(f"Gene {gene_name}, ortholog {i+1}: must be an object")
                valid = False
                continue
            
            # Check required fields
            if "species" not in ortholog:
                logger.error(f"Gene {gene_name}, ortholog {i+1}: missing 'species' field")
                valid = False
            elif not isinstance(ortholog["species"], str) or not ortholog["species"].strip():
                logger.error(f"Gene {gene_name}, ortholog {i+1}: 'species' must be a non-empty string")
                valid = False
            
            # Check that at least one ID is provided
            has_uniprot = ortholog.get("uniprot_id", "").strip()
            has_entrez = ortholog.get("entrez_protein_id", "").strip()
            
            if not has_uniprot and not has_entrez:
                logger.error(f"Gene {gene_name}, ortholog {i+1}: must have either 'uniprot_id' or 'entrez_protein_id'")
                valid = False
    
    if valid:
        logger.info(f"Configuration validation passed: {len(config_data)} genes configured")
    else:
        logger.error("Configuration validation failed")
    
    return valid


def safe_write_file(file_path: Union[str, Path], content: str, backup: bool = True) -> bool:
    """Safely write content to a file with backup and error handling.
    
    Args:
        file_path: Path to write to
        content: Content to write
        backup: Whether to create a backup of existing file
        
    Returns:
        True if successful, False if failed
    """
    try:
        file_path = Path(file_path)
        
        # Create parent directories if needed
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Create backup if file exists and backup requested
        if backup and file_path.exists():
            backup_path = file_path.with_suffix(f"{file_path.suffix}.backup")
            try:
                import shutil
                shutil.copy2(file_path, backup_path)
                logger.debug(f"Created backup: {backup_path}")
            except Exception as e:
                logger.warning(f"Failed to create backup of {file_path}: {e}")
        
        # Write content
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(content)
        
        logger.debug(f"Successfully wrote {len(content)} characters to {file_path}")
        return True
        
    except PermissionError:
        logger.error(f"Permission denied writing to file: {file_path}")
        return False
    except OSError as e:
        logger.error(f"OS error writing to file {file_path}: {e}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error writing to file {file_path}: {e}")
        return False


