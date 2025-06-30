import httpx
import asyncio
import logging

logger = logging.getLogger(__name__)

class ClinVarClient:
    """Minimal ClinVar client using NCBI E-utilities JSON API"""
    
    def __init__(self):
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        
    async def search_variants_by_gene(self, gene_name: str, limit: int = 20):
        """Search for variants by gene name"""
        # Add delay to avoid rate limiting
        await asyncio.sleep(0.3)
        
        search_url = f"{self.base_url}/esearch.fcgi"
        params = {
            "db": "clinvar",
            "term": f"{gene_name}[gene]",
            "retmode": "json",
            "retmax": limit
        }
        
        async with httpx.AsyncClient() as client:
            response = await client.get(search_url, params=params)
            response.raise_for_status()
            data = response.json()
            
            if "esearchresult" in data and "idlist" in data["esearchresult"]:
                return data["esearchresult"]["idlist"]
            return []
    
    async def get_variant_summaries(self, variant_ids: list):
        """Get summary data for variant IDs"""
        if not variant_ids:
            return []
        
        # Add delay to avoid rate limiting
        await asyncio.sleep(0.5)
            
        summary_url = f"{self.base_url}/esummary.fcgi"
        params = {
            "db": "clinvar",
            "id": ",".join(variant_ids),
            "retmode": "json"
        }
        
        async with httpx.AsyncClient() as client:
            response = await client.get(summary_url, params=params)
            response.raise_for_status()
            data = response.json()
            
            variants = []
            if "result" in data:
                for vid in variant_ids:
                    if vid in data["result"]:
                        variant_data = data["result"][vid]
                        
                        # Extract clinical significance from nested structure
                        clinical_sig = ""
                        if "germline_classification" in variant_data:
                            clinical_sig = variant_data["germline_classification"].get("description", "")
                        
                        # Extract gene symbol from genes array
                        gene_symbol = ""
                        if "genes" in variant_data and variant_data["genes"]:
                            gene_symbol = variant_data["genes"][0].get("symbol", "")
                        
                        variants.append({
                            "id": vid,
                            "title": variant_data.get("title", ""),
                            "clinical_significance": clinical_sig,
                            "gene_symbol": gene_symbol,
                            "variation_type": variant_data.get("obj_type", ""),
                            "protein_change": variant_data.get("protein_change", "")
                        })
            return variants
    
    async def get_gene_variants(self, gene_name: str, limit: int = 20):
        """Get variants for a gene with basic info"""
        logger.info(f"Fetching ClinVar variants for {gene_name}")
        
        # Search for variants
        variant_ids = await self.search_variants_by_gene(gene_name, limit)
        if not variant_ids:
            logger.warning(f"No variants found for {gene_name}")
            return []
        
        # Get summary data
        variants = await self.get_variant_summaries(variant_ids)
        logger.info(f"Found {len(variants)} variants for {gene_name}")
        
        return variants