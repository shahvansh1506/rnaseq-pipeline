"""
geo_fetcher.py — Downloads RNA-seq datasets from NCBI GEO.

What this file does:
- Connects to NCBI GEO database
- Downloads study metadata (sample names, conditions, organism)
- Downloads supplementary count matrix files
- Uses streaming download for large files (memory efficient)

Usage:
    fetcher = GEOFetcher("GSE60450")
    fetcher.download_metadata()
    fetcher.download_supplementary()
"""

import requests
import GEOparse
import logging
from pathlib import Path
from tqdm import tqdm
from src.config import RAW_DIR, LOG_FORMAT, LOG_LEVEL

# ── Setup logger for this module ──────────────────────────────
logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT)
logger = logging.getLogger(__name__)


class GEOFetcher:
    """
    Downloads GEO datasets and their associated files.

    GEO (Gene Expression Omnibus) is NCBI's public database
    for storing gene expression data. Each study has a GSE
    accession number e.g. GSE60450.

    Each study contains multiple GSM samples e.g. GSM1480291
    (one GSM per biological sample/replicate).
    """

    def __init__(self, geo_accession: str):
        """
        Args:
            geo_accession: GEO study ID e.g. "GSE60450"
        """
        self.accession = geo_accession
        self.gse       = None        # Will hold the full GEO study object
        self.metadata  = {}          # Will hold per-sample metadata

    # ──────────────────────────────────────────────────────────
    # PUBLIC METHODS
    # ──────────────────────────────────────────────────────────

    def download_metadata(self) -> dict:
        """
        Fetch and parse metadata for all samples in a GEO study.

        Metadata includes:
        - Sample title
        - Biological condition / characteristics
        - Organism
        - Sequencing instrument

        Returns:
            dict: {sample_name -> {title, characteristics, organism, instrument}}
        """
        logger.info(f"Fetching metadata for {self.accession}...")

        # GEOparse downloads the soft file and parses it automatically
        self.gse = GEOparse.get_GEO(
            geo     = self.accession,
            destdir = str(RAW_DIR),   # Save raw files here
            silent  = True,           # Suppress GEOparse's own output
        )

        # Loop through each sample (GSM) in the study
        for gsm_name, gsm in self.gse.gsms.items():
            self.metadata[gsm_name] = {
                "title"          : gsm.metadata.get("title", [""])[0],
                "characteristics": gsm.metadata.get("characteristics_ch1", []),
                "organism"       : gsm.metadata.get("organism_ch1", [""])[0],
                "instrument"     : gsm.metadata.get("instrument_model", [""])[0],
                "library_layout" : gsm.metadata.get("library_layout", [""])[0],
                "gsm_accession"  : gsm_name,
            }

        logger.info(f"Found {len(self.metadata)} samples in {self.accession}")
        return self.metadata

    def download_supplementary(self) -> list:
        """
        Download supplementary files attached to each sample.
        These are often pre-made count matrices (CSV/TXT files)
        which we can use directly without re-aligning reads.

        Returns:
            list: Paths to all downloaded files
        """
        if self.gse is None:
            raise RuntimeError("Call download_metadata() first!")

        downloaded = []

        for gsm_name, gsm in self.gse.gsms.items():
            # Each sample can have multiple supplementary files
            for key in ["supplementary_file_1", "supplementary_file_2"]:
                urls = gsm.metadata.get(key, [])

                for url in urls:
                    # Convert ftp:// to https:// for easier downloading
                    if url.startswith("ftp://"):
                        url = url.replace("ftp://", "https://")

                    # Build destination path
                    dest = RAW_DIR / f"{gsm_name}_{Path(url).name}"

                    # Skip if already downloaded
                    if dest.exists():
                        logger.info(f"Already exists, skipping: {dest.name}")
                        downloaded.append(dest)
                        continue

                    # Download the file
                    try:
                        self._stream_download(url, dest)
                        downloaded.append(dest)
                    except Exception as e:
                        logger.warning(f"Could not download {url}: {e}")

        logger.info(f"Downloaded {len(downloaded)} supplementary files.")
        return downloaded

    def print_summary(self) -> None:
        """Print a human-readable summary of the study metadata."""
        if not self.metadata:
            print("No metadata loaded. Run download_metadata() first.")
            return

        print(f"\n{'='*55}")
        print(f"  GEO Study: {self.accession}")
        print(f"  Total samples: {len(self.metadata)}")
        print(f"{'='*55}")

        for gsm, info in self.metadata.items():
            print(f"\n  Sample : {gsm}")
            print(f"  Title  : {info['title']}")
            print(f"  Organism: {info['organism']}")
            print(f"  Instrument: {info['instrument']}")
            print(f"  Characteristics:")
            for char in info["characteristics"]:
                print(f"    - {char}")

        print(f"\n{'='*55}\n")

    # ──────────────────────────────────────────────────────────
    # PRIVATE METHODS
    # ──────────────────────────────────────────────────────────

    def _stream_download(self, url: str, dest_path: Path) -> None:
        """
        Download a file in small chunks (8KB at a time).

        WHY streaming?
        Normal download = loads entire file into RAM first.
        Streaming download = loads 8KB at a time, writes to disk.
        This means a 10GB FASTQ file uses almost NO RAM!

        Args:
            url      : URL to download from
            dest_path: Where to save the file
        """
        logger.info(f"Downloading: {dest_path.name}")

        # stream=True tells requests NOT to load everything into memory
        response = requests.get(url, stream=True, timeout=60)
        response.raise_for_status()  # Raise error if download fails

        # Get total file size for the progress bar
        total_size = int(response.headers.get("content-length", 0))

        # Write file chunk by chunk with a progress bar
        with open(dest_path, "wb") as f, tqdm(
            total     = total_size,
            unit      = "B",
            unit_scale= True,
            desc      = dest_path.name,
        ) as progress_bar:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                progress_bar.update(len(chunk))

        logger.info(f"Saved: {dest_path}")
