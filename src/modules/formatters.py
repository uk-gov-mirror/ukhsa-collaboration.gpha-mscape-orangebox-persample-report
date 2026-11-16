import pandas as pd
import numpy as np
import logging

logger = logging.getLogger(__name__)

# add in the metadata fields
METADATA_FIELDS = {
    "climb_id": "CLIMB ID",
    "site": "Site", 
    "published_date": "Published date", 
    "collection_date": "Collection date", 
    "sample_source": "Sample source", 
    "sample_type": "Sample type"
}

HCID_COUNTS_COLUMNS = {"name", "taxon_id"}

class FormattedSection:
    def __init__(self, status: str, table: pd.DataFrame | None = None, message: str = ""):
        self.status = status
        if table is None:
            table = pd.DataFrame()  # a new empty table for each section
        self.table = table
        self.message = message

def get_single_row(data: pd.DataFrame | None, section: str) -> pd.Series | None:
    """
    Return the only row of a per-sample DataFrame, None if there is no data, or raise if there are several lines.
    Enforces the assumption that each section of the report is per-sample, so there should be only one row of data.
    """
    if data is None or data.empty:
        return None
    if len(data) != 1:
        raise ValueError(f"{section}: expected data for one sample, got {len(data)} rows")
    return data.iloc[0]

def is_blank(value) -> bool:
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip() == ""
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):  # lists and other containers
        return False

def format_date(value) -> str:
    """Format a date as DD/MM/YYYY. Raises ValueError if it cannot be parsed."""
    return pd.to_datetime(value, errors="raise").strftime("%d/%m/%Y")

def format_hcid(counts: pd.DataFrame | None, warnings:list[dict]) -> FormattedSection:
    """
    Input func deets here pls

    1. If any warning.json present -> flagged, list taxa and report. 
    2. Else hcid.counts.csv found -> screened and not flagged
    3. Else (nothing found) - Error or HCID check not run, reusult not available. 

    Args: 
        counts: the hcid.counts.csv table, or none if not found 
        warnings: the contents of each warning json
    """
    if counts is not None and not counts.empty:
        missing_columns = HCID_COUNTS_COLUMNS - set(counts.columns)
        if missing_columns:
            raise ValueError(f"Missing expected columns in hcid.counts.csv: {missing_columns}")

    # 1. Flaggeed
    if warnings:
        return format_hcid_flagged(counts, warnings)

    # 2. Screened and not flagged - phew
    if counts is not None and not counts.empty:
        logger.info("HCID check completed: sample screened and not flagged.")
        return format_hcid_screened(counts)

    # 3. Not run or an error raised
    logger.warning("HCID check not run or no results available.")
    return FormattedSection(status=STATUS_MISSING, message="HCID screening result not available for this sample.")

def format_hcid_flagged(counts: pd.DataFrame | None, warnings: list[dict]) -> FormattedSection:
    """
    Build the flagged message and table for the HCID section. 
    Look up taxon counts in the counts table if available, otherwise just report the taxon_id.
    """
    names = {} # taxon_id -> HCID name; such as 33278 = AI-A infection H7N9. 

    # find taxon and name for each warning and add to the table. 
    if counts is not None and not counts.empty: 
        names = dict(zip(counts["taxon_id"], counts["name"]))
    else:
        logger.error("HCID counts table not available, will report taxon_id only for flagged taxa.")

    rows = []
    for warning in warnings:
        if "taxon_id" not in warning:
            raise ValueError(f"Warning JSON missing 'taxon_id': {warning}")
        taxon_id = warning["taxon_id"]
        name = names.get(taxon_id)

        if name is None:
            # register a positive result, we can't find the parsewd info
            logger.error(f"HCID counts table does not contain taxon_id {taxon_id} from warning JSON.")
            name = f"Unknown taxon_id {taxon_id}"
        rows.append({"HCID taxon": taxon_id, "name": name})

    table = pd.DataFrame(rows)
    message = f"HCID check completed: sample flagged for {len(warnings)} taxon(s)."
    return FormattedSection(status=STATUS_OK, table=table, message=message)

def format_metadata(data: pd.DataFrame) -> FormattedSection:
        """
        Function to format the metadata of the run. Included are the following fields:
        •	Climb-id
        •	Site
        •	Published date
        •	Collection date
        •	Sample source
        •	Sample type


        Args: 
            data (pd.DataFrame): dataframe containing the flattened json metadata.

        Returns: 
            formatted section containing all of the sample metadata. 
        """
    if data is not None and not data.empty:
        row = get_single_row(data, "metadata")
    else:
        logger.warning("metadata is not present in the sample output")
        return FormattedSection(status = STATUS_MISSING, message="Sample metadata not found for this sample")

    rows = []
    # retrieve the metadata fields
    for field_name, label in METADATA_FIELDS.items():
        value = row.get(field_name)
        if is_blank(value):
            # record blank in logger
            logger.info(f"Metadata field {field_name} is missing or blank.")
            display = NOT_RECORDED
        elif field_name in METADATA_DATE_FIELDS:
            display = format_date(value)
        else:
            display = str(value)

    return FormattedSection(status=STATUS_OK, table=pd.DataFrame(rows, columns=["Field", "Value"]))
    
def format_kraken_genera(data: pd.DataFrame) -> FormattedSection:
    """
    Format the Kraken genera output for HTML and plotly output.

    Args:
        data (pd.DataFrame): pandas dataframe of kraken outputs
    """
    if data is None or data.empty:
        logger.warning(f"Kraken genus results not available")
        return FormattedSection(
            status = STATUS_MISSING, 
            message = "Kraken genus results not available for this sample."
        )
    
    missing_columns = KRAKEN_GENRA_COLUMNS - set(data.columns)
    if missing_columns:
        raise ValueError(f"Expected Kraken columns are missing")

    genera = data.














