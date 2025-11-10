"""
Utilities for converting Roofer CityJSONL exports into Tyler-compatible assets.

From a CityJSONL file produced by Roofer using default settings for --[no-]split-cjseq (false) and --[no]-omit-metadata (false), this script extracts the metadata record stored in the first CityObject in the CityObjects dictionary, writes it to a `metadata.city.json` file, and writes individual feature records in separate JSONL files an output directory at the location of the metadata file.

Example (PowerShell):
    python .\roofer2tyler.py "D:\tyler\data-jasper\Jasper.city.jsonl"
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, Tuple


def iter_cityjsonl(jsonl_path: Path) -> Iterator[Tuple[Dict[str, Any], str]]:
    """
    Yield parsed JSON objects and their raw lines from a CityJSONL file.

    Args:
        jsonl_path (Path): Source CityJSONL file path.

    Yields:
        Iterator[Tuple[Dict[str, Any], str]]: Parsed JSON payload and the raw line.
    """

    with jsonl_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            yield json.loads(stripped), stripped


REQUIRED_METADATA_KEYS = {"type", "version", "transform", "metadata", "CityObjects", "vertices"}
REQUIRED_FEATURE_KEYS = {"CityObjects", "id", "type", "vertices"}


def validate_metadata(payload: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate the metadata payload using lightweight checks.

    Args:
        payload (Dict[str, Any]): Raw metadata record.

    Returns:
        Dict[str, Any]: Validated metadata mapping.

    Raises:
        ValueError: If required keys are missing or have incompatible types.
    """

    missing = REQUIRED_METADATA_KEYS - payload.keys()
    if missing:
        raise ValueError(f"Metadata record is missing required keys: {sorted(missing)}")

    if not isinstance(payload["CityObjects"], dict):
        raise ValueError("Metadata 'CityObjects' must be a dictionary.")
    if not isinstance(payload["transform"], dict):
        raise ValueError("Metadata 'transform' must be a dictionary.")
    if not isinstance(payload["metadata"], dict):
        raise ValueError("Metadata 'metadata' must be a dictionary.")
    if not isinstance(payload["vertices"], list):
        raise ValueError("Metadata 'vertices' must be a list.")

    return payload


def validate_feature(payload: Dict[str, Any]) -> Dict[str, Any]:
    """
    Validate a feature payload using lightweight checks.

    Args:
        payload (Dict[str, Any]): Raw feature record.

    Returns:
        Dict[str, Any]: Validated feature mapping.

    Raises:
        ValueError: If required keys are missing or have incompatible types.
    """

    missing = REQUIRED_FEATURE_KEYS - payload.keys()
    if missing:
        raise ValueError(f"Feature record is missing required keys: {sorted(missing)}")

    if not isinstance(payload["CityObjects"], dict):
        raise ValueError("Feature 'CityObjects' must be a dictionary.")
    if not isinstance(payload["vertices"], list):
        raise ValueError("Feature 'vertices' must be a list.")
    if not isinstance(payload["id"], str) or not payload["id"]:
        raise ValueError("Feature 'id' must be a non-empty string.")

    return payload


def write_metadata(metadata: Dict[str, Any], destination: Path) -> None:
    """
    Persist the metadata record to disk using Tyler's expected structure.

    Args:
        metadata (Dict[str, Any]): Validated metadata record.
        destination (Path): Output file path.
    """

    destination.write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def write_feature(feature: Dict[str, Any], raw_line: str, features_dir: Path) -> None:
    """
    Write a single feature record to a JSONL file.

    Args:
        feature (Dict[str, Any]): Validated feature record.
        raw_line (str): Original JSON line to persist.
        features_dir (Path): Target directory for feature files.
    """

    feature_id = feature.get("id")
    if not feature_id:
        raise ValueError("Feature record is missing an 'id' key.")

    feature_path = features_dir / f"{feature_id}.jsonl"
    feature_path.write_text(f"{raw_line}\n", encoding="utf-8")


def process_cityjsonl(jsonl_path: Path) -> Tuple[Path, Path]:
    """
    Convert a Roofer CityJSONL file into Tyler-compatible metadata and feature assets.

    Args:
        jsonl_path (Path): Source CityJSONL file path.

    Returns:
        Tuple[Path, Path]: Paths to the generated metadata file and features directory.
    """

    if not jsonl_path.exists():
        raise FileNotFoundError(f"Input file {jsonl_path} does not exist.")

    record_iter = iter_cityjsonl(jsonl_path)
    try:
        metadata_payload, _ = next(record_iter)
    except StopIteration as exc:
        raise ValueError(f"{jsonl_path} is empty; expected metadata line.") from exc

    metadata = validate_metadata(metadata_payload)

    base_dir = jsonl_path.parent
    metadata_path = base_dir / "metadata.city.json"
    features_dir = base_dir / "features"

    if features_dir.exists():
        shutil.rmtree(features_dir)
    features_dir.mkdir(parents=True, exist_ok=True)

    feature_count = 0
    for feature_payload, raw_line in record_iter:
        feature = validate_feature(feature_payload)
        write_feature(feature, raw_line, features_dir)
        feature_count += 1

    if feature_count == 0:
        raise ValueError(f"{jsonl_path} does not contain any CityObject feature records.")

    write_metadata(metadata, metadata_path)
    return metadata_path, features_dir


def parse_args(args: Iterable[str] | None = None) -> argparse.Namespace:
    """
    Parse CLI arguments.

    Args:
        args (Iterable[str] | None): Optional arguments for testing.

    Returns:
        argparse.Namespace: Namespace containing parsed arguments.
    """

    parser = argparse.ArgumentParser(description="Convert Roofer CityJSONL exports into Tyler assets.")
    # The source JSONL must contain several CityObject items to satisfy Tyler's feature batching model.
    parser.add_argument("jsonl_file", type=Path, help="Path to the Roofer CityJSONL file.")
    return parser.parse_args(args=args)


def main() -> None:
    """
    Entry point for the CLI script.
    """

    namespace = parse_args()
    try:
        metadata_path, features_dir = process_cityjsonl(namespace.jsonl_file)
    except (json.JSONDecodeError, ValueError, FileNotFoundError) as exc:
        raise SystemExit(str(exc))

    print(f"metadata: {metadata_path}")
    print(f"features: {features_dir}")


if __name__ == "__main__":
    main()

