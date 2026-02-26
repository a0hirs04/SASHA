from __future__ import annotations

import logging
import math
import re
import struct
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

try:
    from lxml import etree as XML  # type: ignore[import-not-found]
except ImportError:  # pragma: no cover
    import xml.etree.ElementTree as XML  # type: ignore[no-redef]

import numpy as np

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    pd = None


LOGGER = logging.getLogger(__name__)


@dataclass
class SimulationMetrics:
    total_tumor_cells: int
    total_stromal_cells: int
    live_tumor_cells: int
    activated_cafs: int
    mesenchymal_tumor_cells: int
    mean_ecm_density: float
    max_ecm_density: float
    mean_tumor_drug_sensitivity: float
    tumor_extent: float
    stroma_barrier_score: float
    drug_penetration: float
    hypoxic_fraction: float

    # Useful bookkeeping for downstream debugging / plotting.
    time: float = math.nan
    source_file: str = ""


class OutputParser:
    def __init__(self, output_dir: Path):
        candidate = Path(output_dir).expanduser().resolve()
        if not candidate.exists():
            raise FileNotFoundError(f"Output directory does not exist: {candidate}")

        # Accept either run_dir or run_dir/output.
        nested_output = candidate / "output"
        if nested_output.is_dir() and not list(candidate.glob("*.xml")):
            candidate = nested_output

        self.output_dir = candidate

    def parse_final_state(self) -> SimulationMetrics:
        final_xml = self._find_final_snapshot_xml()
        snapshot = self._read_physicell_xml(final_xml)
        return self._compute_metrics(snapshot, source_file=final_xml)

    def parse_timeseries(self):
        if pd is None:  # pragma: no cover
            raise ImportError("pandas is required for parse_timeseries()")

        xml_files = self._list_timeseries_xmls()
        if not xml_files:
            raise FileNotFoundError(f"No PhysiCell XML snapshots found in: {self.output_dir}")

        rows: List[Dict[str, Any]] = []
        for xml_file in xml_files:
            snapshot = self._read_physicell_xml(xml_file)
            metrics = self._compute_metrics(snapshot, source_file=xml_file)
            row = {
                "time": metrics.time,
                "tumor_count": metrics.total_tumor_cells,
                "stroma_count": metrics.total_stromal_cells,
                "live_tumor_cells": metrics.live_tumor_cells,
                "activated_cafs": metrics.activated_cafs,
                "mesenchymal_tumor_cells": metrics.mesenchymal_tumor_cells,
                "mean_ecm": metrics.mean_ecm_density,
                "max_ecm": metrics.max_ecm_density,
                "mean_tumor_drug_sensitivity": metrics.mean_tumor_drug_sensitivity,
                "tumor_extent": metrics.tumor_extent,
                "stroma_barrier_score": metrics.stroma_barrier_score,
                "drug_penetration": metrics.drug_penetration,
                "hypoxic_fraction": metrics.hypoxic_fraction,
                "source_file": metrics.source_file,
            }
            rows.append(row)

        df = pd.DataFrame(rows)
        if "time" in df.columns:
            df = df.sort_values("time").reset_index(drop=True)
        return df

    def _read_physicell_xml(self, filepath) -> dict:
        filepath = Path(filepath).expanduser().resolve()
        tree = XML.parse(str(filepath))
        root = tree.getroot()

        time_text = root.findtext(".//metadata/current_time")
        current_time = float(time_text) if time_text is not None else math.nan

        labels: Dict[int, Dict[str, Any]] = {}
        name_to_label: Dict[str, Dict[str, Any]] = {}
        for label_node in root.findall(".//cellular_information//simplified_data/labels/label"):
            index = self._safe_int(label_node.get("index"))
            size = self._safe_int(label_node.get("size"), default=1)
            name = (label_node.text or "").strip()
            if index is None or not name:
                continue
            entry = {"index": index, "size": max(1, size), "name": name}
            labels[index] = entry
            name_to_label[name] = entry

        cell_type_names: Dict[int, str] = {}
        for type_node in root.findall(".//cellular_information//cell_types/type"):
            type_value = self._safe_int(type_node.get("type"))
            if type_value is None:
                type_value = self._safe_int(type_node.get("ID"))
            type_name = (type_node.text or "").strip()
            if type_value is None or not type_name:
                continue
            cell_type_names[type_value] = type_name

        cell_data_filename = root.findtext(".//cellular_information//simplified_data/filename")
        if not cell_data_filename:
            raise ValueError(f"Cell data filename not found in {filepath}")
        cell_data_path = (filepath.parent / cell_data_filename).resolve()
        cell_matrices = self._read_mat_v4(cell_data_path)
        cell_matrix = self._pick_matrix(cell_matrices, preferred_names=("cells",))

        micro_data_filename = root.findtext(".//microenvironment//domain/data/filename")
        if not micro_data_filename:
            raise ValueError(f"Microenvironment data filename not found in {filepath}")
        micro_data_path = (filepath.parent / micro_data_filename).resolve()
        micro_matrices = self._read_mat_v4(micro_data_path)
        micro_matrix = self._pick_matrix(
            micro_matrices,
            preferred_names=("multiscale_microenvironment", "microenvironment"),
        )

        variables: List[Tuple[int, str]] = []
        for variable_node in root.findall(".//microenvironment/domain/variables/variable"):
            var_id = self._safe_int(variable_node.get("ID"))
            var_name = (variable_node.get("name") or "").strip()
            if var_id is None or not var_name:
                continue
            variables.append((var_id, var_name))
        variables.sort(key=lambda x: x[0])

        var_row_map = self._build_variable_row_map(variables, micro_matrix.shape[0])
        coords = micro_matrix[0:3, :].T if micro_matrix.shape[0] >= 3 else np.empty((0, 3), dtype=float)
        micro_values: Dict[str, np.ndarray] = {}
        for _, var_name in variables:
            row_idx = var_row_map.get(var_name)
            if row_idx is None or row_idx >= micro_matrix.shape[0]:
                continue
            micro_values[var_name] = micro_matrix[row_idx, :]

        return {
            "filepath": filepath,
            "time": current_time,
            "cell_matrix": cell_matrix,
            "label_map": labels,
            "label_name_map": name_to_label,
            "cell_type_names": cell_type_names,
            "micro_matrix": micro_matrix,
            "micro_coords": coords,
            "micro_values": micro_values,
        }

    def _compute_stroma_barrier_score(self, cell_positions, ecm_field):
        cell_positions = np.asarray(cell_positions, dtype=float)
        if cell_positions.size == 0:
            return math.nan

        coords = np.asarray(ecm_field.get("coords", []), dtype=float)
        values = np.asarray(ecm_field.get("values", []), dtype=float)
        if coords.size == 0 or values.size == 0:
            return math.nan

        centroid = np.nanmean(cell_positions, axis=0)
        radii = np.linalg.norm(coords - centroid, axis=1)
        annulus_mask = (radii >= 100.0) & (radii <= 300.0)
        if not np.any(annulus_mask):
            return math.nan

        return float(np.nanmean(values[annulus_mask]))

    def _compute_metrics(self, snapshot: Dict[str, Any], source_file: Path) -> SimulationMetrics:
        cell_matrix: np.ndarray = snapshot["cell_matrix"]
        label_name_map: Dict[str, Dict[str, Any]] = snapshot["label_name_map"]
        micro_values: Dict[str, np.ndarray] = snapshot["micro_values"]
        micro_coords: np.ndarray = snapshot["micro_coords"]

        cell_type = self._get_row(cell_matrix, label_name_map, "cell_type")
        if cell_type is None:
            cell_type_int = np.zeros(cell_matrix.shape[1], dtype=int)
        else:
            cell_type_int = np.rint(cell_type).astype(int)

        tumor_mask = cell_type_int == 0
        stroma_mask = cell_type_int == 1

        dead = self._get_row(cell_matrix, label_name_map, "dead")
        death_model = self._get_row(cell_matrix, label_name_map, "current_death_model")
        dead_mask = (dead > 0.5) if dead is not None else np.zeros_like(tumor_mask, dtype=bool)
        apoptotic_mask = (np.rint(death_model).astype(int) == 100) if death_model is not None else np.zeros_like(
            tumor_mask, dtype=bool
        )
        live_tumor_mask = tumor_mask & ~(dead_mask | apoptotic_mask)

        is_activated = self._get_row(cell_matrix, label_name_map, "is_activated")
        if is_activated is None:
            is_activated = self._get_row(cell_matrix, label_name_map, "ACTA2")
        activated_cafs = int(np.sum(stroma_mask & (is_activated > 0.5))) if is_activated is not None else 0

        is_mesenchymal = self._get_row(cell_matrix, label_name_map, "is_mesenchymal")
        mesenchymal_tumor_cells = int(np.sum(tumor_mask & (is_mesenchymal > 0.5))) if is_mesenchymal is not None else 0

        drug_sensitivity = self._get_row(cell_matrix, label_name_map, "drug_sensitivity")
        if drug_sensitivity is not None and np.any(tumor_mask):
            mean_tumor_drug_sensitivity = float(np.nanmean(drug_sensitivity[tumor_mask]))
        else:
            mean_tumor_drug_sensitivity = math.nan

        hif1a = self._get_row(cell_matrix, label_name_map, "HIF1A")
        if hif1a is not None and np.any(tumor_mask):
            hypoxic_fraction = float(np.mean(hif1a[tumor_mask] > 0.5))
        else:
            hypoxic_fraction = math.nan

        positions = self._get_positions(cell_matrix, label_name_map)
        tumor_positions = positions[tumor_mask, :] if positions.size else np.empty((0, 3), dtype=float)

        ecm = micro_values.get("ecm_density")
        drug = micro_values.get("drug")

        mean_ecm_density = float(np.nanmean(ecm)) if ecm is not None and ecm.size else math.nan
        max_ecm_density = float(np.nanmax(ecm)) if ecm is not None and ecm.size else math.nan

        tumor_extent = self._max_pairwise_distance(tumor_positions)
        stroma_barrier_score = self._compute_stroma_barrier_score(
            tumor_positions,
            {"coords": micro_coords, "values": ecm if ecm is not None else np.array([])},
        )
        drug_penetration = self._compute_drug_penetration(
            tumor_positions,
            micro_coords,
            drug if drug is not None else np.array([]),
        )

        metrics = SimulationMetrics(
            total_tumor_cells=int(np.sum(tumor_mask)),
            total_stromal_cells=int(np.sum(stroma_mask)),
            live_tumor_cells=int(np.sum(live_tumor_mask)),
            activated_cafs=activated_cafs,
            mesenchymal_tumor_cells=mesenchymal_tumor_cells,
            mean_ecm_density=mean_ecm_density,
            max_ecm_density=max_ecm_density,
            mean_tumor_drug_sensitivity=mean_tumor_drug_sensitivity,
            tumor_extent=tumor_extent,
            stroma_barrier_score=stroma_barrier_score,
            drug_penetration=drug_penetration,
            hypoxic_fraction=hypoxic_fraction,
            time=float(snapshot.get("time", math.nan)),
            source_file=str(source_file),
        )
        return metrics

    def _compute_drug_penetration(
        self,
        tumor_positions: np.ndarray,
        voxel_coords: np.ndarray,
        drug_values: np.ndarray,
    ) -> float:
        if tumor_positions.size == 0 or voxel_coords.size == 0 or drug_values.size == 0:
            return math.nan

        centroid = np.nanmean(tumor_positions, axis=0)
        tumor_radius = float(np.max(np.linalg.norm(tumor_positions - centroid, axis=1)))
        # Add a small shell so edge voxels are included even for sparse tumor layouts.
        region_radius = max(tumor_radius + 20.0, 20.0)

        voxel_dist = np.linalg.norm(voxel_coords - centroid, axis=1)
        mask = voxel_dist <= region_radius
        if not np.any(mask):
            return math.nan
        return float(np.nanmean(drug_values[mask]))

    def _max_pairwise_distance(self, points: np.ndarray) -> float:
        points = np.asarray(points, dtype=float)
        n = points.shape[0]
        if n < 2:
            return 0.0

        max_d2 = 0.0
        chunk = 256
        for i in range(0, n, chunk):
            p = points[i : i + chunk]  # [c, 3]
            diff = p[:, None, :] - points[None, :, :]  # [c, n, 3]
            d2 = np.einsum("ijk,ijk->ij", diff, diff, optimize=True)
            local_max = float(np.max(d2))
            if local_max > max_d2:
                max_d2 = local_max
        return float(math.sqrt(max_d2))

    def _get_positions(self, cell_matrix: np.ndarray, label_name_map: Dict[str, Dict[str, Any]]) -> np.ndarray:
        label = label_name_map.get("position")
        if label is None:
            return np.empty((0, 3), dtype=float)

        start = int(label["index"])
        size = int(label["size"])
        if size < 3:
            return np.empty((0, 3), dtype=float)
        stop = min(start + 3, cell_matrix.shape[0])
        if stop - start < 3:
            return np.empty((0, 3), dtype=float)
        return cell_matrix[start:stop, :].T

    def _get_row(
        self,
        cell_matrix: np.ndarray,
        label_name_map: Dict[str, Dict[str, Any]],
        label_name: str,
    ) -> Optional[np.ndarray]:
        label = label_name_map.get(label_name)
        if label is None:
            return None
        idx = int(label["index"])
        if idx < 0 or idx >= cell_matrix.shape[0]:
            return None
        return cell_matrix[idx, :]

    def _build_variable_row_map(self, variables: Sequence[Tuple[int, str]], micro_rows: int) -> Dict[str, int]:
        if not variables:
            return {}
        var_ids = [vid for vid, _ in variables]
        # In PhysiCell MAT exports this is typically: x,y,z,voxel_volume + variables.
        base = micro_rows - len(variables)
        row_map: Dict[str, int] = {}
        contiguous_ids = (min(var_ids) == 0) and (sorted(var_ids) == list(range(max(var_ids) + 1)))
        if contiguous_ids:
            for var_id, var_name in variables:
                row_map[var_name] = base + var_id
        else:
            for i, (_, var_name) in enumerate(sorted(variables, key=lambda x: x[0])):
                row_map[var_name] = base + i
        return row_map

    def _read_mat_v4(self, filepath: Path) -> Dict[str, np.ndarray]:
        if not filepath.exists():
            raise FileNotFoundError(f"MAT file not found: {filepath}")

        buf = filepath.read_bytes()
        offset = 0
        matrices: Dict[str, np.ndarray] = {}

        while offset + 20 <= len(buf):
            type_code, mrows, ncols, imagf, namelen = struct.unpack_from("<iiiii", buf, offset)
            offset += 20

            if namelen <= 0 or offset + namelen > len(buf):
                break

            name_bytes = buf[offset : offset + namelen]
            offset += namelen
            name = name_bytes.split(b"\x00", 1)[0].decode("utf-8", errors="ignore")
            if not name:
                name = f"matrix_{len(matrices)}"

            if imagf != 0:
                raise ValueError(f"Complex MAT v4 arrays are not supported: {filepath}")
            if type_code != 0:
                raise ValueError(f"Unsupported MAT v4 type code {type_code} in {filepath}")

            count = mrows * ncols
            bytes_needed = count * 8
            if count < 0 or offset + bytes_needed > len(buf):
                raise ValueError(f"Corrupt MAT v4 payload in {filepath}")

            matrix = np.frombuffer(buf, dtype="<f8", count=count, offset=offset)
            matrix = np.array(matrix, copy=True).reshape((mrows, ncols), order="F")
            matrices[name] = matrix
            offset += bytes_needed

        if not matrices:
            raise ValueError(f"No matrices parsed from {filepath}")
        return matrices

    def _pick_matrix(self, matrices: Dict[str, np.ndarray], preferred_names: Sequence[str]) -> np.ndarray:
        for name in preferred_names:
            if name in matrices:
                return matrices[name]
        return next(iter(matrices.values()))

    def _find_final_snapshot_xml(self) -> Path:
        numbered = self._numbered_output_xmls()
        if numbered:
            return numbered[-1]

        final_xml = self.output_dir / "final.xml"
        if final_xml.exists():
            return final_xml

        xml_candidates = sorted(self.output_dir.glob("*.xml"))
        if not xml_candidates:
            raise FileNotFoundError(f"No XML snapshots found in {self.output_dir}")
        return xml_candidates[-1]

    def _list_timeseries_xmls(self) -> List[Path]:
        numbered = self._numbered_output_xmls()
        if numbered:
            # Include final.xml if present for end-of-run diagnostics.
            final_xml = self.output_dir / "final.xml"
            if final_xml.exists():
                return numbered + [final_xml]
            return numbered

        fallback = []
        for name in ("initial.xml", "final.xml"):
            path = self.output_dir / name
            if path.exists():
                fallback.append(path)
        if fallback:
            return fallback
        return sorted(self.output_dir.glob("*.xml"))

    def _numbered_output_xmls(self) -> List[Path]:
        pattern = re.compile(r"^output(\d+)\.xml$")
        numbered: List[Tuple[int, Path]] = []
        for path in self.output_dir.glob("output*.xml"):
            match = pattern.match(path.name)
            if not match:
                continue
            numbered.append((int(match.group(1)), path))
        numbered.sort(key=lambda x: x[0])
        return [p for _, p in numbered]

    def _safe_int(self, value: Optional[str], default: Optional[int] = None) -> Optional[int]:
        if value is None:
            return default
        try:
            return int(value)
        except (TypeError, ValueError):
            return default


def metrics_to_dict(metrics: SimulationMetrics) -> Dict[str, Any]:
    return asdict(metrics)
